function [H, R, active_UEs] = function_Setup(square_length,L, Kdrop, B, noiseVardBm, Kmax, f, M, polarization,carrier_freq, num_subcarriers)

%% Set Network Layout

% Number of channel realizations
%num_subcarriers = 400;

% Generate Simulation parameters
s = qd_simulation_parameters;
s.sample_density = 1;
s.use_absolute_delays = 1;
s.show_progress_bars = 0;

center_frequency = carrier_freq;
wavelength = 3e8/center_frequency;
s.center_frequency = center_frequency;

% New Quadriga Layout
l = qd_layout(s);

% Generate Base Stations
l.no_tx = L;

% Base stations heights
l.tx_position(3,:) = 25;

%Set the length in meters of the total square area
squareLength = square_length;

%Number of BSs per dimension
BSsPerDim = sqrt(L);

%% Base Station Deployment

% Minimum distance between BSs and UEs
min_distance = 35;

%Distance between BSs in vertical/horizontal direction
interSiteDistance = squareLength/BSsPerDim;

% Deploy BSs on the grid
horizontal_axis = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[BSsPerDim 1]);
vertical_axis = horizontal_axis';
BS_positions = horizontal_axis(:) + 1i*vertical_axis(:);

% Give BS postions to Quadriga layout to set the positions
for i = 1:length(BS_positions)
    
    l.tx_position(1:2,i) = [real(BS_positions(i)) ; imag(BS_positions(i))];
        
end

%% Creating Antenna array (Planner) for each BS

%Number of vertical antennas
Mv = 10;

%Number of antennas in horizontal direction
Mh = M/Mv;

%Antenna Spacing 
ant_spacing = 1/2;

%Polarization (dual polarized)
polarization_Indicator = polarization;
    
%Go through all BSs

for i = 1:L
    
    l.tx_array(i).generate('3gpp-3d', 1,Mv, Mh, center_frequency, polarization_Indicator, 0, ant_spacing);
     
    for j = 1:Mv
        for k = 1:Mh
            
            index = (j-1)*Mh + k;
            l.tx_array(i).element_position(2,index) = mod(index-1,Mh)*ant_spacing*wavelength;
            l.tx_array(i).element_position(3,index) = floor(index-1/Mh)*ant_spacing*wavelength;
            
        end
    end
    
end
%% Drop UEs in the cells

user_positions = zeros(Kdrop,L);
users_perBS = zeros(L,1);

% Maximum cell distance
maxCell_distance = interSiteDistance; 

% Go through all cells
for i=1:L
    
    % Put users randomly in the cell
    while users_perBS(i)<Kdrop
        
        %Drop users
        UEremaining = Kdrop-users_perBS(i);
        X_coord = rand(UEremaining,1)*maxCell_distance - maxCell_distance/2;
        Y_coord = rand(UEremaining,1)*maxCell_distance - maxCell_distance/2;
        XY_coord = X_coord + 1i*Y_coord;
        
        %Keep those UEs that are at minimum distance (35 meter) away from
        %BS
        XY_coord = XY_coord(abs(XY_coord)>=min_distance);
        
        %Store new UEs
        user_positions(users_perBS(i)+1:users_perBS(i)+length(XY_coord),i) = XY_coord + BS_positions(i);
        users_perBS(i) = users_perBS(i)+length(XY_coord);
        
    end
end

%% Receiver Characteristics

% Total number of UEs
UE_total = Kdrop*L;

% UE height (meters)
UE_height = 1.5;

% Generate UEs in the layout
l.no_rx = UE_total;

% UEs antenna type
l.rx_array.generate('omni');

% Extract UE positions from user_positions matrix
reshape_user_positions = reshape(user_positions,[UE_total 1]);

% Get all users positions
for i = 1:length(reshape_user_positions)
    
    l.rx_position(1:3,i) = [real(reshape_user_positions(i)) ; imag(reshape_user_positions(i)) ; UE_height];
    
end

% UE tracks

for j=1:UE_total
    
    l.track(j).generate('linear',2,0); %Define a linear track of 2m and east direction
    l.track(j).scenario = '3GPP_38.901_UMi_NLOS'; %Select the Urban Microcell NLOS scenario
    l.track(j).set_speed(135/3.6)    % Mobility speed 135 km/h
    
end

%Generate pilot patterns
if f == 1
    
    pilotPattern = ones(L,1);
    
elseif f == 2 %Only works for 16 BSs
    
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]);
    
end


%Randomize pilot allocation in each cell
randOrder = zeros(Kmax*f,L);

for j = 1:L
    
    randOrder(1+(pilotPattern(j)-1)*Kmax:pilotPattern(j)*Kmax,j) = randperm(Kmax)+(pilotPattern(j)-1)*Kmax;
    
end

% Noise Variance and standard deviation
noiseVar_dBm = noiseVardBm;  %(dBm)
noise_Var = 10^(noiseVar_dBm/10);
noise_Std = sqrt(noise_Var);
%% Generate channel output and Spatial correlation Matrix

% Intialize variables to store values
H = zeros(M,num_subcarriers,Kmax,L,L);
R = zeros(M,M,Kmax,L,L);
users_perBS = zeros(L,1);
active_UEs = zeros(Kmax,L);



for i = 1:L
       
%Output simulation progress
disp([num2str(i) ' cells generated out of ' num2str(L)]);
    
% Get Channels
ch = l.get_channels;
    
    for j = 1:Kdrop
        
        h_k = zeros(M,num_subcarriers,1,1,L);
        r_k = zeros(M,M,1,1,L);
        
        for k = 1:L
            
            [~,ind] =  min(abs(user_positions(j,:)-BS_positions(k)));
            
            % Find the right user index
            user_ind = find(user_positions(j,ind)==reshape_user_positions);
            
            % Select channel for individual user
            h = ch(user_ind,k);
            
            H_freq = h.fr(B, num_subcarriers);  % Getting frequency response
            
            % Getting Channels of individual user from all BSs
            H_freq = mean(H_freq,4);
            h_k(:,:,1,1,k) = reshape(H_freq,[M num_subcarriers])/noise_Std;
            
            % Getting Spatial Correlation Matrix of individual user from all BSs
            r_k(:,:,1,1,k) = diag(mean(abs(h_k(:,:,1,1,k)).^2,2)/noise_Var);
            
        end
         %Determine which BS should serve the UE
        [max_value,BS_index] = max(mean(sum(abs(h_k(:,:,1,1,:)).^2,1),2));
        
        if users_perBS(BS_index)<Kmax
           
           users_perBS(BS_index) = users_perBS(BS_index) + 1; 
           H(:,:,randOrder(users_perBS(BS_index)+(pilotPattern(BS_index)-1)*Kmax,BS_index),BS_index,:) = h_k;
           R(:,:,randOrder(users_perBS(BS_index)+(pilotPattern(BS_index)-1)*Kmax,BS_index),BS_index,:) = r_k;
           active_UEs(randOrder(users_perBS(BS_index)+(pilotPattern(BS_index)-1)*Kmax,BS_index),BS_index) = 1; 
        end
        
    end
end
end