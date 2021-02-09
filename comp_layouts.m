% function [H, Rest] = setup_quadriga(L, Kdrop, B, noiseVardBm, Kmax, f, M, polarization)

%% Set simulation parameters
close all
clc
% Number of users to be droped in cell 
Kdrop = 5;

% Maximum users served by BS
Kmax = 6;

% Number of arrays on BS antenna
M = 100;

% Number of Cells
L = 16;

% Uplink Power
up_pwr = 100;

% Number of channel realizations
num_subcarriers = 400;

% Number of setups with random UE locations
num_setups = 2;

% Transmitter power
txPwr = 1000;      %(mW)

% Pilot reuse factor
f = 2;

% Communication Bandwidth (20 MHz)
B = 20e6;

% Effective bandwidth after 5% overhead
effective_BW = B*0.95;

% Length of pilot Signals
tau_p = Kmax*f;

% Length of coherence block
tau_c = 190;

%% Calculate DL and UL throughput of UEs for number of setups

% Initialize variable to store simulations
SE_MR_UL = zeros(tau_p,L,num_setups);
SE_MR_DL_prodSINR = zeros(tau_p,L,num_setups);
%for 500m layout
SE_MR_UL_2 = zeros(tau_p,L,num_setups);
SE_MR_DL_prodSINR_2 = zeros(tau_p,L,num_setups);

SE_RZF_UL = zeros(tau_p,L,num_setups);
SE_RZF_DL_prodSINR = zeros(tau_p,L,num_setups);
%for 500m layout
SE_RZF_UL_2 = zeros(tau_p,L,num_setups);
SE_RZF_DL_prodSINR_2 = zeros(tau_p,L,num_setups);

%SE_RZF_UL = zeros(tau_p,L,num_setups);

num_UEs = zeros(num_setups);
%% Set Network Layout

% Calculating for each setup for random UE locations
for setups = 1:num_setups
  setups
% Generate Simulation parameters
s = qd_simulation_parameters;
s.sample_density = 1;
s.use_absolute_delays = 1;

% Centre Frequency 30 GHz
center_frequency = 30e9;
wavelength = 3e8/center_frequency;
s.center_frequency = center_frequency;

% New Quadriga Layout
l = qd_layout(s);
% Layout for 500m*500m
l2 = qd_layout(s);

% Generate Base Stations
l.no_tx = L;
l2.no_tx = L;

% Base stations heights
l.tx_position(3,:) = 25;
l2.tx_position(3,:) = 25;

%Set the length in meters of the total square area
squareLength = 1000;
squareLength_2 = 500;
%Number of BSs per dimension
BSsPerDim = sqrt(L);

%% Base Station Deployment

% Minimum distance between BSs and UEs
min_distance = 35;

%Distance between BSs in vertical/horizontal direction
interSiteDistance = squareLength/BSsPerDim;
interSiteDistance_2 = squareLength_2/BSsPerDim;

% Deploy BSs on the 1000m grid
horizontal_axis = repmat(interSiteDistance/2:interSiteDistance:squareLength-interSiteDistance/2,[BSsPerDim 1]);
vertical_axis = horizontal_axis';
BS_positions = horizontal_axis(:) + 1i*vertical_axis(:);

% Deploy BSs on the 500m grid
horizontal_axis_2 = repmat(interSiteDistance_2/2:interSiteDistance_2:squareLength_2-interSiteDistance_2/2,[BSsPerDim 1]);
vertical_axis_2 = horizontal_axis_2';
BS_positions_2 = horizontal_axis_2(:) + 1i*vertical_axis_2(:);

% Give BS postions to Quadriga layout to set the positions
for i = 1:length(BS_positions)
    
    l.tx_position(1:2,i) = [real(BS_positions(i)) ; imag(BS_positions(i))];
end

% Give BS postions to Quadriga 500m layout to set the positions
for i = 1:length(BS_positions_2)
    
    l2.tx_position(1:2,i) = [real(BS_positions_2(i)) ; imag(BS_positions_2(i))];
end
%% Creating Antenna array (Planner) for each BS

M_V = 5; %Number of vertical antennas
M_H = M/M_V; %Number of antennas on each horizontal circle

%Number of vertical antennas
Mv = 10;

%Number of antennas in horizontal direction
Mh = M/Mv;

%Antenna Spacing 
ant_spacing = 1/2;

%Compute height of array
arrayHeight = (M_V-1)*ant_spacing*3e8/center_frequency;

%Polarization (dual polarized)
polarization_Indicator = 3;

    circumference = M_H/2*ant_spacing*3e8/center_frequency;
    radius = circumference/(2*pi);
    delta_angle = 2*pi/(M_H/2);
    
    %Go through all BSs
%     for b = 1:L
%         
%         %Create rectangular array of size M_V x M_H/2
%         l.tx_array(b).generate('3gpp-3d', 1, M_V, M_H/2, center_frequency, polarization_Indicator, 0, ant_spacing);
%         
%         %Place antennas on a circle and rotate radiation patters (while keeping co-located antennas together)
%         for i = 1:M_V
%             for j = 1:M_H/2
%                 indices = (i-1)*M_H + 2*j-1 : (i-1)*M_H + 2*j;
%                 angle = (j-1)*delta_angle;
%                 l.tx_array(b).element_position(1, indices) = radius*cos(angle);
%                 l.tx_array(b).element_position(2, indices) = radius*sin(angle);
%                 l.tx_array(b).element_position(3, indices) = (i-1)*ant_spacing*3e8/center_frequency - arrayHeight/2;
%                 l.tx_array(b).rotate_pattern(rad2deg(angle), 'z', indices, 0);
%             end
%         end
%     end
% Go through all cells

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

% Antenna array for 500m layout
for i = 1:L
    
    l2.tx_array(i).generate('3gpp-3d', 1,Mv, Mh, center_frequency, polarization_Indicator, 0, ant_spacing);
 
    for j = 1:Mv
        for k = 1:Mh
            index_2 = (j-1)*Mh + k;
            l2.tx_array(i).element_position(2,index_2) = mod(index_2-1,Mh)*ant_spacing*wavelength;
            l2.tx_array(i).element_position(3,index_2) = floor(index_2-1/Mh)*ant_spacing*wavelength;
        end
    end
    
end

l.set_scenario('3GPP_38.901_UMa_NLOS');
l2.set_scenario('3GPP_38.901_UMa_NLOS');
%% Drop UEs in the cells

user_positions = zeros(Kdrop,L);
users_perBS = zeros(L,1);
user_positions_2 = zeros(Kdrop,L);
users_perBS_2 = zeros(L,1);

% Maximum cell distance
maxCell_distance = interSiteDistance;
% Maximum cell distance for 500m layout
maxCell_distance_2 = interSiteDistance_2;

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

% Go through all cells (Drop UEs in 500m layout)
for i=1:L
    
    % Put users randomly in the cell
    while users_perBS_2(i)<Kdrop
        
        %Drop users
        UEremaining_2 = Kdrop-users_perBS_2(i);
        X_coord_2 = rand(UEremaining_2,1)*maxCell_distance_2 - maxCell_distance_2/2;
        Y_coord_2 = rand(UEremaining_2,1)*maxCell_distance_2 - maxCell_distance_2/2;
        XY_coord_2 = X_coord_2 + 1i*Y_coord_2;
        
        %Keep those UEs that are at minimum distance (35 meter) away from
        %BS
        XY_coord_2 = XY_coord_2(abs(XY_coord_2)>=min_distance);
        
        %Store new UEs
        user_positions_2(users_perBS_2(i)+1:users_perBS_2(i)+length(XY_coord_2),i) = XY_coord_2 + BS_positions_2(i);
        users_perBS_2(i) = users_perBS_2(i)+length(XY_coord_2);
        
    end
end

%% Receiver Characteristics

% Total number of UEs
UE_total = Kdrop*L;

% UE height (meters)
UE_height = 1.5;

% Generate UEs in the layout
l.no_rx = UE_total;
l2.no_rx = UE_total;

% UEs antenna type
l.rx_array.generate('omni');
l2.rx_array.generate('omni');


% Extract UE positions from user_positions matrix
reshape_user_positions = reshape(user_positions,[UE_total 1]);
% Extract UE positions from user_positions matrix for 500m layout
reshape_user_positions_2 = reshape(user_positions_2,[UE_total 1]);


% Get all users positions
for i = 1:length(reshape_user_positions)
    
    l.rx_position(1:3,i) = [real(reshape_user_positions(i)) ; imag(reshape_user_positions(i)) ; UE_height];
end

% Get all users positions for 500m layout
for i = 1:length(reshape_user_positions_2)
    
    l2.rx_position(1:3,i) = [real(reshape_user_positions_2(i)) ; imag(reshape_user_positions_2(i)) ; UE_height];
end

% UE tracks

for j=1:UE_total
    
    l.track(j).generate('linear',2,0); %Define a linear track of 2m and east direction
    l.track(j).scenario = '3GPP_38.901_UMa_NLOS'; %Select the Urban Macrocell NLOS scenario
    %t(j).movement_profile = [ 0,0.5 ; 5,0 ; 6,0 ; 10,2  ]';     % Generate movement profile
    
end

% UE tracks for 500m layout

for j=1:UE_total
    
    l2.track(j).generate('linear',2,0); %Define a linear track of 2m and east direction
    l2.track(j).scenario = '3GPP_38.901_UMa_NLOS'; %Select the Urban Microcell NLOS scenario
    %t(j).movement_profile = [ 0,0.5 ; 5,0 ; 6,0 ; 10,2  ]';     % Generate movement profile
    
end
%t.interpolate_positions(s.samples_per_meter);


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
noiseVar_dBm = -90;  %(dBm)
noise_Var = 10^(noiseVar_dBm/10);
noise_Std = sqrt(noise_Var);
%% Generate channel output and Spatial correlation Matrix

% Intialize variables to store values
H = zeros(M,num_subcarriers,Kmax,L,L);
R = zeros(M,M,Kmax,L,L);
users_perBS = zeros(L,1);
active_UEs = zeros(Kmax,L);



for i = 1:L
    i
    % Get Channels
    ch = l.get_channels;
    %dist = t.interpolate_movement( 1e-3 );                  % Get snapshot positions
    %ch = cn.interpolate(dist);
    
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

%% Generate channel output and Spatial correlation Matrix for 500m layout

% Intialize variables to store values
H_2 = zeros(M,num_subcarriers,Kmax,L,L);
R_2 = zeros(M,M,Kmax,L,L);
users_perBS_2 = zeros(L,1);
active_UEs_2 = zeros(Kmax,L);



for i = 1:L
    i
    % Get Channels
    ch_2 = l2.get_channels;
    %dist = t.interpolate_movement( 1e-3 );                  % Get snapshot positions
    %ch = cn.interpolate(dist);
    
    for j = 1:Kdrop
        
        h_k_2 = zeros(M,num_subcarriers,1,1,L);
        r_k_2 = zeros(M,M,1,1,L);
        
        
        for k = 1:L
            
            [~,ind_2] =  min(abs(user_positions_2(j,:)-BS_positions_2(k)));
            
            % Find the right user index
            user_ind_2 = find(user_positions_2(j,ind_2)==reshape_user_positions_2);
            
            % Select channel for individual user
            h_2 = ch_2(user_ind_2,k);
            
            H_freq_2 = h_2.fr(B, num_subcarriers);  % Getting frequency response
            
            % Getting Channels of individual user from all BSs
            H_freq_2 = mean(H_freq_2,4);
            h_k_2(:,:,1,1,k) = reshape(H_freq_2,[M num_subcarriers])/noise_Std;
            
            % Getting Spatial Correlation Matrix of individual user from all BSs
            r_k_2(:,:,1,1,k) = diag(mean(abs(h_k_2(:,:,1,1,k)).^2,2)/noise_Var);

        end
         %Determine which BS should serve the UE
        [max_value_2,BS_index_2] = max(mean(sum(abs(h_k_2(:,:,1,1,:)).^2,1),2));
        
        if users_perBS_2(BS_index_2)<Kmax
           
           users_perBS_2(BS_index_2) = users_perBS_2(BS_index_2) + 1; 
           H_2(:,:,randOrder(users_perBS_2(BS_index_2)+(pilotPattern(BS_index_2)-1)*Kmax,BS_index_2),BS_index_2,:) = h_k_2;
           R_2(:,:,randOrder(users_perBS_2(BS_index_2)+(pilotPattern(BS_index_2)-1)*Kmax,BS_index_2),BS_index_2,:) = r_k_2;
           active_UEs(randOrder(users_perBS(BS_index)+(pilotPattern(BS_index)-1)*Kmax,BS_index),BS_index) = 1; 
        end
    end
end

%% Channel Estimation for Downlink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H,R,num_subcarriers,M,tau_p,L,up_pwr);
% for 500m layout
[Hhat_LS_2,C_LS_2] = functionChannelEstimates_LS(H_2,R_2,num_subcarriers,M,tau_p,L,up_pwr);

%% Compute SINR for Downlink

% Active UEs
num_UEs(setups) = sum(active_UEs(:)); 

[signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE, interf_MMMSE] = functionComputeSINR_DL(H,Hhat_LS,C_LS,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr);
% for 500m layout
[signal_MR_2,interf_MR_2,signal_RZF_2,interf_RZF_2,signal_MMMSE_2, interf_MMMSE_2] = functionComputeSINR_DL(H_2,Hhat_LS_2,C_LS_2,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr);


prelogFactor = 2/3*(1-tau_p/tau_c);
%% Compute Spectral Efficiency

% Spectral Efficiency using Product SINR power control

SE_MR_DL_prodSINR(:,:,setups) = functionPowerOptimization_prodSINR(signal_MR ,interf_MR, txPwr, prelogFactor);

SE_RZF_DL_prodSINR(:,:,setups) = functionPowerOptimization_prodSINR( signal_RZF ,interf_RZF, txPwr,prelogFactor);

% Spectral Efficiency using Product SINR power control for 500m layout

SE_MR_DL_prodSINR_2(:,:,setups) = functionPowerOptimization_prodSINR(signal_MR_2 ,interf_MR_2, txPwr, prelogFactor);

SE_RZF_DL_prodSINR_2(:,:,setups) = functionPowerOptimization_prodSINR( signal_RZF_2 ,interf_RZF_2, txPwr,prelogFactor);

%SE_MMMSE_DL_prodSINR(:,:,setups) = functionPowerOptimization_prodSINR( signal_MMMSE ,interf_MMMSE, txPwr,prelogFactor);

% Spectral Efficiency without power control

% SE_MR_DL(:,:,setups) = function_compute_SE_without_powerAllocation(txPwr,signal_MR ,interf_MR, prelogFactor);
% 
% SE_RZF_DL(:,:,setups) = function_compute_SE_without_powerAllocation( txPwr,signal_RZF ,interf_RZF,prelogFactor);
% 
% SE_MMMSE_DL(:,:,setups) = function_compute_SE_without_powerAllocation( txPwr,signal_MMMSE ,interf_MMMSE,prelogFactor);

%% Computing Uplink Results

clear Hhat_LS C_LS

% 20 dB power control
Pwr_diff_dB = 20;

H_scaled = H;
R_scaled = R;

% Go through all Base stations
for i = 1:L
    
    % Computing average channel gains in dB
    betas = 10*log10(squeeze(sum(sum(R(:,:,active_UEs(:,i)==1,i,i),1),2)/M));
    
    % minimum from all betas
    betaj_min = min(betas);
    
    % beta(j,j,k)/beta(j,min)
    difference = betas - betaj_min;
    
    backoff = difference - Pwr_diff_dB;
    
    % find indices of active UEs
    activeUEs_index = find(active_UEs(:,i));
    
    backoff(backoff<0) = 0;
    
    for j = 1:length(activeUEs_index)
        
        H_scaled(:,:,activeUEs_index(j),i,:) = H(:,:,activeUEs_index(j),i,:)/ 10^(backoff(j)/20);
        R_scaled(:,:,activeUEs_index(j),i,:) = R(:,:,activeUEs_index(j),i,:)/ 10^(backoff(j)/10);
    end
    
end


%% Computing Uplink Results for 500m layout

clear Hhat_LS_2 C_LS_2

% 20 dB power control
Pwr_diff_dB_2 = 20;

H_scaled_2 = H_2;
R_scaled_2 = R_2;

% Go through all Base stations
for i = 1:L
    
    % Computing average channel gains in dB
    betas_2 = 10*log10(squeeze(sum(sum(R_2(:,:,active_UEs(:,i)==1,i,i),1),2)/M));
    
    % minimum from all betas
    betaj_min_2 = min(betas_2);
    
    % beta(j,j,k)/beta(j,min)
    difference_2 = betas_2 - betaj_min_2;
    
    backoff_2 = difference_2 - Pwr_diff_dB_2;
    
    % find indices of active UEs
    activeUEs_index_2 = find(active_UEs(:,i));
    
    backoff_2(backoff_2<0) = 0;
    
    for j = 1:length(activeUEs_index_2)
        
        H_scaled_2(:,:,activeUEs_index_2(j),i,:) = H_2(:,:,activeUEs_index_2(j),i,:)/ 10^(backoff_2(j)/20);
        R_scaled_2(:,:,activeUEs_index_2(j),i,:) = R_2(:,:,activeUEs_index_2(j),i,:)/ 10^(backoff_2(j)/10);
    end
    
end
%% Compute Channel Estimates for Uplink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H_scaled,R_scaled,num_subcarriers,M,tau_p,L,up_pwr);
% for 500m layout
[Hhat_LS_2,C_LS_2] = functionChannelEstimates_LS(H_scaled_2,R_scaled_2,num_subcarriers,M,tau_p,L,up_pwr);


%% Compute SE for Uplink

%[SE_MR_UL,SE_RZF_UL] = functionComputeSE_UL(Hhat_LS,C_LS,R_scaled,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr);
[SE_MR,SE_RZF] = functionComputeSE_UL_impairments(H_scaled,Hhat_LS,C_LS,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr,1,1);
% for 500m layout
[SE_MR_2,SE_RZF_2] = functionComputeSE_UL_impairments(H_scaled_2,Hhat_LS_2,C_LS_2,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr,1,1);

SE_MR_UL(:,:,setups) = 1/3*SE_MR;
SE_RZF_UL(:,:,setups) = 1/3*SE_RZF;
% for 500m layout
SE_MR_UL_2(:,:,setups) = 1/3*SE_MR_2;
SE_RZF_UL_2(:,:,setups) = 1/3*SE_RZF_2;

end


%% Plotting Simulation results for Downlink

total_users = sum(num_UEs(:));

% % Sorting results in ascending order without power control
% SE_MR_DL = sort(SE_MR_DL(SE_MR_DL>0));
% SE_RZF_DL = sort(SE_RZF_DL(SE_RZF_DL>0));
% SE_MMMSE_DL = sort(SE_MMMSE_DL(SE_MMMSE_DL>0));

% Sorting results in ascending order for product SINR
SE_MR_DL_prodSINR = sort(SE_MR_DL_prodSINR(SE_MR_DL_prodSINR>0));
SE_RZF_DL_prodSINR = sort(SE_RZF_DL_prodSINR(SE_RZF_DL_prodSINR>0));

% Sorting results in ascending order for product SINR for 500m layout
SE_MR_DL_prodSINR_2 = sort(SE_MR_DL_prodSINR_2(SE_MR_DL_prodSINR_2>0));
SE_RZF_DL_prodSINR_2 = sort(SE_RZF_DL_prodSINR_2(SE_RZF_DL_prodSINR_2>0));



% Product SINR Results

figure(2)
plot((effective_BW/1e6)*SE_MR_DL_prodSINR,linspace(0,1,total_users),'r','LineWidth',1);

hold on

% SE_RZF_DL(numel(SE_MR_DL)) = 0;
plot((effective_BW/1e6)*SE_RZF_DL_prodSINR,linspace(0,1,total_users),'k--','LineWidth',1);

hold on

plot((effective_BW/1e6)*SE_MR_DL_prodSINR_2,linspace(0,1,total_users),'b:','LineWidth',1);

hold on

plot((effective_BW/1e6)*SE_RZF_DL_prodSINR_2,linspace(0,1,total_users),'m-.','LineWidth',1);



legend('SINR MR_1000mby1000m layout (M=100 for Planner array)','SINR RZF_1000mby1000m layout (M=100 for Planner array)','SINR MR_500m by 500m layout','SINR RZF_500m by 500m layout','Location','Southeast');
xlabel('Downlink throghout per UE [Mbits/s]');
ylabel('CDF');
title('DL throughput per UE using product SINR power allocation with 5 user per cell (1000m by 1000m vs 500m by 500m)(Freq=30 GHz)')
grid on

%% Plotting Simulation results for Uplink

SE_MR_UL = sort(SE_MR_UL(SE_MR_UL>0));
SE_RZF_UL = sort(SE_RZF_UL(SE_RZF_UL>0));
% for 500m layout
SE_MR_UL_2 = sort(SE_MR_UL_2(SE_MR_UL_2>0));
SE_RZF_UL_2 = sort(SE_RZF_UL_2(SE_RZF_UL_2>0));


%Uplink spectral efficeiency results
figure(1)
%SE_MR_UL(numel(SE_MR_DL_prodSINR)) = 0;
plot((effective_BW/1e6)*SE_MR_UL,linspace(0,1,total_users),'r','LineWidth',1);

hold on

%SE_RZF_UL(numel(SE_RZF_DL_prodSINR)) = 0;
plot((effective_BW/1e6)*SE_RZF_UL,linspace(0,1,total_users),'k--','LineWidth',1);

hold on

plot((effective_BW/1e6)*SE_MR_UL_2,linspace(0,1,total_users),'b:','LineWidth',1);

hold on

plot((effective_BW/1e6)*SE_RZF_UL_2,linspace(0,1,total_users),'m-.','LineWidth',1);


xlim([0 50])

legend('SINR MR_1000m by 1000m layout (M=100 for Planner array)','SINR RZF_1000m by 1000m layout (M=100 for Planner array)','SINR MR_500m by 500m layout (M=100 for Planner array)','SINR RZF_500m by 500m layout (M=100 for Planner array)','Location','Southeast');
xlabel('Uplink throghout per UE [Mbits/s]');
ylabel('CDF');
title('UL throughput per UE with 5 user per cell (1000m by 1000m vs 500m by 500m)(Freq=30 GHz)')
grid on


%% Visualize antenna patterns

% close all
% 
% set(0,'defaultTextFontSize', 18)                      	% Default Font Size
% set(0,'defaultAxesFontSize', 18)                     	% Default Font Size
% set(0,'defaultAxesFontName','Times')               	    % Default Font Type
% set(0,'defaultTextFontName','Times')                 	% Default Font Type
% set(0,'defaultFigurePaperPositionMode','auto')       	% Default Plot position
% set(0,'DefaultFigurePaperType','<custom>')             	% Default Paper Type
% set(0,'DefaultFigurePaperSize',[14.5 7.3])            	% Default Paper Size
% 
% [ map,x_coords,y_coords] = l2.power_map( '3GPP_38.901_UMa_NLOS','quick',5,0,500,0,500,1.5,30);
% P = 10*log10(sum( abs( cat(3,map{:}) ).^2 ,3));         % Total received power
% dim = length(x_coords);
% P = reshape(P(:,:,1,1),[dim dim]);
% 
% l2.visualize([],[],0);                                   % Show BS and MT positions on the map
% hold on
% imagesc( x_coords, y_coords, P );                     % Plot the received power
% hold on
% yline(250);
% hold on
% yline(500);
% hold on
% yline(750);
% hold on
% xline(250);
% hold on
% xline(500);
% hold on
% xline(750);
% axis([0 500 0 500])                               % Plot size
% caxis( max(P(:)) + [-20 0] )                            % Color range 
% colmap = colormap;
% colormap( colmap*0.5 + 0.5 );                           % Adjust colors to be "lighter"
% set(gca,'layer','top')                                  % Show grid on top of the map
% title('Antenna orientation');                 % Set plot title
