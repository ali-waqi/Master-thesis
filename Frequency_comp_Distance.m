%% Set simulation parameters
close all
clc

% Random seed generator
rng('shuffle');

% Number of users to be droped in cell 
Kdrop = 5;

% Maximum users served by BS
K_max = 6;

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
txPwr = [1000,1000,1000];      %(mW)

% Pilot reuse factor
f = 2;

% Communication Bandwidths (20 MHz, 40 MHz, 100 MHz)
B = [20e6,40e6,100e6];

% Effective bandwidth after 5% overhead
effective_BW = B(1)*0.95;
effective_BW_6GHz = B(2)*0.95;
effective_BW_30GHz = B(3)*0.95;


% Length of pilot Signals
tau_p = K_max*f;

% Length of coherence block
tau_c = 190;
tau_c_30GHz = 1428;

% Noise variance in dBm
noise = -94;

% Polarization (dual polarized)
polarization = 3;

% Square Length
square_length = 4000;

% Distances vector
distance_steps = 170:10:230;

% Carrier frequencies (1 GHz, 6 GHz, 30 GHz)
carrier_freq = [1e9,6e9,30e9];

%% Calculate DL and UL throughput of UEs for number of setups

% Intialize to store different runs values
sum_SE_MR_UL = zeros(length(distance_steps),num_setups);
sum_SE_MR_DL_prodSINR = zeros(length(distance_steps),num_setups);
sum_SE_MR_UL_6GHz = zeros(length(distance_steps),num_setups);
sum_SE_MR_DL_prodSINR_6GHz = zeros(length(distance_steps),num_setups);
sum_SE_MR_UL_30GHz = zeros(length(distance_steps),num_setups);
sum_SE_MR_DL_prodSINR_30GHz = zeros(length(distance_steps),num_setups);


sum_SE_RZF_UL = zeros(length(distance_steps),num_setups);
sum_SE_RZF_DL_prodSINR = zeros(length(distance_steps),num_setups);
sum_SE_RZF_UL_6GHz = zeros(length(distance_steps),num_setups);
sum_SE_RZF_DL_prodSINR_6GHz = zeros(length(distance_steps),num_setups);
sum_SE_RZF_UL_30GHz = zeros(length(distance_steps),num_setups);
sum_SE_RZF_DL_prodSINR_30GHz = zeros(length(distance_steps),num_setups);


%% Set Network Layout
for setups = 1:num_setups

%Output simulation progress
disp([num2str(setups) ' setups out of ' num2str(num_setups)]);


% Calculating for each setup for random UE locations
for distance = 1:length(distance_steps)

%Output simulation progress
disp([num2str(distance) ' BS-UE distance out of ' num2str(length(distance_steps))]);


%% Compute Throughput for 1 GHz frequency

% Set Quadriga layout for 1 GHz frequency
[H,R,active_UEs] = function_Setup_Distance(distance_steps(distance),square_length,L,Kdrop,B(1),noise,K_max,f,M,polarization,carrier_freq(1),num_subcarriers); 

%% Channel Estimation for Downlink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H,R,num_subcarriers,M,tau_p,L,up_pwr);

%% Compute SINR for Downlink
 
% Compute signal and interference for SE computation
[signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE, interf_MMMSE] = functionComputeSINR_DL(H,Hhat_LS,C_LS,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr);
 
prelogFactor = 2/3*(1-tau_p/tau_c);
%% Compute Spectral Efficiency

% clear large matrices
clear Hhat_LS C_LS signal_MMMSE interf_MMMSE 

% Spectral Efficiency using Product SINR power control
 
 SE_MR_DL = functionPowerOptimization_prodSINR(signal_MR ,interf_MR, txPwr(1), prelogFactor);
 sum_SE_MR_DL_prodSINR(distance,setups) = sum(rmmissing(SE_MR_DL(:)));
 
 SE_RZF_DL = functionPowerOptimization_prodSINR( signal_RZF ,interf_RZF, txPwr(1),prelogFactor);
 sum_SE_RZF_DL_prodSINR(distance,setups) = sum(rmmissing(SE_RZF_DL(:)));

%% Heristic power control

% 20 dB power control
delta_dB = 20;

[H_scaled,R_scaled] = function_powerControl_UL(L,M,delta_dB,H,R,active_UEs);


%% Compute Channel Estimates for Uplink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H_scaled,R_scaled,num_subcarriers,M,tau_p,L,up_pwr);

%% Compute SE for Uplink

[SE_MR,SE_RZF] = functionComputeSE_UL_impairments(H_scaled,Hhat_LS,C_LS,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr,1,1);
 
SE_MR_UL = 1/3*SE_MR;
sum_SE_MR_UL(distance,setups) = sum(SE_MR_UL(:));

SE_RZF_UL = 1/3*SE_RZF;
sum_SE_RZF_UL(distance,setups) = sum(SE_RZF_UL(:));


%% Compute Throughput for 6 GHz frequency

% clear large matrices
clear R H Hhat_LS C_LS

% Set Quadriga layout for 6 GHz frequency
[H,R,active_UEs] = function_Setup_Distance(distance_steps(distance),square_length,L,Kdrop,B(2),noise,K_max,f,M,polarization,carrier_freq(2),num_subcarriers); 

%% Channel Estimation for Downlink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H,R,num_subcarriers,M,tau_p,L,up_pwr);

%% Compute SINR for Downlink
 
% Compute signal and interference for SE computation
[signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE, interf_MMMSE] = functionComputeSINR_DL(H,Hhat_LS,C_LS,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr);
 
prelogFactor = 2/3*(1-tau_p/tau_c);
%% Compute Spectral Efficiency

% clear large matrices
clear Hhat_LS C_LS signal_MMMSE interf_MMMSE 

% Spectral Efficiency using Product SINR power control
 
 SE_MR_DL_6GHz =  functionPowerOptimization_prodSINR(signal_MR ,interf_MR, txPwr(2), prelogFactor);
 sum_SE_MR_DL_prodSINR_6GHz(distance,setups) = sum(rmmissing(SE_MR_DL_6GHz(:)));
 
 SE_RZF_DL_6GHz = functionPowerOptimization_prodSINR( signal_RZF ,interf_RZF, txPwr(2),prelogFactor);
 sum_SE_RZF_DL_prodSINR_6GHz(distance,setups) = sum(rmmissing(SE_RZF_DL_6GHz(:)));

%% Heristic power control

% 20 dB power control
delta_dB = 20;

[H_scaled,R_scaled] = function_powerControl_UL(L,M,delta_dB,H,R,active_UEs);


%% Compute Channel Estimates for Uplink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H_scaled,R_scaled,num_subcarriers,M,tau_p,L,up_pwr);

%% Compute SE for Uplink

[SE_MR,SE_RZF] = functionComputeSE_UL_impairments(H_scaled,Hhat_LS,C_LS,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr,1,1);
 
SE_MR_UL_6GHz = 1/3*SE_MR;
sum_SE_MR_UL_6GHz(distance,setups) = sum(SE_MR_UL_6GHz(:));

SE_RZF_UL_6GHz =1/3*SE_RZF;
sum_SE_RZF_UL_6GHz(distance,setups) = sum(SE_RZF_UL_6GHz(:));

%% Compute Throughput for 30 GHz frequency

% Clear large matrices
clear R H Hhat_LS C_LS

% Set Quadriga layout for 30 GHz frequency
[H,R,active_UEs] = function_Setup_Distance(distance_steps(distance),square_length,L,Kdrop,B(3),noise,K_max,f,M,polarization,carrier_freq(3),num_subcarriers); 

%% Channel Estimation for Downlink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H,R,num_subcarriers,M,tau_p,L,up_pwr);

%% Compute SINR for Downlink
 
% Compute signal and interference for SE computation
[signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE, interf_MMMSE] = functionComputeSINR_DL(H,Hhat_LS,C_LS,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr);
 
prelogFactor = 2/3*(1-tau_p/tau_c);
%% Compute Spectral Efficiency

% clear large matrices
clear Hhat_LS C_LS signal_MMMSE interf_MMMSE 

% Spectral Efficiency using Product SINR power control
 
 SE_MR_DL_30GHz = functionPowerOptimization_prodSINR(signal_MR ,interf_MR, txPwr(3), prelogFactor);
 sum_SE_MR_DL_prodSINR_30GHz(distance,setups) = sum(rmmissing(SE_MR_DL_30GHz(:)));
 
 SE_RZF_DL_30GHz = functionPowerOptimization_prodSINR( signal_RZF ,interf_RZF, txPwr(3),prelogFactor);
 sum_SE_RZF_DL_prodSINR_30GHz(distance,setups) = sum(rmmissing(SE_RZF_DL_30GHz(:)));

%% Heristic power control

% 20 dB power control
delta_dB = 20;

[H_scaled,R_scaled] = function_powerControl_UL(L,M,delta_dB,H,R,active_UEs);


%% Compute Channel Estimates for Uplink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H_scaled,R_scaled,num_subcarriers,M,tau_p,L,up_pwr);

%% Compute SE for Uplink

[SE_MR,SE_RZF] = functionComputeSE_UL_impairments(H_scaled,Hhat_LS,C_LS,tau_c,tau_p,num_subcarriers,M,tau_p,L,up_pwr,1,1);

% clear large matrices
clear Hhat_LS C_LS R H 
 
SE_MR_UL_30GHz = 1/3*SE_MR;
sum_SE_MR_UL_30GHz(distance,setups) = sum(SE_MR_UL_30GHz(:));

SE_RZF_UL_30GHz = 1/3*SE_RZF;
sum_SE_RZF_UL_30GHz(distance,setups) = sum(SE_RZF_UL_30GHz(:));

% Setups for end
end

% distance for end
end
%% Plotting Simulation results for Downlink

% Product SINR Results
figure(2)

% plot(distance_steps,mean(sum_SE_MR_DL_prodSINR,2),'r-','LineWidth',1);
%  
% hold on
%  
% plot(distance_steps,mean(sum_SE_RZF_DL_prodSINR,2),'r--','LineWidth',1);
%  
% %Product SINR Results for 6GHZ freq
% hold on
% 
% plot(distance_steps,mean(sum_SE_MR_DL_prodSINR_6GHz,2),'k-','LineWidth',1);
%  
% hold on
% 
% plot(distance_steps,mean(sum_SE_RZF_DL_prodSINR_6GHz,2),'k--','LineWidth',1);
%  
% %Product SINR Results for 30GHz freq
% hold on

plot(distance_steps,mean(sum_SE_MR_DL_prodSINR_30GHz,2),'b-','LineWidth',1);
 
hold on
 
plot(distance_steps,mean(sum_SE_RZF_DL_prodSINR_30GHz,2),'b--','LineWidth',1);
 
 
legend('MR (1 GHz)','RZF(1 GHz) ','MR(6 GHz) ','RZF(6 GHz) ','MR(30 GHz) ','RZF(30 GHz)','Location','Southeast');
xlabel('Distance between UEs and BS');
ylabel('Average sum spectral efficiency');
grid on

%% Plotting Simulation results for Uplink

% Uplink spectral efficeiency results
figure(1)

% plot(distance_steps,mean(sum_SE_MR_UL,2),'r-','LineWidth',1);
% 
% hold on
% 
% plot(distance_steps,mean(sum_SE_RZF_UL,2),'r--','LineWidth',1);
% 
% hold on
% 
% plot(distance_steps,mean(sum_SE_MR_UL_6GHz,2),'k-','LineWidth',1);
% 
% hold on
% 
% plot(distance_steps,mean(sum_SE_RZF_UL_6GHz,2),'k--','LineWidth',1);
% 
% hold on

plot(distance_steps,mean(sum_SE_MR_UL_30GHz,2),'b-','LineWidth',1);

hold on

plot(distance_steps,mean(sum_SE_RZF_UL_30GHz,2),'b--','LineWidth',1);


legend('MR(1 GHz) ','RZF(1 GHz) ','MR(6 GHz) ','RZF(6 GHz) ','MR(30 GHz) ','RZF(30 GHz) ','Location','Southwest');
xlabel('Distance between UEs and BS');
ylabel('Average sum spectral efficiency');
grid on

