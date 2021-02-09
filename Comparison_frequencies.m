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
num_subcarriers = [400,800,1200];

% Number of setups with random UE locations
num_setups = 30;

% Transmitter power
txPwr = [1000,2000,5000];      %(mW)

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
tau_c = [1140,190,38];

% Noise variance in dBm
noise = [-94,-91,-87];

% Polarization (dual polarized)
polarization = 3;

% Square Length
square_length = 1000;

% Carrier frequencies (1 GHz, 6 GHz, 30 GHz)
carrier_freq = [1e9,6e9,30e9];

%% Calculate DL and UL throughput of UEs for number of setups

% Initialize variable to store simulations
SE_MR_UL = zeros(tau_p,L,num_setups);
SE_MR_DL_prodSINR = zeros(tau_p,L,num_setups);
SE_MR_UL_6GHz = zeros(tau_p,L,num_setups);
SE_MR_DL_prodSINR_6GHz = zeros(tau_p,L,num_setups);
SE_MR_UL_30GHz = zeros(tau_p,L,num_setups);
SE_MR_DL_prodSINR_30GHz = zeros(tau_p,L,num_setups);


SE_RZF_UL = zeros(tau_p,L,num_setups);
SE_RZF_DL_prodSINR = zeros(tau_p,L,num_setups);
SE_RZF_UL_6GHz = zeros(tau_p,L,num_setups);
SE_RZF_DL_prodSINR_6GHz = zeros(tau_p,L,num_setups);
SE_RZF_UL_30GHz = zeros(tau_p,L,num_setups);
SE_RZF_DL_prodSINR_30GHz = zeros(tau_p,L,num_setups);

num_UEs = zeros(num_setups);
num_UEs_6GHz = zeros(num_setups);
num_UEs_30GHz = zeros(num_setups);
%% Set Network Layout

% Calculating for each setup for random UE locations
for setups = 1:num_setups

%Output simulation progress
disp([num2str(setups) ' setups out of ' num2str(num_setups)]);

%% Compute Throughput for 1 GHz frequency

% Set Quadriga layout for 1 GHz frequency
[H,R,active_UEs] = function_Setup(square_length,L,Kdrop,B(1),noise(1),K_max,f,M,polarization,carrier_freq(1),num_subcarriers(1)); 

% Active UEs
num_UEs(setups) = sum(active_UEs(:));
%% Channel Estimation for Downlink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H,R,num_subcarriers(1),M,tau_p,L,up_pwr);

%% Compute SINR for Downlink
 
% Compute signal and interference for SE computation
[signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE, interf_MMMSE] = functionComputeSINR_DL(H,Hhat_LS,C_LS,tau_c(1),tau_p,num_subcarriers(1),M,tau_p,L,up_pwr);
 
prelogFactor = 2/3*(1-tau_p/tau_c(1));
%% Compute Spectral Efficiency

% clear large matrices
clear Hhat_LS C_LS signal_MMMSE interf_MMMSE 

% Spectral Efficiency using Product SINR power control
 
 SE_MR_DL_prodSINR(:,:,setups) = functionPowerOptimization_prodSINR(signal_MR ,interf_MR, txPwr(1), prelogFactor);
 
 SE_RZF_DL_prodSINR(:,:,setups) = functionPowerOptimization_prodSINR( signal_RZF ,interf_RZF, txPwr(1),prelogFactor);

%% Heristic power control

% 20 dB power control
delta_dB = 20;

[H_scaled,R_scaled] = function_powerControl_UL(L,M,delta_dB,H,R,active_UEs);


%% Compute Channel Estimates for Uplink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H_scaled,R_scaled,num_subcarriers(1),M,tau_p,L,up_pwr);

%% Compute SE for Uplink

[SE_MR,SE_RZF] = functionComputeSE_UL_impairments(H_scaled,Hhat_LS,C_LS,tau_c(1),tau_p,num_subcarriers(1),M,tau_p,L,up_pwr,1,1);
 
SE_MR_UL(:,:,setups) = 1/3*SE_MR;
SE_RZF_UL(:,:,setups) = 1/3*SE_RZF;



%% Compute Throughput for 6 GHz frequency

% clear large matrices
clear R H Hhat_LS C_LS

% Set Quadriga layout for 6 GHz frequency
[H,R,active_UEs] = function_Setup(square_length,L,Kdrop,B(2),noise(2),K_max,f,M,polarization,carrier_freq(2),num_subcarriers(2)); 

% Active UEs
num_UEs_6GHz(setups) = sum(active_UEs(:));
%% Channel Estimation for Downlink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H,R,num_subcarriers(2),M,tau_p,L,up_pwr);

%% Compute SINR for Downlink
 
% Compute signal and interference for SE computation
[signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE, interf_MMMSE] = functionComputeSINR_DL(H,Hhat_LS,C_LS,tau_c(2),tau_p,num_subcarriers(2),M,tau_p,L,up_pwr);
 
prelogFactor = 2/3*(1-tau_p/tau_c(2));
%% Compute Spectral Efficiency

% clear large matrices
clear Hhat_LS C_LS signal_MMMSE interf_MMMSE 

% Spectral Efficiency using Product SINR power control
 
 SE_MR_DL_prodSINR_6GHz(:,:,setups) = functionPowerOptimization_prodSINR(signal_MR ,interf_MR, txPwr(2), prelogFactor);
 
 SE_RZF_DL_prodSINR_6GHz(:,:,setups) = functionPowerOptimization_prodSINR( signal_RZF ,interf_RZF, txPwr(2),prelogFactor);

%% Heristic power control

% 20 dB power control
delta_dB = 20;

[H_scaled,R_scaled] = function_powerControl_UL(L,M,delta_dB,H,R,active_UEs);


%% Compute Channel Estimates for Uplink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H_scaled,R_scaled,num_subcarriers(2),M,tau_p,L,up_pwr);

%% Compute SE for Uplink

[SE_MR,SE_RZF] = functionComputeSE_UL_impairments(H_scaled,Hhat_LS,C_LS,tau_c(2),tau_p,num_subcarriers(2),M,tau_p,L,up_pwr,1,1);
 
SE_MR_UL_6GHz(:,:,setups) = 1/3*SE_MR;
SE_RZF_UL_6GHz(:,:,setups) = 1/3*SE_RZF;


%% Compute Throughput for 30 GHz frequency

% Clear large matrices
clear R H Hhat_LS C_LS

% Set Quadriga layout for 30 GHz frequency
[H,R,active_UEs] = function_Setup(square_length,L,Kdrop,B(3),noise(3),K_max,f,M,polarization,carrier_freq(3),num_subcarriers(3)); 

% Active UEs
num_UEs_30GHz(setups) = sum(active_UEs(:));
%% Channel Estimation for Downlink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H,R,num_subcarriers(3),M,tau_p,L,up_pwr);

%% Compute SINR for Downlink
 
% Compute signal and interference for SE computation
[signal_MR,interf_MR,signal_RZF,interf_RZF,signal_MMMSE, interf_MMMSE] = functionComputeSINR_DL(H,Hhat_LS,C_LS,tau_c(3),tau_p,num_subcarriers(3),M,tau_p,L,up_pwr);
 
prelogFactor = 2/3*(1-tau_p/tau_c(3));
%% Compute Spectral Efficiency

% clear large matrices
clear Hhat_LS C_LS signal_MMMSE interf_MMMSE 

% Spectral Efficiency using Product SINR power control
 
 SE_MR_DL_prodSINR_30GHz(:,:,setups) = functionPowerOptimization_prodSINR(signal_MR ,interf_MR, txPwr(3), prelogFactor);
 
 SE_RZF_DL_prodSINR_30GHz(:,:,setups) = functionPowerOptimization_prodSINR( signal_RZF ,interf_RZF, txPwr(3),prelogFactor);

%% Heristic power control

% 20 dB power control
delta_dB = 20;

[H_scaled,R_scaled] = function_powerControl_UL(L,M,delta_dB,H,R,active_UEs);


%% Compute Channel Estimates for Uplink

[Hhat_LS,C_LS] = functionChannelEstimates_LS(H_scaled,R_scaled,num_subcarriers(3),M,tau_p,L,up_pwr);

%% Compute SE for Uplink

[SE_MR,SE_RZF] = functionComputeSE_UL_impairments(H_scaled,Hhat_LS,C_LS,tau_c(3),tau_p,num_subcarriers(3),M,tau_p,L,up_pwr,1,1);

% clear large matrices
clear Hhat_LS C_LS R H 
 
SE_MR_UL_30GHz(:,:,setups) = 1/3*SE_MR;
SE_RZF_UL_30GHz(:,:,setups) = 1/3*SE_RZF;


end
%% Plotting Simulation results for Downlink

total_users = sum(num_UEs(:));
total_users_6GHz = sum(num_UEs_6GHz(:));
total_users_30GHz = sum(num_UEs_30GHz(:));

% Sorting results in ascending order for product SINR
SE_MR_DL_prodSINR = sort(SE_MR_DL_prodSINR(SE_MR_DL_prodSINR>0));
SE_RZF_DL_prodSINR = sort(SE_RZF_DL_prodSINR(SE_RZF_DL_prodSINR>0));

% Sorting results in ascending order for product SINR for 6GHz freq
SE_MR_DL_prodSINR_6GHz = sort(SE_MR_DL_prodSINR_6GHz(SE_MR_DL_prodSINR_6GHz>0));
SE_RZF_DL_prodSINR_6GHz = sort(SE_RZF_DL_prodSINR_6GHz(SE_RZF_DL_prodSINR_6GHz>0));

% Sorting results in ascending order for product SINR for 30GHz freq
SE_MR_DL_prodSINR_30GHz = sort(SE_MR_DL_prodSINR_30GHz(SE_MR_DL_prodSINR_30GHz>0));
SE_RZF_DL_prodSINR_30GHz = sort(SE_RZF_DL_prodSINR_30GHz(SE_RZF_DL_prodSINR_30GHz>0));

% Product SINR Results
figure(2)
plot((effective_BW/1e6)*SE_MR_DL_prodSINR,linspace(0,1,total_users),'r-','LineWidth',1);
 
 hold on
 
 % SE_RZF_DL(numel(SE_MR_DL)) = 0;
 plot((effective_BW/1e6)*SE_RZF_DL_prodSINR,linspace(0,1,total_users),'r--','LineWidth',1);
 
 % Product SINR Results for 6GHZ freq
 hold on

% SE_MR_DL_prodSINR_6GHz(length(SE_RZF_DL_prodSINR)+1:end) = [];
 plot((effective_BW_6GHz/1e6)*SE_MR_DL_prodSINR_6GHz,linspace(0,1,total_users_6GHz),'k-','LineWidth',1);
 
 hold on
 
% SE_RZF_DL_prodSINR_6GHz(length(SE_RZF_DL_prodSINR)+1:end) = [];
 plot((effective_BW_6GHz/1e6)*SE_RZF_DL_prodSINR_6GHz,linspace(0,1,total_users_6GHz),'k--','LineWidth',1);
 
% Product SINR Results for 30GHz freq
 hold on
% SE_MR_DL_prodSINR_30GHz(length(SE_RZF_DL_prodSINR)+1:end) = [];
 plot((effective_BW_30GHz/1e6)*SE_MR_DL_prodSINR_30GHz,linspace(0,1,total_users_30GHz),'b-','LineWidth',1);
 
 hold on
 
% SE_RZF_DL_prodSINR_30GHz(length(SE_RZF_DL_prodSINR)+1:end) = [];
 plot((effective_BW_30GHz/1e6)*SE_RZF_DL_prodSINR_30GHz,linspace(0,1,total_users_30GHz),'b--','LineWidth',1);
 
 
 legend('MR (1 GHz)','RZF(1 GHz) ','MR(6 GHz)','RZF(6 GHz)','MR(30 GHz)','RZF(30 GHz)','Location','Southeast','FontSize',20);
 ax = gca;
 ax.FontSize = 26;
 xlabel('DL throughput per UE [Mbits/s]','FontSize',32);
 ylabel('CDF','FontSize',32);
 grid on

%% Plotting Simulation results for Uplink


SE_MR_UL = sort(SE_MR_UL(SE_MR_UL>0));
SE_RZF_UL = sort(SE_RZF_UL(SE_RZF_UL>0));

%for 6GHz freq
SE_MR_UL_6GHz = sort(SE_MR_UL_6GHz(SE_MR_UL_6GHz>0));
SE_RZF_UL_6GHz = sort(SE_RZF_UL_6GHz(SE_RZF_UL_6GHz>0));

%for 30GHz freq
SE_MR_UL_30GHz = sort(SE_MR_UL_30GHz(SE_MR_UL_30GHz>0));
SE_RZF_UL_30GHz = sort(SE_RZF_UL_30GHz(SE_RZF_UL_30GHz>0));


% Uplink spectral efficeiency results
figure(1)
%SE_MR_UL(numel(SE_MR_DL_prodSINR)) = 0;
plot((effective_BW/1e6)*SE_MR_UL,linspace(0,1,total_users),'r-','LineWidth',1);

hold on

%SE_RZF_UL(numel(SE_RZF_DL_prodSINR)) = 0;
plot((effective_BW/1e6)*SE_RZF_UL,linspace(0,1,total_users),'r--','LineWidth',1);

hold on

plot((effective_BW_6GHz/1e6)*SE_MR_UL_6GHz,linspace(0,1,total_users_6GHz),'k-','LineWidth',1);


hold on

%SE_RZF_UL(numel(SE_RZF_DL_prodSINR)) = 0;
plot((effective_BW_6GHz/1e6)*SE_RZF_UL_6GHz,linspace(0,1,total_users_6GHz),'k--','LineWidth',1);

hold on

plot((effective_BW_30GHz/1e6)*SE_MR_UL_30GHz,linspace(0,1,total_users_30GHz),'b-','LineWidth',1);

hold on

%SE_RZF_UL(numel(SE_RZF_DL_prodSINR)) = 0;
plot((effective_BW_30GHz/1e6)*SE_RZF_UL_30GHz,linspace(0,1,total_users_30GHz),'b--','LineWidth',1);


legend('MR(1 GHz)','RZF(1 GHz) ','MR(6 GHz)','RZF(6 GHz)','MR(30 GHz)','RZF(30 GHz)','Location','Southeast','FontSize',20);
ax = gca;
ax.FontSize = 26;
xlabel('UL throughput per UE [Mbits/s]','FontSize',32);
ylabel('CDF','FontSize',32);
grid on

