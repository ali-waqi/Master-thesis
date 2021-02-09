function SE = function_compute_SE( rho ,signal, interference,prelog_factor)

% This fuction computes the spectral efficiency using Theorum (4.6) of book
% Massive MIMO Networks

% INPUT:

%signal       = K x L matrix where element (k,j) is a_jk in (7.2)
%interference = K x L x K x L matrix where (l,i,j,k) is b_lijk in (7.3)
%Pmax         = Maximum transmit power per BS
%prelogFactor = Prelog factor

% OUTPUT:

%SE           = K x L matrix where element (k,j) is the downlink SE of UE k in cell j
%               using the max-min power allocation solution




% Getting number of users
K = size(signal,1);

% Getting number of cells
L = size(signal,2);

% To store Spectral efficiency matrix
SE = zeros(K,L);

% Go through all cells
for j = 1:L
    
    % Go through all users
    for k = 1:K
        
        SE(k,j) = prelog_factor * log2(1+(rho(k,j) * signal(k,j)/sum(sum(rho .* interference(:,:,k,j)))+1));
    end
end
end