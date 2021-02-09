function [G] = estimateChannel(Y,phi)
% Estimation of the channel 
% Input:    Y - the recorded pilot from the BS (MxN)
%           phi - pilot (NxK)
% Output:   G - channel matrix (MxK)
% where M=number of BSA, N=number of samples, K=number of UTs. 

    G = Y*phi; 

end

