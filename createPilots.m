function [phi] = createPilots(numUT,numSamples)
% Create pilot matrix
% Input:    numUT - number of UTs (K)
%           numSamples - number of samples 
% Output:   phi - pilot matrix (NxK)
% where N=number of samples and K=number of UTs. 

    %oneSine = sin(2*pi*(1:numSamples)/numSamples);
    phi = eye(numUT,numSamples)';
    
    %phi = sqrt(numUT*numSamples)*phi; % Normalizing phi
    
end