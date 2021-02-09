function[H_scaled, R_scaled] = function_powerControl_UL(L,M,delta_dB,H,R,active_UEs)

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
    
    backoff = difference - delta_dB;
    
    % find indices of active UEs
    activeUEs_index = find(active_UEs(:,i));
    
    backoff(backoff<0) = 0;
    
    for j = 1:length(activeUEs_index)
        
        H_scaled(:,:,activeUEs_index(j),i,:) = H(:,:,activeUEs_index(j),i,:)/ 10^(backoff(j)/20);
        R_scaled(:,:,activeUEs_index(j),i,:) = R(:,:,activeUEs_index(j),i,:)/ 10^(backoff(j)/10);
    end
    
end

end