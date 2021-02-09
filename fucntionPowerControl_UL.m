function [] = fucntionPowerControl_UL(R,active_UEs,M,K,L,uplink_power)

for i = 1:L
    betas = trace(R(:,:,active_UEs(:,i),:,i))/M;
end
end