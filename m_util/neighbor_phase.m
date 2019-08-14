function NFC = neighbor_phase(data, fs, winLen, stepLen)
%NFC = neighbor_phase(data, fs, winLen, stepLen)

%%
lag = 1;
[F Ft Fa] = stft([],data, fs, winLen, stepLen, 400); %prime window size
ph = F./abs(F); %normalize fourier amplitude
df = fs/winLen;
%dph = F(1:end-lag,:,:)./F(1+lag:end,:,:);
%keyboard
for i=1:(200/df)+1
    dph(i,:) = mean((ph(i+1,:,:)./ph(i,:,:)),3);   
end
NFC = abs(dph(1:end-1,:))+abs(dph(2:end,:));