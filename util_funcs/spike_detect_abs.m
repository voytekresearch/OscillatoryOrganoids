function [spikes, spk_cnt] = spike_detect_abs(data,fs, std_thr)
%function spikes = spike_detect_abs(data,fs)
% data: time X chan (filtered)
% std_thr: std threshold, usually 5 or 5.5 (peak and valley)
mpd = ceil(0.002*fs); %min peak distance 2ms

[len, numchan] = size(data);
spikes = cell(numchan,1);
spk_cnt = zeros(numchan,1);
for chan = 1:numchan
    %robust estimation by R.Q.Quiroga
    thresh = median(abs(data(:,chan))/0.6745)*std_thr;              
    if max(abs(data(:,chan)))>thresh
        [p, spikes{chan,1}] = findpeaks(abs(data(:,chan)), 'minpeakdistance', mpd, 'minpeakheight', thresh);
        spk_cnt(chan) = length(p);        
    else        
        disp(sprintf('No spikes on channel: %i',chan));
    end
end
end