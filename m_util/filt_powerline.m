function data = filt_powerline(data, fs)

%line_freq = [60, 180, 300, 454, 540, 1260, 1363];
line_freq = [58 62; 178 182; 298 302; 538 542; 1256 1262; 1358 1366; 1616 1622; 2263 2277];
for i=1:size(line_freq,1)
    stopband = [line_freq(i,1),line_freq(i,2)];
    [b, a]=butter(3,2*stopband/fs, 'stop');
    data = filtfilt(b,a,data);    
end

end