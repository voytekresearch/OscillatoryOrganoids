function Sp = filt_spike(data,fs)

spike_freq = [300, 2800];

[b,a] = butter(10, 2*spike_freq/fs, 'bandpass');
Sp = filtfilt(b,a,data);

end