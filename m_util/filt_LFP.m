function LFP = filt_LFP(data,fs)

LFP_cut = 300;

[b,a] = butter(10, 2*LFP_cut/fs, 'low');
LFP = filtfilt(b,a,data);

end