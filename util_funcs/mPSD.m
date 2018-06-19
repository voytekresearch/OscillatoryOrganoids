function [PSD, Fa] = mPSD(data, Fs, winLen, stepLen, end_freq)
%function mPSD(data, Fs, winLen, stepLen, end_freq)
% computes median of stft, colloquially referred to as median welch.
% Basically it does the exact same thing as welch's method, but takes
% median instead of the mean, because median is robust to huge outliers and
% spectral data is non-uniformly distributed
%
% data: time series
% Fs: sampling frequency of the data
% winLen: window length in number of samples (use same number as Fs)
% stepLen: step length in number of samples (use about 1/50 as Fs (or 10 samples)
% end_freq (optional): cut of frequency of stft, default is fs/2
% PSD: returns median PSD
% Fa: returns frequency axis

if size(data,1)<size(data,2)
    data = data';
end

[F, Ft, Fa] = stft([], data, Fs, winLen, stepLen, end_freq);
PSD = median(abs(F).^2,3);