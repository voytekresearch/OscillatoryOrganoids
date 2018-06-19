function dataF = butterpass(data,fs, passband, order)
%dataF = butterpass(data,fs, passband, order)
%bandpass filter data matrix using butterworth filter
%outputs:
%   dataF: output data matrix
%inputs:
%   data: data matrix, [sample x channel]
%   fs: sampling rate
%   passband: [lowcut highcut]
%   (optional)order: filter order, default = 2

if nargin==3
    order=2;
end

[b,a] = butter(order, 2*passband/fs, 'bandpass');
dataF = filtfilt(b,a,data);

end