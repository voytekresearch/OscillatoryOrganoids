function [cf, bw, pw] = inst_cfpw(data, fs, osc_band, winLen, stepLen)
%[cf bw pw] = inst_cfpw(data, fs, osc_band, winLen, stepLen)
% inputs:
%   data - data, will automatically rotate 2D matrix to be compatible
%   fs - sampling frequency
%   osc_band - [lowcut highcut] of band pass filter
%   winLen - window length (samples) to do the median averaging, 0 for raw
%   stepLen - step length to move the window for median averaging
% outputs:
%   cf - center frequency; instantaneous center freq if no winLen specified
%   bw - bandwidth, [] if no winLen specified
%   pw - instantaneous power

% rotate data to column vector
if size(data,2)>size(data,1)
    data=data';
end
%size(data)
disp('Filtering...')
osc = eegfilt(data',fs, osc_band(1), osc_band(2))';
IF = (diff(unwrap(angle(hilbert(osc))))*fs/2/pi);
PW = log10(abs((hilbert(osc))).^2);

% not doing median stacking
if winLen ==0
    cf = IF;
    bw = [];
    pw = PW(1:end-1);
    return
end

num_chan = size(data,2);
cf = zeros(ceil((length(IF)-winLen)/stepLen), num_chan);
bw = cf;
pw = cf;
disp('Stacking...')
for i=1:num_chan    
    cf_stack = sliding_win(IF(:,i),winLen,stepLen);
    cf(:,i) = median(cf_stack);
    bw(:,i) = iqr(cf_stack);
    pw(:,i) = median(sliding_win(PW(2:end,i),winLen,stepLen));
end
end

