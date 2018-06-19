function ac = autocorr(data, maxlag)
% function ac = autocorr(data, maxlag)
%   computes column-wise autocorrelation of data with specified maxlag
%   data: [time x chan] matrix
%   maxlag: maximum lag to compute autocorrelation

[nsamples, nchan] = size(data);
if nsamples<maxlag
    % data not long enough for requested window
    disp('Not enough samples for requested lag. Exiting.')
    ac = [];
else
    ac = zeros(maxlag*2+1, nchan);
    for chan = 1:nchan
        % compute autocorr using xcor
        ac(:,chan) = xcorr(data(:,chan), maxlag);
    end
end