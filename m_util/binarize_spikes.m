function bin_sp = binarize_spikes(t,spk_fs,spikes,bin_rate)
%function bin_sp = binarize_spikes(t,spk_fs,spikes,bin_rate)
%   t: ORIGINAL time vector (just need last time value)
%   spk_fs: sampling rate of spikes
%   spikes: cell of spikes times (well x chan)
%   bin_rate: rate to binarize 
[numwell, numchan, ~] = size(spikes);
len = ceil(t(end)*bin_rate);
bin_sp = zeros(numwell, numchan, len);
for well=1:numwell
    for chan=1:numchan
        if any(spikes{well,chan})
            try
                bin_sp(well,chan,:) = bin_spikes(spikes{well,chan}, [spk_fs, len/bin_rate], bin_rate);
            catch
                keyboard
            end
        end
    end
end
end

