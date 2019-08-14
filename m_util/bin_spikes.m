function bin_sp = bin_spikes(spike_times, t_orig, binned_rate)
%function bin_sp = bin_spikes(spike_times, t_orig, binned_rate)
    % spike_times: vector of spike time stamps in index number
    % t_orig = [spk_rate t_end]: original sampled rate, and end time in s
    % binned_rate: bin rate, usually 1000
    %
    % don't call this for cell data, call binarize_spikes
    
    spk_rate = t_orig(1);
    t_end = t_orig(2);    
    bin_sp = zeros(ceil(t_end*binned_rate),1);
    spk_inds = round(spike_times/spk_rate*binned_rate);
    if length(unique(spk_inds)) == length(spk_inds)
        %if no repeating indices, just index it
        bin_sp(spk_inds(spk_inds>0)) = 1;
    else
        for idx=spk_inds
            bin_sp(idx)= bin_sp(idx)+1;
        end
    end
end