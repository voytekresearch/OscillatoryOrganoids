function [burst_labels, b_inds, b_times] = burst_detect(spk_times, min_num, burst_int)
%[burst_labels, b_inds, b_times] = find_bursts(spk_times, min_num, burst_int)
%   spike_times: timestamp of spike train
%   min_num: minimum number of spikes required to be labeled burst (>=2)
%   burst_int: time difference between consecutive spikes to be considered

burst_labels = zeros(size(spk_times));

bursting = 0;
cur_burst = 1;
for i=1:length(spk_times)-min_num+1
    
    %scroll through windows
    cur_inds = i+(1:min_num)-1;
    cur_win = spk_times(cur_inds);
    
    if all(diff(cur_win)<burst_int)
        %group in burst, label by burst number        
        burst_labels(cur_inds)=cur_burst;
        bursting=1;
    else
        %group not bursting, add burst counter        
        if bursting
            cur_burst = cur_burst+1;
            bursting = 0;
        end
    end
    
end
%get interburst interval and count
b_inds = zeros(max(burst_labels),2);
for i=1:max(burst_labels)    
   b_inds(i,1)=find(burst_labels==i,1,'first'); 
   b_inds(i,2)=find(burst_labels==i,1,'last');    
end
b_times=spk_times(b_inds);


% figure
% plot(spk_times,1,'bo')
% hold on
% plot(spk_times,burst_labels>0, 'r')
% hold off