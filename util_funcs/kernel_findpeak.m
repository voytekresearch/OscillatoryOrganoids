function pktimes = kernel_findpeak(x, MPH, MPH_ratio, MPD, toss_thresh_ratio)
% returns pktimes, which has 3 columns: 
% [first crossing, peak, last crossing]
%% iterative peak sifting
% 1. find all peaks (above some optional minimum value and dist)
% 2. sort peak indices by peak amplitude
% while (a peak was eliminated on last iter)
%   for all peak indices
%       find all indices flanking the first peak that are above certain threshold
%       find all detected peak indices that fall within those flanking indices
%           if set is not empty
%               throw out peak indices, peak eliminated = true

% whether to plot the spikes
do_viz = 0;

% first, find all peaks satisfying criterion of minpeakheight and
% minpeakdistance
[PK,IND] =findpeaks(x, ...
    'minpeakheight', max(max(x)*MPH_ratio, MPH), ...
    'minpeakdistance',MPD);

% sort peaks descending
[PK_sorted, sorted_ind] = sort(PK, 'descend');
IND_sorted = IND(sorted_ind);

% set sifting threshold
thresh_val = max(PK_sorted)*toss_thresh_ratio;

% sifting
peaks_eliminated = 1;
did_elim = 0;
while peaks_eliminated
    peaks_eliminated = 0; % reset flag
    cur_ind = 1;
    while cur_ind<=length(PK_sorted)
        % loop through all peak inds, use while loop because number of
        % peaks dynamically decreases       
        
        % find the last index before the peak that met threshold
        i1 = find(x(1:IND_sorted(cur_ind))<thresh_val, 1, 'last');
        % find the first index after the peak that met threshold
        i2 = find(x(IND_sorted(cur_ind)+1:end)<thresh_val, 1, 'first')+IND_sorted(cur_ind);
        
        % if either i1 or i2 is empty, need to resolve conflict
        % set them as beginning and end of recording        
        if isempty(i1); i1=1; end
        if isempty(i2); i2=length(x); end
        % OR just toss that peak
        
        % peak to toss: indices that fell within (and is not) the current one
        toss_inds = find(IND_sorted>i1 & IND_sorted<i2 & IND_sorted~=IND_sorted(cur_ind));                        
        
        if ~isempty(toss_inds)
            if do_viz
                viz_sifting(x, PK_sorted, IND_sorted, thresh_val, toss_inds, i1, i2)
            end
            
            % set the tossed inds as 0
            PK_sorted(toss_inds)=[];
            IND_sorted(toss_inds)=[];
            peaks_eliminated = 1;
            did_elim = 1;
        end
        cur_ind=cur_ind+1;
    end
end

% done sifting, collect peak times
if ~isempty(IND_sorted)
    pktimes = zeros(length(IND_sorted),3);
    % peak times
    pktimes(:,2) = IND_sorted;
    for cur_ind=1:length(IND_sorted)
         % find the last index before the peak that met threshold
        i1 = find(x(1:IND_sorted(cur_ind))<thresh_val, 1, 'last');
        % find the first index after the peak that met threshold
        i2 = find(x(IND_sorted(cur_ind)+1:end)<thresh_val, 1, 'first')+IND_sorted(cur_ind);       
        % if either i1 or i2 is empty, need to resolve conflict
        % set them as beginning and end of recording 
        if isempty(i1); i1=1; end
        if isempty(i2); i2=length(x); end
        pktimes(cur_ind,1) = i1;
        pktimes(cur_ind,3) = i2;
    end
    % sort in time
    pktimes = sort(pktimes,1);
else
    pktimes = [];
end

if did_elim & do_viz
    viz_sifting(x, PK_sorted, IND_sorted, thresh_val, toss_inds, pktimes(:,1), pktimes(:,3))    
    title('DONE')
    pause
end

function viz_sifting(x, PK_sorted, IND_sorted, thresh_val, toss_inds, i1, i2)
% visualize the sifting
plot(x, 'k-')
hold on
plot(IND_sorted, PK_sorted, 'or')
plot(IND_sorted(toss_inds), PK_sorted(toss_inds), 'xb')
plot([i1, i1], ylim, 'k-')
plot([i2, i2], ylim, 'k-')
plot(xlim, [thresh_val thresh_val], 'k--')
hold off
%pause

%% 
% function [kernels, kerTimes] = kernel_findpeak(pop_spk, ker_win, MPH_ratio, MPD)
% % [kernels, kerTimes] = kernel_findpeak(pop_spk, ker_win, MPH_ratio, MPD)
% %baked in parameter to eliminate "peaks" that are not at least 1 spk tall
% min_event_spike = 1;
% 
% %find peaks
% [PK,IND] =findpeaks(pop_spk, ...
%     'minpeakheight', max(max(pop_spk)*MPH_ratio, min_event_spike), ...
%     'minpeakdistance',MPD);
% 
% %collect kernels and output kernel times
% kernels = collect_spikes(pop_spk,[],IND,ker_win);
% kerTimes = IND;



