function spikes = collect_spikes(data, peaks, loc, winlen)
%spikes = collect_spikes(data, peaks, loc, winlen)
%   data: single channel data
%   peaks: peak value, only applies for spike info
%   loc: trigger index
%   winlen: [W] for -W:W, or [A, B] for A:B


%sort peaks to descending
if ~isempty(peaks)
    [p, ind] = sort(peaks);
    loc = loc(ind);
end
if length(winlen)>1
    %non-symmetric window
    window = winlen(1):winlen(2);
else
    %symmetric window
    window = -winlen:winlen;
end
%grab spikes
spikes = zeros(length(loc), length(window));
for i=1:length(loc)  
    try
        spikes(i,:) = data(loc(i)+window);
    catch
        disp(sprintf('Spike window #%i outside of data range', i))
    end
end
%remove empty rows
spikes(all(spikes==0,2),:)=[];
spikes=spikes';