function spike_raster(t,spikes,well, th_mode)
% function spike_raster(t,spikes,well, th_mode)

%figure
if nargin==2
    %plot binarized data    
    colors = get(gca,'colororder');
    %hold on
    for i=1:size(spikes,1)       
        if any(find(spikes(i,:)))               
            %line([t(find(spikes(i,:)));t(find(spikes(i,:)))], [i-0.2;i+0.2]+20, 'color', colors(2,:), 'linewidth',2);
            %line([t(find(spikes(i,:)));t(find(spikes(i,:)))], [i-0.2;i+0.2], 'color', 'k', 'linewidth',1)
            plot([t(find(spikes(i,:)));t(find(spikes(i,:)))], [i-0.2;i+0.2], 'k.')
        end
    end
    hold off       
else
    %plot time data in cell form
    if nargin<4
        th_mode = 1;
    end    
    hold on
    for i=1:size(spikes,2)
        if ~isempty(spikes{well,i})            
            %line([t(spikes{well,i,th_mode}) t(spikes{well,i,th_mode})]', [i-0.2 i+0.2], 'color', 'k')
            %line([t(spikes{well,i}) t(spikes{well,i})]', [i-0.2 i+0.2], 'color', 'k')
            %plot(t(spikes{well,i,th_mode}), i, 'k.')            
            scatter(t(spikes{well,i,th_mode}), i*ones(1,length(spikes{well,i,th_mode})), 5, 'k.')
        end
    end
    hold off
    ylim([0 size(spikes,2)+1])
end
xlabel('Time')
ylabel('Channel')