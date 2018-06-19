function plot_filled(x,y,yb,c)
%function plot_filled(x,y,yb,c)
% utility function to plot mean +/- std
% 
% x: x-coors
% y: y-coors
% y: width around y (std)
% c: color, can be string code (like 'k') or RGB
if isempty(c)
    c = [0 0 0];
end
if isempty(x)
    x = 1:length(y);
end
if size(x,2)==1; x = x'; end
if size(y,2)==1; y = y'; end
if size(yb,2)==1; yb = yb'; end
h = fill([x fliplr(x)],[y+yb fliplr(y-yb)], c, 'edgecolor', 'none');
set(h,'facealpha',.2)
set(h,'edgealpha',.2)
hold on
plot(x,y,'color',c,'linewidth',1)
%plot(x,y, 'o-','color',c, 'markerfacecolor',c, 'markersize',2)
hold off
if x(1)<x(end)
    xlim([x(1) x(end)])
else
    xlim([x(end) x(1)])
end
