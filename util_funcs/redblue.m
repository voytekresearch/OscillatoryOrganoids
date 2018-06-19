function C = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

colors = get(gca,'colororder');
BB = colors(1,:);
RR = colors(2,:);
m1 = m/2;
C = ones(m,3);
for i=1:m/2
    C(i,:) = ([1,1,1]-RR)*(i-m/2)/(m/2)+[1,1,1];
    C(m-i+1,:) = ([1,1,1]-BB)*(i-m/2)/(m/2)+[1,1,1];
end
