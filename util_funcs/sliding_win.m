function stacked = sliding_win(data, winLen, stepLen)

%rotate
if size(data,1)==1
    data = data';
end

max_iter = ceil((size(data,1)-winLen)/stepLen);
stacked = zeros(winLen, max_iter);
for i=1:max_iter 
    stacked(:,i) = data((1:winLen)+(i-1)*stepLen,:);
end