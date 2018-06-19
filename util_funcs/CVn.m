function [cv_all isi_avg]= CVn(isi, n)
%function [cv_all isi_avg]= CVn(isi, n)
% isi: interspike interval of the spiketrain (i.e. diff(spike_time))
% n: nth order CV, calculates spread over n adjacent spikes
% cv_all: output vector of cv2 between pairs of isi (see HOLT CV2 paper)
if length(isi)<n+1
    cv_all=[];
    isi_avg=[];
    return
end
if size(isi,1)==1
    %turn into column vector
    isi=isi';
end
isi_n = [isi(1:end-(n-1)) isi(n:end)];
isi_avg = mean(isi_n,2);
cv_all = abs(diff(isi_n')')./isi_avg;
if any(isnan(isi_avg))
    keyboard
end