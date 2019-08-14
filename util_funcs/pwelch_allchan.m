function P = pwelch_allchan(data, fs, overlap)
    Nchan = size(data,2);
    P = zeros(fs/2+1, Nchan);
    for j=1:Nchan
        P(:,j) = pwelch(data(:,j), fs, overlap, fs);
    end
end