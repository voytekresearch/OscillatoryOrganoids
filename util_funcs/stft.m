function [output, output_time, output_f] = stft(timestamp, data, Fs, winLen, stepLen, end_freq)
%function [output, output_time, output_f] = stft(timestamp, data, Fs, winLen, stepLen, end_freq)
% performs short-time windowed fourier transform on data, with a hamming
% window applied. Basically does the same thing as spectrogram.m
% 
% timestamp: leave blank []
% data: time series
% Fs: sampling frequency of the data
% winLen: window length in number of samples (use same number as Fs)
% stepLen: step length in number of samples (use about 1/50 as Fs (or 10 samples)
% end_freq (optional): cut of frequency of stft, default is fs/2

if isempty(timestamp)
   timestamp = linspace(1/Fs,size(data,1)/Fs, size(data,1)); 
end
%winLen=Fs;
max_iter = ceil((size(data,1)-winLen)/stepLen);
%preallocate output 
output_time = zeros(max_iter,1);
output_f = linspace(0,Fs-Fs/winLen, winLen);
if isempty(end_freq)
    end_ind = length(output_f);
else
    end_ind = find(output_f>end_freq,1)-1; %find index
end
output_f=output_f(1:end_ind);
    
output = zeros(end_ind, size(data,2), max_iter);
H = hamming(winLen);
HAM = repmat(H./sqrt(sum(H.^2)),1,size(data,2));

%%%stepping through
for i=1:max_iter
    %dc offset    
    cur_window = data((1:winLen)+(i-1)*stepLen,:)/sqrt(winLen);
    cur_window = (cur_window - repmat(mean(cur_window,1),winLen,1)).*HAM;         
    output_time(i)=timestamp(winLen+(i-1)*stepLen);    
    F = fft(cur_window);
    %use multitaper
    %F = pmtm(cur_window,4,winLen,Fs);
    output(:,:,i) = (F(1:end_ind,:));
    %output(:,:,i) = abs(F(1:end_ind,:));

    %output(:,i)=step_features(cur_window, FREQS);
end
end


function C = step_features(data_window, FREQS)
eeg = data_window;
F = fft(eeg);
%find frequency with highest power
%[temp maxF] = max(log(abs(F(FREQS(2):55,:)))); %cut at 55Hz to avoid powerline
%F3=(maxF+FREQS(2)-1)';

F = F(1:FREQS(3),:); %lowpass
F(60:62,:)=[]; %remove powerline
F=F(FREQS(1):end,:);
spect=log(abs(F));
F1 = reshape(spect,[],1);

D=size(F,2); %dimension (num of channels)
F=F./abs(F); %normalize magnitude

%pairwise phase difference between every channel
C=zeros(size(F,1), sum(1:(D-1)));
ind=1;
for j=1:D
    for k=j+1:D
        C(:,ind)=F(:,j)./F(:,k);
        ind=ind+1;
    end
end
F2 = -squeeze(sum(abs(angle(C)),2))/size(C,2);
C = [F1;F2;F3];
end

function C = phase_cor(data_window, FREQS)
eeg = data_window;
F = fft(eeg);
F=F(1:101,:); %up to 100hz

D=size(F,2); %dimension (num of channels)
F=F./abs(F); %normalize magnitude

%pairwise phase difference between every channel
C=zeros(size(F,1), sum(1:(D-1)));
ind=1;
for j=1:D
    for k=j+1:D
        C(:,ind)=F(:,j)./F(:,k);
        ind=ind+1;
    end
end
F2 = -squeeze(sum(abs(angle(C)),2))/size(C,2);
C = F2;
end