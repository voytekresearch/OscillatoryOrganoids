%% Negraes, Gao, Trujillo et al. 2018
%% master script for processing organoid spikes and LFPs 

%% loading data
cd /Users/rdgao/Documents/data/Muotri/Pri_Corticoids
%load('/Users/rdgao/Documents/data/Muotri/InfantEEGFeatures/preterm_features.mat')
load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/aggregate.mat', 'nws_smo', 'T')
load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/names.mat')

%% 
% grab date information and set parameters
wells = 5:12;
%exclude recordings not done after exchanging media, or pharmacology
exclusion = [5 9 12 15 25 26 34 35 36 40];
recs = setdiff(1:length(dates),exclusion);
dayVec = zeros(1,length(recs));
for day = 1:length(recs)
    %parse numerical date
    date = dates(recs(day)).name(5:end);
    disp(date)
    dayVec(day) = datenum(str2num(date(5:6)), str2num(date(1:2)), str2num(date(3:4)));
end
recDays = dayVec-dayVec(1)+1;
dayVec = recDays/7+7; %use weeks instead

%% spike pre-processing
% from raw spiketime data: 
%   - binarize spikes and get network spiking


%% peak finding
% from network spikes (nws_smo), find peaks
fs = 1000;
MPH = 1;
MPH_ratio = 0.8; 
MPD = 200; 
toss_thresh_ratio = 0.001;
% call script
warning('off')
disp('Peak finding & Event Features...')
ctc_findpeaks % sub-script
save('event_timing.mat', 'pktimes', 'IEI*', 'ED*', 'EPH')

%% grab features for baby eeg classifier analysis
ctc_eegfeatures % sub-script
save('ctc_EMAfeatures.mat', 'dayVec', 'ctc_EMAfeatures')

%% grab extra organoid features for self-prediction