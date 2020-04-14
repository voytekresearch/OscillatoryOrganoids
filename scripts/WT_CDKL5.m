%% CDKL5 vs. WT organoids (see daily17_05_11.m for original code)
set(0,'defaultfigurecolor',[1 1 1])
load('/Users/rdgao/Documents/data/Muotri/CDKL5_WT/aggregate.mat', 'nws_smo', 'ker_win')
load('/Users/rdgao/Documents/data/Muotri/CDKL5_WT/CDKL5_dates.mat')
%% peak finding
% from network spikes (nws_smo), find peaks
fs = 1000;

disp('Peak finding & Event Features...')
warning('off')
fp_params.MPH = 1;
fp_params.MPH_ratio = 0.8; 
fp_params.MPD = 200; 
fp_params.toss_thresh_ratio = 0.001;
%[event_feats, pk_times] = compute_peakfeats(nws_smo, fs, fp_params);

%% some housekeeping stuff
dayVec=1:24;
for i=1:24
    dates{i}=num2str(i);
end
nws_smo{24}(:,12)=0;

%% recompute kernel locations with new code
kernels = cell(24,1);
for r_idx = 1:length(nws_smo)
    for well = 1:12
        %pktimes = kernel_findpeak(nws_smo, fp_params.MPH, fp_params.MPH_ratio, fp_params.MPD, fp_params.toss_thresh_ratio);
        pktimes{r_idx}{well} = kernel_findpeak(nws_smo{r_idx}(:,well), fp_params.MPH, fp_params.MPH_ratio, fp_params.MPD, fp_params.toss_thresh_ratio);
        if ~isempty(pktimes{r_idx}{well})
            kernels{r_idx}{well} = collect_spikes(nws_smo{r_idx}(:,well),[],pktimes{r_idx}{well}(:,2),ker_win);
        else
            kernels{r_idx}{well} = [];
        end
    end
end

%% kernel aggregate stats
nsubs = zeros(length(kernels),12);
peakAmp = zeros(length(kernels),12);
FWHM = zeros(length(kernels),12);
total_spikes = zeros(length(kernels),12);
for day = 1:length(kernels)
    disp(dates{day})
    for well = 1:12
        [nsubs_, peakAmp_,FWHM_, subpks_, subpksT_] = kernel_pkfinding(kernels{day}{well});
        nsubs(day, well) = mean(nsubs_);
        peakAmp(day, well) = mean(peakAmp_);
        FWHM(day, well) = mean(FWHM_);
        total_spikes(day, well) = mean(sum(kernels{day}{well},1));
    end
end
%% load WT data
WT_data = load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/aggregate.mat', 'kernels');
load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/names.mat');
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
recDays = dayVec-dayVec(1);
dayVec = recDays/7+7; %use weeks instead

kernels_WT = WT_data.kernels;
nsubs_WT = zeros(length(kernels_WT),8);
peakAmp_WT = zeros(length(kernels_WT),8);
FWHM_WT = zeros(length(kernels_WT),8);
total_spikes_WT = zeros(length(kernels_WT),8);
for day = 1:length(kernels_WT)
    disp(day)
    for well = 5:12
        [nsubs_, peakAmp_,FWHM_, subpks_, subpksT_] = kernel_pkfinding(kernels_WT{day}{well});
        nsubs_WT(day, well-4) = mean(nsubs_);
        peakAmp_WT(day, well-4) = mean(peakAmp_);
        FWHM_WT(day, well-4) = mean(FWHM_);
        total_spikes_WT(day, well-4) = mean(sum(kernels_WT{day}{well},1));
    end
end


%%
conds = {[1 2 5 6 9 10], [3 4 7 8 11 12]};
plot_feat = {nsubs, nsubs_WT};
figure
plot_filled(CDKL5_dayVec, mean(plot_feat{1}(:,conds{1}),2), std(plot_feat{1}(:,conds{1}),1,2), 'r')
hold on
plot_filled(CDKL5_dayVec, mean(plot_feat{1}(:,conds{2}),2), std(plot_feat{1}(:,conds{2}),1,2), 'k')
plot_filled(dayVec, mean(plot_feat{2}(recs,:),2), std(plot_feat{2}(recs,:),1,2), 'b')
title('Oscillation Subpeaks')
ylabel('# Peaks')
xlabel('Recording')
hold off

figure
hold on

plot(CDKL5_dayVec, plot_feat{1}(:,conds{1}), 'or-')
plot(CDKL5_dayVec, plot_feat{1}(:,conds{2}), 'ok-')
plot(dayVec, plot_feat{2}(recs,:), 'ob-')
hold off
%%

%% mean % of subpeaks
conds = {[1 2 5 6 9 10], [3 4 7 8 11 12]};
figure
hold on
%plot(dayVec,nsubs(:,conds{1}), '.-', 'color', [0 0 0 0.4])
plot(dayVec,mean(nsubs(:,conds{1}),2), 'o-', 'color', [0 0 0 0.4])
%plot_filled(dayVec, mean(nsubs(:,conds{1}),2), std(nsubs(:,conds{1}),2), 'k')

%plot(dayVec,nsubs(:,conds{2}), '.-', 'color', [1 0 0 0.4])
plot(dayVec,mean(nsubs(:,conds{2}),2), 'o-', 'color', [1 0 0 0.4])
hold off
xlabel('DIV')
ylabel('# of subpeaks')
