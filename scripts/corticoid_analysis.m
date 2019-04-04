% network spike aggregate analysis over all recordings
%cd /Users/rgao/Documents/Data/Muotri/Pri_Corticoids/ORGANOIDS
%cd('/Volumes/My Passport for Mac/Dish/CTC')
ker_win = [-500,2500];
smo_len = 100;
LFP_smo = {};
%gaussian smoothing window to avoid harmonics in PSD
smo_mask = gausswin(smo_len)./sum(gausswin(smo_len));
fs = 12500;
fs_ds = 1000;
min_event_spike = 1; % peak is at least 1spk tall to be considered event
%recs = [1:25 27:33 36];
wells = 1:12;

%% collect aggregate data
% features:
%   - event driven population spiking, LFP kernel
%   - aggregate PSDs

%F = dir('CTC_*');
%reorder folders because 2017 recordings are in front
%F = [F(15:end);F(1:14)];

for f=1:length(F)
    cd(F(f).name)
    disp(F(f).name)
    load LFP_Sp.mat LFP spikes spike_cnt t_s t_ds
    
    %recording time
    T{f} = t_ds(end);
    
    %binarize spikes & compute network spiking
    disp('Bin Spikes..')
    bsp = binarize_spikes(ceil(t_ds(end)), fs,spikes,fs_ds);
    nws = squeeze(sum(bsp,2))'; nws = nws(1:length(t_ds),:);
    
    disp('Kernel Calculations..')
    for well=wells
        if isempty(LFP{well})
            %skip well if there's no LFP
            disp(sprintf('Well %i skipped, no data.', well))
            continue
        end
        %smooth network spike vector
        nws_smo{f}(:,well) = conv(nws(:,well), smo_mask,'same');
        
        %find network bursts based on peaks in smoothed well spikes 
        [PK,IND] =findpeaks(nws_smo{f}(:,well), 'minpeakheight', max(max(nws_smo{f}(:,well))*0.75,min_event_spike), 'minpeakdistance',fs_ds);
        
        % well 10 in some recordings have missing channels, normalize to 64
        adj = sum(1-all(LFP{well}==0))/64;
        if adj~=1
            disp(sprintf('Well %i kernel adjusted.',well))
        end
        kernels{f}{well} = collect_spikes(nws_smo{f}(:,well),[],IND,ker_win)/adj;
        
        %record spike stamps
        kerTimes{f}{well} = IND./fs_ds;
        
        %collect LFP kernel at those time points
        for chan = 1:64     
            if all(LFP{well}(:,chan)==0)
                %LFP being all zero screws up the event detection, so just
                %manually fill in all 0s for the event
                %this is a patchy solution but will do for now
                fill_shape = size(kernels{f}{well});
                LFP_ker_smo{f}{well}(:,:,chan) = zeros(fill_shape(1:2));
                LFP_ker{f}{well}(:,:,chan) = zeros(fill_shape(1:2));
            else
                %make smoothed LFP
                LFP_smo = conv(LFP{well}(:,chan), smo_mask,'same');
                % get LFP kernel
                LFP_ker_smo{f}{well}(:,:,chan) = collect_spikes(LFP_smo,[],IND,ker_win);
                LFP_ker{f}{well}(:,:,chan) = collect_spikes(LFP{well}(:,chan),[],IND,ker_win);
            end
        end
        
        %calculate PSD
        PSDw{f}{well} = pwelch(LFP{well},fs_ds*2,fs_ds*1.5,fs_ds*2,fs_ds);
        %PSDm{f}{well} = mPSD(LFP{well}, fs_ds, fs_ds*2, fs_ds/2, fs_ds/2);
    end    
    cd ..    
end    
%save aggregate.mat ker_win T nws_smo kernels kerTimes LFP_ker LFP_ker_smo PSDw PSDm -v7.3
save aggregate.mat ker_win T nws_smo kernels kerTimes LFP_ker LFP_ker_smo PSDw -v7.3

%% secondary processing on event kernel
cd /Users/rgao/Documents/Data/Muotri/Pri_Corticoids
load aggregate.mat
%% event stats & LFP self similarity
wells = 5:12;
for well = wells
    for rec = 1:length(kernels)
        n_events = size(kernels{rec}{well},2);
        
        % mean and std of inter event latency
        ker_latM(rec,well-4) = mean(diff(kerTimes{rec}{well}));
        ker_latS(rec,well-4) = std(diff(kerTimes{rec}{well}));
        
        % event count
        ker_cnt(rec, well-4) = n_events;
        
        % mean self correlation
        if n_events<2
            %only 1 event, don't compute
            self_cor(rec,well-4) = NaN;
            self_cor_LFP(rec,well-4,:) = NaN;
        else
            %more than 1 event, compute pair-wise correlation and average
            self_cor(rec,well-4) = sum(sum(triu(corr(kernels{rec}{well}),1)))/nchoosek(n_events,2);
            
            % LFP kernel correlation
            for chan = 1:64
                self_cor_LFP(rec,well-4, chan) = sum(sum(triu(corr(LFP_ker_smo{rec}{well}(:,:,chan)),1)))/nchoosek(n_events,2);                
            end
        end
        
        % total spiking under the curve
        ker_totalM(rec, well-4) = mean(sum(kernels{rec}{well}));
        
        % initial peak amplitude
        peakAmp = max(kernels{rec}{well});
        ker_peakM(rec, well-4) = mean(peakAmp);
        ker_peakS(rec, well-4) = std(peakAmp);                      
    end    
end

%% event subpeaks analyses & LFP spatial corr
wells = 5:12;
subpks = {};
subpksT = {};
FWHM = {};
nsubpks = zeros(length(kernels),8);
LFP_spatialcorr = zeros(length(kernels),8);
sp2take=5;
p2P = zeros(length(kernels),8,sp2take);
p2T = zeros(length(kernels),8,sp2take);
for well = wells
    for rec = 1:length(kernels)

        n_events = size(kernels{rec}{well},2);        
        % subpeak processing        
        peakAmp = max(kernels{rec}{well});
        
        %peak finding        
        subpks{rec,well-4} = nan(10,n_events);
        subpksT{rec,well-4} = nan(10,n_events);
        nsubs = zeros(1,n_events);
        FWHM{rec,well-4} = zeros(1,n_events);
        
        temp_corr = zeros(64,64,n_events);
        for event = 1:n_events                        
            %sub-peak finding
            %normalize by first peak            
            [pks, pkts] = findpeaks(kernels{rec}{well}(:,event)/peakAmp(event),...
                'minpeakheight', 1/4, ...
                'minpeakwidth', 50, ...
                'npeaks',10, ...
                'minpeakdistance',200, ...
                'minpeakprominence',1/peakAmp(event));            
             
            subpks{rec,well-4}(1:length(pks), event) = pks; %subpeak height
            subpksT{rec,well-4}(1:length(pks), event) = pkts; %subpeak time
            nsubs(event) = length(pks); %number of subpeaks
            
            %FWHM of initial peak            
            if ~isempty(pkts)
                onset = find(kernels{rec}{well}(1:500,event)/peakAmp(event)<0.5, 1, 'last');
                offset = pkts(1)+find(kernels{rec}{well}(pkts(1)+1:end,event)/peakAmp(event)<0.5, 1, 'first');
                FWHM{rec,well-4}(event) = offset-onset;                        
            end
            
            % LFP spatial corr
            temp_corr(:,:,event) = abs(corr(squeeze(LFP_ker_smo{rec}{well}(:,event,:))));
        end
        
        nsubpks(rec,well-4) = mean(nsubs);
        for pk = 1:sp2take
            p2P(rec, well-4, pk) = nanmedian(subpks{rec,well-4}(pk,:));
            p2T(rec, well-4, pk) = nanmedian(subpksT{rec,well-4}(pk,:));                        
        end
        
        %spatial correlation
        LFP_spatialcorr(rec, well-4) = (sum(sum(mean(temp_corr,3)))-64)/(64^2-64);                   
    end
end

%% PSD analysis of LFP kernel
dF = 0.5;
osc_freq = 2.5:dF:4.5;
osc_idx = osc_freq/dF+1;
norm_freq = (0.5:dF:10);
%norm_freq = [2, 5];
norm_idx= norm_freq/dF+1;
fit_freq = 0.5:dF:20;
slope_range = {1:dF:10, 10:dF:30, 30:dF:50};
warning('off')
osc_power = zeros(length(kernels),length(wells),64);
osc_power_foof = zeros(length(kernels),length(wells),64);
slopes = zeros(length(kernels),length(wells),64,length(slope_range));

for well = wells
    disp(well)
    for rec = 1:length(kernels)
        
        osc_power(rec, well-4,:) = mean(log10(PSDw{rec}{well}(osc_idx,:)))-mean(log10(PSDw{rec}{well}(norm_idx,:)));       
        
        sfit = robfitPSD(PSDw{rec}{well},fit_freq,dF);
        osc_power_foof(rec, well-4,:) = mean(log10(PSDw{rec}{well}(osc_idx,:))) - (mean(log10(osc_freq)'*sfit(:,2)')+sfit(:,1)');
        
%         for slp = 1:length(slope_range)
%             sfit = robfitPSD(PSDw{rec}{well},slope_range{slp},dF);
%             slopes(rec,well-4,:,slp) = sfit(:,2);
%         end
    end
end

%% per well regression
% fit order 2 polynomial on 6 features, or a subset of the 6
load names.mat
%exclude recordings not done after exchanging media, or pharmacology
%exclusion = [2 5 9 12 15 25 26 34 35 39];
exclusion = [5 9 12 15 25 26 34 35 36 40];
recs = setdiff(1:length(dates),exclusion);
dayVec = zeros(1,length(recs));
for day = 1:length(recs)
    %parse numerical date
    date = dates(recs(day)).name(5:end);    
    dayVec(day) = datenum(str2num(date(5:6)), str2num(date(1:2)), str2num(date(3:4)));
    disp(date)
end
recDays = dayVec-dayVec(1)+1;
%%
features = (cat(3, ker_latM, ker_latS./ker_latM, nsubpks, ker_peakM, self_cor, mean(osc_power_foof,3)));
features = features(:,:,[1 4 6]);
features = cat(3,features,features.^2);
SSq = zeros(8, size(features,3)+1);
PV = zeros(8, size(features,3)+1);
FITNESS = zeros(8,2);
for well = 1:8
    trainX = squeeze(features(recs,well,:));   
    mdl = fitlm(trainX, recDays);    
    FITNESS(well,:) = [mdl.Rsquared.Adjusted, sqrt(mdl.MSE)];
    A = anova(mdl);
    PV(well,:) = A.pValue;
    SSq(well,:) = A.MeanSq;
end
disp('---Within-Well Regression---')
disp(sprintf('Mean R^2: %f; Mean RMSE: %f',mean(FITNESS)));
%% cross-well classifier analysis
%features = (cat(3, ker_latM, ker_latS./ker_latM, nsubpks, ker_peakM, self_cor, mean(osc_power_foof,3)));
%features = cat(3,features(:,:,[1 4 6]));
%features = cat(3,features,features.^2);
N_train = 7; %6:2 split
N_test = 8-N_train;
Ycat = repmat(recDays,1,N_train);
groups = nchoosek(1:8,N_train);
train_fitness = zeros(size(groups,1),2);
test_fitness = zeros(size(groups,1),2);
SSq = zeros(size(groups,1),size(features,3)+1);
for comb = 1:size(groups,1)
    trainX = reshape(features(recs,groups(comb,:),:),length(recs)*N_train,[]);
    testX = reshape(features(recs,setdiff(1:8,groups(comb,:)),:),length(recs)*N_test,[]);
    mdl = fitlm(trainX, Ycat);        
    Ypred = predict(mdl,testX);
    good_idx = find(1-isnan(Ypred));    
    train_fitness(comb,:) = [mdl.Rsquared.Adjusted, sqrt(mdl.MSE)];
    test_fitness(comb,:) = [corr(Ycat(good_idx)', Ypred(good_idx)).^2, sqrt(mean((Ycat(good_idx)'-Ypred(good_idx)).^2))];    
    AN = anova(mdl);
    SSq(comb,:) = AN.MeanSq;    
end
disp('---Cross-Well Regression---')
disp(sprintf('Mean R^2: %f; Mean RMSE: %f',mean(test_fitness)));

