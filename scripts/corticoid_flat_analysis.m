fig_folder =  '~/Dropbox/Research/Reports/Muotri/CorticoidFigs/';
cd /Users/rdgao/Documents/data/Muotri/Pri_Corticoids
load aggregate.mat
figure
colors = get(gca,'colororder');
close all
%% event stats & LFP self similarity
wells = 5:12;
load names.mat
%exclude recordings not done after exchanging media, or pharmacology days
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

%% interevent interval - SI figure
recDays_flat = cell(1,length(recs));
ker_lat_flat = cell(1,length(recs));
for r_idx = 1:length(recs)
    rec = recs(r_idx);
    for well=wells
        
        % hack to compare the results for peak finding
        if ~isempty(pktimes{rec}{well})
            kerTimes{rec}{well} = pktimes{rec}{well}(:,2)/1000;
        else
            kerTimes{rec}{well} = [];
        end
        
        n_events = length(kerTimes{rec}{well});        
        recDays_flat{r_idx} = [recDays_flat{r_idx} recDays(r_idx).*ones(1,n_events-1)];        
        %inter event latency
        ker_lat_flat{r_idx} = [ker_lat_flat{r_idx} diff(kerTimes{rec}{well})'];        
        
        ker_latM(r_idx,well-4) = mean(diff(kerTimes{rec}{well}));        
        ker_latS(r_idx,well-4) = std(diff(kerTimes{rec}{well}));
    end
end
mdl = LinearModel.fit([recDays_flat{:}], [ker_lat_flat{:}]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
%figure
hold on
X = dayVec;
Y = cellfun(@mean,ker_lat_flat);
Y_s = cellfun(@std,ker_lat_flat)./sqrt(cellfun(@length,ker_lat_flat));
plot_filled(X,Y,Y_s, 'r')
xlabel('Weeks')
set(gca,'ytick',0:50:150)
ylabel('Inter-Event Interval (s)')
xlim([8 40])
set(gca,'xtick',10:10:40)
text(10,120, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
%nice_figure(gcf, [fig_folder 'SI_IEI'],[3 3])

%% analysis for first peak width, height, and sum
recDays_flat = cell(1,length(recs));
ker_total_flat = cell(1,length(recs)); %total spiking
ker_pk_flat = cell(1,length(recs)); %peak 1 amplitude
ker_width_flat = cell(1,length(recs)); %peak 1 FWHM
for r_idx = 1:length(recs)
    rec = recs(r_idx);
    for well = wells
        n_events = size(kernels{rec}{well},2);
        recDays_flat{r_idx} = [recDays_flat{r_idx} recDays(r_idx).*ones(1,n_events)];        
        % 1e: total spikes under event
        ker_total_flat{r_idx} = [ker_total_flat{r_idx} sum(kernels{rec}{well})];
        
        % 1f top: first peak amplitude
        peakAmp = kernels{rec}{well}(501,:);
        ker_pk_flat{r_idx} = [ker_pk_flat{r_idx} peakAmp];
        ker_peakM(r_idx,well-4) = mean(peakAmp);        
        
        % 1f bottom: first peak width
        for event = 1:n_events
            onset = find(kernels{rec}{well}(1:500,event)/peakAmp(event)<0.5, 1, 'last');
            offset = 501+find(kernels{rec}{well}(501+1:end,event)/peakAmp(event)<0.5, 1, 'first');
            ker_width_flat{r_idx} = [ker_width_flat{r_idx} offset-onset];            
        end                
    end
end
%% IEI CV - figure 3
X = dayVec;
Y_flat = reshape(ker_latS./ker_latM,1,[]);
Y = nanmean(ker_latS./ker_latM,2);
Y_s = nanstd(ker_latS./ker_latM,1,2)./sqrt(size(ker_latM,2));
mdl = LinearModel.fit(repmat(recDays,1,8), Y_flat);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
figure
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
ylabel('Inter-event Interval CV')
text(10,0.8, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
xlim([8 40])
set(gca,'xtick',10:10:40)
nice_figure(gcf, [fig_folder '3_IEI_CV'],[3 3])


%% 1st peak amplitude - SI figure
figure
%P1 amplitude
Y = cellfun(@mean,ker_pk_flat);
Y_s = cellfun(@std,ker_pk_flat)./sqrt(cellfun(@length,ker_pk_flat));
mdl = LinearModel.fit([recDays_flat{:}], [ker_pk_flat{:}], [0;1;2;]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
ylabel('Peak 1 Amplitude')
text(10,23, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
xlim([8 40])
set(gca,'xtick',10:10:40)
nice_figure(gcf, [fig_folder 'SI_p1amp'],[3 3])

%% subpeak analysis
subpks = {};
subpksT = {};
sp2take=5;
event_end_thresh = 0.01;
recDays_flat = cell(1,length(recs));
nsubs_flat = cell(1,length(recs)); %number of subpeaks
subpk_amp_flat = cell(1,length(recs)); %subpeak amplitude
subpk_delay_flat = cell(1,length(recs)); %subpeak time
for r_idx = 1:length(recs)
    rec = recs(r_idx);
    for well = wells
        n_events = size(kernels{rec}{well},2);
        recDays_flat{r_idx} = [recDays_flat{r_idx} recDays(r_idx).*ones(1,n_events)];
        peakAmp = kernels{rec}{well}(501,:);
        
        %peak finding        
        subpks{rec,well-4} = nan(10,n_events);
        subpksT{rec,well-4} = nan(10,n_events);
        
        %temp_corr = zeros(64,64,n_events);
        for event = 1:n_events                        
            %sub-peak finding
            %normalize by first peak            
            [pks, pkts] = findpeaks(kernels{rec}{well}(:,event)/peakAmp(event),...
                'minpeakheight', 1/4, ...
                'minpeakwidth', 50, ...
                'npeaks',10, ...
                'minpeakdistance',200, ...
                'minpeakprominence',1/peakAmp(event));
            
            %check if event drops below threshold between subpeaks
            %if so, delete subsequent subpeaks
            for i = 1:length(pkts)-1
               if min(kernels{rec}{well}(pkts(i):pkts(i+1),event))<event_end_thresh                   
                   pkts(i+1:end)=[];
                   pks(i+1:end)=[];                   
                  break                  
               end
            end
             
            subpks{rec,well-4}(1:length(pks), event) = pks; %subpeak height
            subpksT{rec,well-4}(1:length(pks), event) = pkts; %subpeak time                        
            nsubs_flat{r_idx} = [nsubs_flat{r_idx} length(pks)];                        
        end        
        subpk_amp_flat{r_idx} = [subpk_amp_flat{r_idx} subpks{rec,well-4}(1:sp2take,:)];        
        subpk_delay_flat{r_idx} = [subpk_delay_flat{r_idx} subpksT{rec,well-4}(1:sp2take,:)];
    end
end
%% number of subpeaks - figure 3
figure
Y = cellfun(@mean,nsubs_flat);
Y_s = cellfun(@std,nsubs_flat)./sqrt(cellfun(@length,nsubs_flat));
mdl = LinearModel.fit([recDays_flat{:}], [nsubs_flat{:}], [0;1;2;]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
plot_filled(X,Y,Y_s, 'k')
set(gca,'ytick',0:1:6)
xlabel('Weeks')
ylabel('# of Subpeaks')
text(10,4.5, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
xlim([8 40])
set(gca,'xtick',10:10:40)
nice_figure(gcf, [fig_folder '3_Nsubpks'],[3 3])

%% subpeak amplitude - figure 3
X = dayVec;
Y = zeros(sp2take,length(subpk_amp_flat));
Y_s = zeros(sp2take,length(subpk_amp_flat));
for idx = 1:length(subpk_amp_flat)
   Y(:,idx) = nanmean(subpk_amp_flat{idx},2);
   Y_s(:,idx) = nanstd(subpk_amp_flat{idx},1,2)./sqrt(sum(1-isnan(subpk_amp_flat{idx}),2));
end
plot(X,Y(2:4,:))
xlim([8 40])
set(gca,'xtick',10:10:40)
ylim([0.2 0.8])
ylabel('Subpeak Amp. (Rel. to PK1)')
legend({'PK2', 'PK3', 'PK4'}, 'location', 'northwest')
xlabel('Weeks')
nice_figure(gcf, [fig_folder '3_subpk_amp'],[3 3])
%% subpeak delay - SI figure
figure
Y = zeros(sp2take-1,length(subpk_delay_flat));
for idx = 1:length(subpk_delay_flat)
   Y(:,idx) = nanmean(diff(subpk_delay_flat{idx}),2);   
end
plot(X,Y(1:3,:)')
xlim([8 40])
set(gca,'xtick',10:10:40)
ylim([200 800])
xlabel('Weeks')
ylabel('Subpeak Latency (ms)')
legend({'P1-P2', 'P2-P3', 'P3-P4'}, 'location', 'northeast')
nice_figure(gcf, [fig_folder 'SI_subpklat'],[3 3])

%% spatial & temporal correlation analysis
recDays_flat = cell(1,length(recs));
spk_corr_flat = cell(1,length(recs));
lfp_meancorr_flat = cell(1,length(recs));
lfp_maxcorr_flat = cell(1,length(recs));
lfpsp_corr_flat = cell(1,length(recs));
for r_idx = 1:length(recs)
    rec = recs(r_idx);
    for well = wells
        n_events = size(kernels{rec}{well},2);
        recDays_flat{r_idx} = [recDays_flat{r_idx} recDays(r_idx).*ones(1,n_events)];        
                
        if n_events<2
            %only 1 event, don't compute
            spk_corr_flat{r_idx} = [spk_corr_flat{r_idx} NaN(1,n_events)];
            lfp_meancorr_flat{r_idx} = [lfp_meancorr_flat{r_idx} NaN(1,n_events)];
            lfp_maxcorr_flat{r_idx} = [lfp_maxcorr_flat{r_idx} NaN(1,n_events)];
        else
            %more than 1 event, compute pair-wise correlation and average
            spk_corr_flat{r_idx} = [spk_corr_flat{r_idx} (sum(corr(kernels{rec}{well}),1)-1)/(n_events-1)];
            
            % LFP kernel correlation
            temp_corr = zeros(n_events,64);
            for chan = 1:64
                temp_corr(:,chan) = (sum(corr(LFP_ker_smo{rec}{well}(:,:,chan)),1)-1)/(n_events-1);                
            end
            lfp_meancorr_flat{r_idx} = [lfp_meancorr_flat{r_idx} mean(temp_corr,2)'];
            lfp_maxcorr_flat{r_idx} = [lfp_maxcorr_flat{r_idx} max(temp_corr,[],2)'];            
        end
        
        for event = 1:n_events
            % LFP spatial corr
            lfpsp_corr_flat{r_idx} = [lfpsp_corr_flat{r_idx} mean((sum(abs(corr(squeeze(LFP_ker_smo{rec}{well}(:,event,:)))),1)-1)/64)];            
        end
    end   
end
%% event temporal correlation - figure 3
figure
Y = cellfun(@nanmean,spk_corr_flat);
Y_s = cellfun(@nanstd,spk_corr_flat)./sqrt(cellfun(@length,spk_corr_flat));
mdl = LinearModel.fit([recDays_flat{:}], [spk_corr_flat{:}], [0;1]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
ylabel('Event Temporal Correlation')
text(10,0.7, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
xlim([8 40])
set(gca,'xtick',10:10:40)
nice_figure(gcf, [fig_folder '3_timecorr'],[3 3])

%% event spatial (LFP) correlation - figure 3
figure
Y = cellfun(@nanmean,lfpsp_corr_flat);
Y_s = cellfun(@nanstd,lfpsp_corr_flat)./sqrt(cellfun(@length,lfpsp_corr_flat));
mdl = LinearModel.fit([recDays_flat{:}], [lfpsp_corr_flat{:}], [0;1;2]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
plot_filled(X,Y,Y_s, colors(1,:))
xlabel('Weeks')
ylabel('Event Spatial (LFP) Correlation')
text(10,0.5, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
xlim([8 40])
set(gca,'xtick',10:10:40)
nice_figure(gcf, [fig_folder '3_spatialcorr'],[3 3])

%% fig 4c: PSD analysis of LFP kernel
dF = 0.5;
osc_freq = 2.5:dF:4.5;
osc_idx = osc_freq/dF+1;
norm_freq = (0.5:dF:10);
norm_idx= norm_freq/dF+1;
fit_freq = 0.5:dF:20;
slope_range = {1:dF:10, 10:dF:30, 30:dF:80};
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
        
        for slp = 1:length(slope_range)
            sfit = robfitPSD(PSDw{rec}{well},slope_range{slp},dF);
            slopes(rec,well-4,:,slp) = sfit(:,2);
        end
    end
end
%% oscillatory power - figure 4
figure
X = dayVec;
Y = mean(mean(osc_power_foof(recs,:,:),3),2);
Y_s = std(mean(osc_power_foof(recs,:,:),3),1,2)/sqrt(8);
mdl = LinearModel.fit(repmat(recDays,1,8), reshape(mean(osc_power_foof(recs,:,:),3),1,[]), [0;1;2]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
plot_filled(X, Y, Y_s, colors(1,:))
xlabel('Weeks')
ylabel('Oscillatory Power')
set(gca,'ytick',0:0.05:0.1)
text(10,0.1, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
set(gca,'xtick',10:10:40)
xlim([8 40])
nice_figure(gcf, [fig_folder '4_oscpow_sq'],[4 3])

%% classifier analysis
%within well
features = (cat(3, ker_latM, ker_peakM, 0.01+mean(osc_power_foof(recs,:,:),3)));
features = cat(3,features,features.^0.5);
SSq = zeros(8, size(features,3)+1);
PV = zeros(8, size(features,3)+1);
FITNESS = zeros(8,2);
for well = 1:8
    trainX = squeeze(features(:,well,:));   
    mdl = fitlm(trainX, recDays);    
    FITNESS(well,:) = [mdl.Rsquared.Adjusted, sqrt(mdl.MSE)];
    A = anova(mdl);
    PV(well,:) = A.pValue;
    SSq(well,:) = A.MeanSq;
end
disp('---Within-Well Regression---')
disp(sprintf('Mean R^2: %f; Mean RMSE: %f',mean(FITNESS)));
disp(sprintf('SEM R^2: %f; SEM RMSE: %f',std(FITNESS)./sqrt(8)));
%% cross-well classifier analysis
N_train = 7; %6:2 split
N_test = 8-N_train;
Ycat = repmat(recDays,1,N_train);
groups = nchoosek(1:8,N_train);
train_fitness = zeros(size(groups,1),2);
test_fitness = zeros(size(groups,1),2);
SSq = zeros(size(groups,1),size(features,3)+1);
for comb = 1:size(groups,1)
    trainX = reshape(features(:,groups(comb,:),:),length(recs)*N_train,[]);
    testX = reshape(features(:,setdiff(1:8,groups(comb,:)),:),length(recs)*N_test,[]);
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
disp(sprintf('SEM R^2: %f; SEM RMSE: %f',std(test_fitness)/sqrt(8)));

%% PAC analysis
load('./organoids_processed/CTC_110416/LFP_Sp.mat')
agg_date = 18;
fs = 1000;
win = -500:2500;
num_bins = 21;
bin_edge = linspace(-pi,pi,num_bins);
%%
for well = 5:12
    disp(well)
    evtimes = kerTimes{agg_date}{well}*fs;    
    GM_EV = {};
    GM_NEV = {};
    
    % get event and non-event indices
    ev_inds = reshape((repmat(evtimes,1,length(win))+repmat(win,length(evtimes),1))',1,[]);
    ev_inds = intersect(ev_inds,1:length(LFP{well}(:,1))); % clip any indices outside of range
    nev_inds = setdiff(1:length(LFP{well}(:,1)),ev_inds);
    
    % filtering & binning
    for chan=1:64
        delta = eegfilt(LFP{well}(:,chan)',fs,0,4)';
        delta_ph = angle(hilbert(delta));
        gamma = eegfilt(LFP{well}(:,chan)',fs,100,200)';
        gamma_amp = abs(hilbert(gamma)).^2;
        
        [count, ev_bins] = histc(delta_ph(ev_inds), bin_edge); % event delta phase
        [count, nev_bins] = histc(delta_ph(nev_inds), bin_edge); % non-event delta phase
        gamma_ev = gamma_amp(ev_inds);
        gamma_nev = gamma_amp(nev_inds);
        
        gm_ev = {};
        gm_nev = {};
        for b = 1:length(bin_edge)-1
            gm_ev{b} = gamma_ev(ev_bins==b);
            gm_nev{b} = gamma_nev(nev_bins==b);
        end
        GM_EV{chan} = gm_ev;
        GM_NEV{chan} = gm_nev;
    end
    
    % realign and calculate modulation index
    max_ph = zeros(64,2);
    MI = zeros(64,2);
    N = num_bins-1;
    for chan = 1:64
        [y, ind] = max(cellfun(@mean,GM_EV{chan}));
        %gmap_ev(:,chan) =  circshift(zscore(cellfun(@mean,GM_EV{chan})),-ind+1,2); %zscore
        gmap_ev(:,chan) =  circshift(cellfun(@mean,GM_EV{chan})/sum(cellfun(@mean,GM_EV{chan})),-ind+1,2); %normed (pdf for MI)
        MI(chan,1) = log2(N)+sum(gmap_ev(:,chan).*log2(gmap_ev(:,chan))); %MI calculation
        max_ph(chan,1) = ind;        
        
        [y, ind] = max(cellfun(@mean,GM_NEV{chan}));
        %gmap_nev(:,chan) =  circshift(zscore(cellfun(@mean,GM_NEV{chan})),-ind+1,2);
        gmap_nev(:,chan) =  circshift(cellfun(@mean,GM_NEV{chan})/sum(cellfun(@mean,GM_NEV{chan})),-ind+1,2);
        MI(chan,2) = log2(N)+sum(gmap_nev(:,chan).*log2(gmap_nev(:,chan)));
        max_ph(chan,2) = ind;                
    end    
    MI_all(:,:,well-4) = MI;
end
%%
[h pv] = ttest(squeeze(mean(MI_all(:,1,:),1)),squeeze(mean(MI_all(:,2,:),1)));
disp(pv)

%%
max_ph = zeros(64,2);
MI = zeros(64,2);
N = num_bins-1;
figure
hold on
for chan = 1:64
    [y, ind] = max(cellfun(@mean,GM_EV{chan}));
    %gmap_ev(:,chan) =  circshift(zscore(cellfun(@mean,GM_EV{chan})),-ind+1,2); %zscore
    gmap_ev(:,chan) =  circshift(cellfun(@mean,GM_EV{chan})/sum(cellfun(@mean,GM_EV{chan})),-ind+1,2); %normed (pdf for MI)
    MI(chan,1) = log2(N)+sum(gmap_ev(:,chan).*log2(gmap_ev(:,chan))); %MI calculation
    max_ph(chan,1) = ind;
    plot(gmap_ev(:,chan), 'color', [0 0 0 0.2]);
    
    [y, ind] = max(cellfun(@mean,GM_NEV{chan}));
    %gmap_nev(:,chan) =  circshift(zscore(cellfun(@mean,GM_NEV{chan})),-ind+1,2);
    gmap_nev(:,chan) =  circshift(cellfun(@mean,GM_NEV{chan})/sum(cellfun(@mean,GM_NEV{chan})),-ind+1,2);
    MI(chan,2) = log2(N)+sum(gmap_nev(:,chan).*log2(gmap_nev(:,chan)));
    max_ph(chan,2) = ind;
    plot(gmap_nev(:,chan), 'color', [1 0 0 0.2]);    
       
end
hold off



%% ------------ Extra stuff --------------------
%
%
%
%
%
%%
% spikes under event
X = dayVec;
Y = cellfun(@mean,ker_total_flat);
Y_s = cellfun(@std,ker_total_flat)./sqrt(cellfun(@length,ker_total_flat));
mdl = LinearModel.fit([recDays_flat{:}], [ker_total_flat{:}], [0;1;2]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
figure
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
set(gca,'ytick', 0:10000:20000)
ylim([0 20000])
set(gca,'yticklabel', num2str((0:3)'))
ylabel('Spikes In Event (x10^5)')
text(10,17000, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
%nice_figure(gcf, [fig_folder '3e'],[3 3])
%%
%P1 fullwidth half maximum
figure
Y = cellfun(@mean,ker_width_flat);
Y_s = cellfun(@std,ker_width_flat)./sqrt(cellfun(@length,ker_width_flat));
mdl = LinearModel.fit([recDays_flat{:}], [ker_width_flat{:}], [0;1;2]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
ylabel('Peak 1 Width')
text(10,330, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
% nice_figure(gcf, [fig_folder '3f'],[3 3])

%% LFP temporal correlation
figure
subplot(2,1,1)
Y = cellfun(@nanmean,lfp_meancorr_flat);
Y_s = cellfun(@nanstd,lfp_meancorr_flat)./sqrt(cellfun(@length,lfp_meancorr_flat));
mdl = LinearModel.fit([recDays_flat{:}], [lfp_meancorr_flat{:}], [0;1;2]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
plot_filled(X,Y,Y_s, colors(1,:))
ylabel('LFP Mean Corr.')
set(gca,'xtick', [])
text(10,0.65, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)

subplot(2,1,2)
Y = cellfun(@nanmean,lfp_maxcorr_flat);
Y_s = cellfun(@nanstd,lfp_maxcorr_flat)./sqrt(cellfun(@length,lfp_maxcorr_flat));
mdl = LinearModel.fit([recDays_flat{:}], [lfp_maxcorr_flat{:}], [0;1]);
r2 = mdl.Rsquared.Ordinary;
pv = coefTest(mdl);
plot_filled(X,Y,Y_s, colors(1,:))
ylim([0.6 1])
ylabel('LFP Max Corr.')
xlabel('Weeks')
text(10,0.7, sprintf(['R^2 = %.3f \np < 10^{-7}'], r2), 'fontsize',10)
%nice_figure(gcf, [fig_folder '3j'],[3 3])

