% script for generating schematic figures for corticoid paper
fig_folder =  '~/Dropbox/Research/Reports/2018 - CorticoidOscillation/CorticoidFigs/';
cd ~/Documents/Data/Muotri/Pri_Corticoids/
figure
colors = get(gca,'colororder');
close all
%dates: 0816, 0916, 1018, 1118
%%
%% schematic of data - figure 3
% plot example raster
% recording Oct 18, 2016, corresponds to index 15 in data matrix
load('organoids_processed/CTC_101816/LFP_Sp.mat', 'LFP', 'spikes', 'spike_shape', 't_s', 't_ds')
load('aggregate.mat', 'nws_smo', 'kernels', 'LFP_ker_smo')
load names.mat
pkts = [6.319 43.02];
date = find(strcmp({dates.name}, 'CTC_101816'));
well = 12;
chan = 53;
XL = [0 60];
figure
subplot(2,1,1)
spike_raster(t_s,spikes,well)
xlabel('')
set(gca, 'xtick', [])
set(gca, 'ytick',[0,60])
xlim(XL)
title('Spike Raster')
hold on
for pk = 1:2
    xb = pkts(pk) +[-0.5 2.5];
    yb = ylim;
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
    set(fbox, 'facealpha', 0.1)
end
hold off

% plot network spike
subplot(4,1,3)
plot(t_ds, nws_smo{date}(:,well), 'k', 'linewidth', 1)
set(gca, 'xtick', [])
ylabel('Spikes')
xlim(XL)
ylim([0 20])
box off
title('Population Spiking')
hold on
for pk = 1:2
    xb = pkts(pk) +[-0.5 2.5];
    yb = ylim;
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
    set(fbox, 'facealpha', 0.1)
end
hold off

% plot example LFP
subplot(4,1,4)
plot(t_ds, LFP{well}(:,chan), 'linewidth', 1)
xlim(XL)
set(gca, 'ytick', 0)
set(gca, 'yticklabel', '0')
xlabel('Time (s)')
ylabel('Voltage')
box off
title('Local Field Potential (LFP)')
hold on
for pk = 1:2
    xb = pkts(pk) +[-0.5 2.5];
    yb = ylim;
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
    set(fbox, 'facealpha', 0.1)
end
hold off
%nice_figure(gcf, [fig_folder '3_dataschematic'],[6 6])

%% big time trace of pop spiking and LFP - figure 3
% plot example of spiking kernel
well = 12;
chan = 53;
date = find(strcmp({dates.name}, 'CTC_101816'));
t = (-500:2500)/1000;

figure
subplot(2,2,1)
plot(t,kernels{date}{well},'color', [0 0 0 0.5], 'linewidth', 1)
ylabel('Spikes')
xlabel('Time (s)')
set(gca, 'xtick', [0 2])
set(gca, 'ytick', [0 20])
xlim([-0.5, 2.5])
box off
title('Pop. Spiking Kernel')

% plot example of LFP kernel
subplot(2,2,3)
plot(t,LFP_ker_smo{date}{well}(:,:,chan),'color', [colors(1,:) 0.5], 'linewidth', 1)
xlabel('Time (s)')
ylabel('Voltage (uV)')
set(gca, 'xtick', [0 2])
ylim([-2 2]*1e-5)
set(gca, 'ytick', [-2:2:2]*1e-5)
yticklabels(gca,{'-20', '0', '20'})
%set(gca, 'yticklabel', '0')
xlim([-0.5, 2.5])
box off
title('LFP Kernel')

% sample time trace from different wells on same day
% plot kernels from different wells
wells = [5 6 7 8];
subp_idx = [3 4 7 8];
for well = 1:4   
    subplot(4,4,subp_idx(well))
    plot(t,kernels{date}{wells(well)},'color', [0 0 0 0.5], 'linewidth', 1)
    box off
    xlim([-0.5, 2.5])
    ylim([0 25])
    set(gca, 'xtick', [])
    %set(gca, 'ytick', [])
end

% plot LFP kernels from different wells
chans = [9 7 3 35];
subp_idx = [11 12 15 16];
for well = 1:4 
    subplot(4,4,subp_idx(well))
    plot(t,LFP_ker_smo{date}{wells(well)}(:,:,chans(well)),'color', [colors(1,:) 0.5], 'linewidth', 1)
    box off
    xlim([-0.5, 2.5])    
    ylim([-2 2]*1e-5)
    set(gca, 'xtick', [])
    set(gca, 'ytick', [-2, 0, 2]*1e-5)
    yticklabels(gca,{'-20', '0', '20'})
end
nice_figure(gcf, [fig_folder '3_wellexamples'],[6 6])

%% plot PSD and spectrogram
ctc_osc = load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/organoids_processed/CTC_110416/LFP_Sp.mat');
%%
close all
figure
well = 12;
chan = 20;
date = find(strcmp({dates.name}, 'CTC_101816'));
f_axis = 0:0.5:500;
subplot(1,2,1)
plot(f_axis,PSDw{date}{well}(:,chan), 'k', 'linewidth', 1);
xlim([1,10])
xlabel('Frequency (Hz)')
ylabel('Power')

subplot(1,2,2)
loglog(f_axis,PSDw{date}{well}(:,chan), 'k', 'linewidth', 1);
xlim([1,200])
ylim(1e-14*[0.5, 200])
xlabel('Frequency (Hz)')
ylabel('Log10 Power')
nice_figure(gcf, [fig_folder 'supp_PSDs'],[7 3.5])
%%
close all
fs = 1000;
f_lim = 400;
figure
subplot(2,1,1)
plot(ctc_osc.t_ds, 1e6*ctc_osc.LFP{well}(:,chan), 'k')
ylabel('Voltage (uV)')
xlim([0,70])
title('Local Field Potential')

subplot(2,1,2)
[S, F, T] = spectrogram(ctc_osc.LFP{well}(:,chan), fs*2, round(fs*1.75), fs*2, fs);
imagesc(T,F(1:f_lim),log10((abs(S(1:f_lim,:)))./max(max((abs(S(1:f_lim,:)))))))
title('Spectrogram')
set(gca,'ydir', 'normal')
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
xlim([1 60])
set(gca,'xtick',[30, 60])
h = colorbar('Location','eastoutside');
set(h,'ylim',[-3 0])
set(h,'ytick',-3:0)
ylabel(h, 'Log10 Power')
nice_figure(gcf, [fig_folder 'supp_spg'],[8 6])
% set(h,'ycolor','w')
%%

subplot('position',[0.07,y_pos(i),0.62, 0.16])
[S, F, T] = spectrogram(data_all{i}, FS_all{i}*2, round(FS_all{i}*1.75), FS_all{i}*2, FS_all{i});
imagesc(T,F(1:50),(abs(S(1:50,:)))./max(max((abs(S(1:50,:))))))
xlim([T(1) 60])
ylim([0 10])
set(gca,'xtick',[])
box off
title(labels{i})
if i==1
    h = colorbar('Location','east');
    set(h,'ylim',[0 1])
    set(h,'ytick',[0 0.5 1])
    set(h,'ycolor','w')
end

hold on
xb = zoom_win{i};
yb = ylim+[0.2 -0.2];
fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], 'w');
set(fbox, 'facecolor', 'none')
set(fbox, 'edgecolor', 'w')
hold off
set(gca,'ytick',0:10:20)

if i==1
    ylabel('Frequency (Hz)')
end
if i==4
    set(gca,'xtick', 0:20:60)
    xlabel('Time (s)')
end

subplot('position',[0.8,y_pos(i),0.16, 0.16])
loglog(freq_all{i}, PSD_all{i}(1:length(freq_all{i}),chans(i))./(max(PSD_all{i}(1:length(freq_all{i}),chans(i)))), 'linewidth',1)
xlim([0 40])
%ylim([min(PSD_all{i}(1:length(freq_all{i}),chans(i))), max(PSD_all{i}(1:length(freq_all{i}),chans(i)))])
set(gca,'ytick',[0.01 1])
set(gca,'yticklabel', {'0.01' '1'})
set(gca,'xticklabel',[])
box off
if i==1
    ylabel('Power (uV^2/Hz)')
end
if i==4
    set(gca,'xtick', [1 10])
    set(gca,'xticklabel', {'1' '10'})
    xlabel('Frequency (Hz)')
end

%% sample time trace from different days - figure 3
t = (-500:2500)/1000;
days = [3 9 15 24];
well = 5;
chan = 13;
figure
for day=1:4
    subplot(2,4,day)
    plot(t,kernels{recs(days(day))}{well},'color', [0 0 0 0.3], 'linewidth', 1)
    xlim([-0.5, 2.5])
    ylim([0 35])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    title([num2str(day*2), ' Months'])
    box off
end
subplot(2,4,1)
set(gca, 'ytick', [0 35])
ylabel('Pop. Spike (Hz)')

% LFP kernel from different days
for day=1:4
    subplot(2,4,day+4)
    plot(t,LFP_ker_smo{recs(days(day))}{well}(:,:,chan),'color', [colors(1,:) 0.3], 'linewidth', 1)
    xlim([-0.5, 2.5])    
    ylim([-2.5e-5, 4e-5])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    box off    
end
subplot(2,4,5)
set(gca, 'xtick', [0 2])
set(gca, 'ytick', [-2e-5 0 2e-5 4e-5])
set(gca, 'yticklabel', {'-20' '0' '20' '40'})
xlabel('Time (s)')
ylabel('Voltage (uV)')
%nice_figure(gcf, [fig_folder '3_dayexamples'],[8 3.5])

%% schematic of interpretation - figure 3
figure
days = [6, 15, 24];
well = 12;
chan = 22;
labels = {'Early', 'Middle', 'Late'};
for i=1:3
    subplot(3,1,i)
    plot(LFP_ker_smo{days(i)}{well}(:,:,chan), 'color', [colors(1,:), 0.5], 'linewidth',1)
    xlim([0 2500])
    axis off
    title(labels(i))
end
nice_figure(gcf, [fig_folder '3_interp'],[1.2 5])

%% timeseries & spectral examples - figure 4
ctc_osc = load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/organoids_processed/CTC_110416/LFP_Sp.mat');
ctc_cpx = load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/organoids_processed/CTC_022417/LFP_Sp.mat');
%fetal = load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/fetal_122016/LFP_Sp.mat');
%ecog = load('/Users/rgao/Documents/Data/Kahana/MerkKaha15-selected/data_unpacked/UP043_19May14_1032.mat');
eeg = load('/Users/rdgao/Documents/data/AEDInfantEEG/eeg12466_17089.mat');
adult_eeg = load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/eeg_data.mat');
load('aggregate.mat','PSDw');
%%
eeg.fs = 167;
figure
y_pos = [0.79 0.55 0.3 0.075];
%data_all = {ctc_osc.LFP{12}(:,44), ctc_cpx.LFP{5}(:,15), fetal.LFP{7}(:,2), detrend(eeg.data(1790*eeg.fs:1850*eeg.fs,1),'constant');};
data_all = {ctc_osc.LFP{12}(:,44), ctc_cpx.LFP{5}(:,15),  detrend(eeg.data(1790*eeg.fs:1850*eeg.fs,1),'constant'), double(adult_eeg.data(4,740*256:800*256))'};
t_all = {ctc_osc.t_ds, ctc_cpx.t_ds, (0:length(data_all{3})-1)/eeg.fs, (0:length(data_all{4})-1)/adult_eeg.fs};

for i=1:4
   %data_all{i} = butterpass(data_all{i}, 1./mean(diff(t_all{i})), [0.1 55], 3);
   data_all{i} = eegfilt(data_all{i}', 1./mean(diff(t_all{i})), 0, 55)';
end
zoom_win = {[27 32], [40 45], [33.5 38.5], [13 18]};
labels = {'Organoid', 'Organoid', 'Preterm Neonatal EEG', 'Adult EEG'};
for i=1:4
    subplot('position',[0.08,y_pos(i),0.67, 0.16])   
    plot(t_all{i}, data_all{i}, 'linewidth',1)
    %set(gca,'yticklabel',[])
    xlim([0 20])
    set(gca,'xtick',[])
    box off
    title(labels{i})
    
    %put in label
    switch i
        case 2            
            ylim([-4e-5, 4e-5])
            yticks([-4e-5 0 4e-5])
            set(gca,'yticklabel',{'-40' '0' '40'})
        case 3
            ylim([-200, 200])
            yticks([-200 0 200])
            set(gca,'yticklabel',{'-200' '0' '200'})
        case 4
            ylim([-50 50])
            set(gca,'ytick', -50:50:50)
            set(gca,'yticklabel',{'-50' '0' '50'})
            set(gca,'xtick', 0:20:60)
            xlabel('Time (s)')
            ylabel('Voltage (uV)')
    end
    if i~= 4
        %draw box
        hold on
        xb = zoom_win{i};
        yb = ylim;
        fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
        set(fbox, 'facealpha', 0.1)
        hold off
    end
    set(gca, 'fontsize', 18)
    
    
%     if i==1
%         set(gca,'yticklabel',{'-40' '-20' '0' '20'})
%     end
%     if i==4
%         set(gca,'ytick', -200:200:200)
%         set(gca,'yticklabel',{'-200' '0' '200'})
%         set(gca,'xtick', 0:20:60)
%         xlabel('Time (s)')
%         ylabel('Voltage (uV)')
%     end
    
    %draw 
    subplot('position',[0.8,y_pos(i),0.16, 0.16])
    plot(t_all{i}, data_all{i}, 'linewidth',1)
    set(gca,'yticklabel',[])
    xlim(zoom_win{i})
    set(gca,'xtick',[])
    box off
    if i==4
        set(gca,'xtick', zoom_win{i})
        xlabel('Time (s)')
        set(gca,'xtickLabel', [0 5])
    end
end
%nice_figure(gcf, [fig_folder '4_timeseries'],[7 7])

%% spectrogram & PSD
figure
FS_all = {1000, 1000, 1000, eeg.fs};
PSD_all = {PSDw{18}{12}, ...
    PSDw{38}{5}, ...
    PSDw{38}{5}, ...
    %pwelch(fetal.LFP{7},fetal.fs_ds*2,fetal.fs_ds,fetal.fs_ds*2,fetal.fs_ds),...    
    pwelch(data_all{4},eeg.fs*2,eeg.fs,eeg.fs*2,eeg.fs)};
%%
chans = [44, 15, 2, 1];
freq_all = {0:0.5:45, 0:0.5:45, 0:0.5:45, 0:0.5:45};
osc_range = {[2.5 4.5], [3 5], [6 10], [6 10]};
%%
for i=1:4
    subplot('position',[0.07,y_pos(i),0.62, 0.16])   
    [S, F, T] = spectrogram(data_all{i}, FS_all{i}*2, round(FS_all{i}*1.75), FS_all{i}*2, FS_all{i});
    imagesc(T,F(1:50),(abs(S(1:50,:)))./max(max((abs(S(1:50,:))))))
    xlim([T(1) 60])
    ylim([0 10])
    set(gca,'xtick',[])
    box off
    title(labels{i})
    if i==1
        h = colorbar('Location','east');
        set(h,'ylim',[0 1])
        set(h,'ytick',[0 0.5 1])
        set(h,'ycolor','w')     
    end
    
    hold on    
    xb = zoom_win{i};
    yb = ylim+[0.2 -0.2];
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], 'w');    
    set(fbox, 'facecolor', 'none')
    set(fbox, 'edgecolor', 'w')
    hold off
    set(gca,'ytick',0:10:20)    
    
    if i==1
        ylabel('Frequency (Hz)')
    end
    if i==4
        set(gca,'xtick', 0:20:60)        
        xlabel('Time (s)')        
    end
    
    subplot('position',[0.8,y_pos(i),0.16, 0.16])
    loglog(freq_all{i}, PSD_all{i}(1:length(freq_all{i}),chans(i))./(max(PSD_all{i}(1:length(freq_all{i}),chans(i)))), 'linewidth',1)
    xlim([0 40])
    %ylim([min(PSD_all{i}(1:length(freq_all{i}),chans(i))), max(PSD_all{i}(1:length(freq_all{i}),chans(i)))])
    set(gca,'ytick',[0.01 1])
    set(gca,'yticklabel', {'0.01' '1'})
    set(gca,'xticklabel',[])
    box off
    if i==1
        ylabel('Power (uV^2/Hz)')
    end
    if i==4
        set(gca,'xtick', [1 10])
        set(gca,'xticklabel', {'1' '10'})        
        xlabel('Frequency (Hz)')
    end
end
%nice_figure(gcf, [fig_folder '4_spgpsd'],[6 6])

%% pharmacology figures
figure
colors = get(gca,'colororder');
% rec 25 & 26 are the CT/AP5+CNQX/Bicuculine pairs
% control spiking kernel
ctr_data = {kernels{25}{6},kernels{25}{12},LFP_ker_smo{25}{6}(:,:,28)*1e5,LFP_ker_smo{25}{12}(:,:,30)*1e5};
drug_data = {kernels{26}{6},kernels{26}{12},LFP_ker_smo{26}{6}(:,:,28)*1e5,LFP_ker_smo{26}{12}(:,:,30)*1e5};
ctr_colors = {'k' 'k' colors(1,:) colors(1,:)};
for i = 1:4
    subplot(2,2,i)
    plot_filled((-500:2500)/1000, mean(ctr_data{i},2), std(ctr_data{i},1,2), ctr_colors{i})
    hold on
    plot_filled((-500:2500)/1000, mean(drug_data{i},2), std(drug_data{i},1,2), colors(2,:))    
    hold off
    xlim([-500 2500]/1000)
    set(gca, 'xtick', '')
end
for i=1:2
    subplot(2,2,i)
    ylim([0 30])
    set(gca, 'ytick', [0 30])
end
for i=3:4
    subplot(2,2,i)
    set(gca, 'ytick', ylim)
end
subplot(2,2,1)
title('Control')
ylabel({'Population Spiking';'(Spike Count)'})
axes_h = gca;
legend(axes_h.Children([3 1]), {'Pre' 'Post'}, 'location', 'best')
subplot(2,2,2)
title('Bicuculline')
subplot(2,2,3)
set(gca, 'xtick', [0 1 2])
xlabel('Time (s)')
ylabel({'LFP'; '(Voltage)'})
axes_h = gca;
legend(axes_h.Children([3 1]), {'Pre' 'Post'}, 'location', 'best')
subplot(2,2,4)
set(gca, 'xtick', [0 1 2])
% plot significance
for i=1:4
    for tt=1:3000
        [hh(tt) pv(tt)] = ttest2(ctr_data{i}(tt,:),drug_data{i}(tt,:));
    end
    subplot(2,2,i)
    YL = ylim;
    hold on
    plot((find(pv<0.01)-500)/fs, ones(1,length(find(pv<0.01)))*YL(2), 'k.')
    hold off
end
subplot(2,2,1); ylim([0 31])
subplot(2,2,2); ylim([0 31])
subplot(2,2,3); ylim([-1.5 1.1])
subplot(2,2,4); ylim([-3 4.2])

%% pharmacology figures
figure
colors = get(gca,'colororder');
% rec 25 & 26 are the CT/AP5+CNQX/Bicuculine pairs
% control spiking kernel
ctr_data = {kernels{25}{6},kernels{25}{12}, kernels{34}{12}, LFP_ker_smo{25}{6}(:,:,28)*1e6,LFP_ker_smo{25}{12}(:,:,30)*1e6, LFP_ker_smo{34}{12}(:,:,30)*1e6};
drug_data = {kernels{26}{6},kernels{26}{12},kernels{35}{12}, LFP_ker_smo{26}{6}(:,:,28)*1e6,LFP_ker_smo{26}{12}(:,:,30)*1e6, LFP_ker_smo{35}{12}(:,:,30)*1e6};
ctr_colors = {'k' 'k' 'k' colors(1,:) colors(1,:) colors(1,:)};
for i = 1:6
    subplot(2,3,i)
    plot_filled((-500:2500)/1000, mean(ctr_data{i},2), std(ctr_data{i},1,2), ctr_colors{i})
    hold on
    plot_filled((-500:2500)/1000, mean(drug_data{i},2), std(drug_data{i},1,2), colors(2,:))    
    hold off
    xlim([-500 2500]/1000)
    set(gca, 'xtick', '')
end
for i=1:3
    subplot(2,3,i)
    ylim([0 30])
    set(gca, 'ytick', [0 30])
end
for i=4:6
    subplot(2,3,i)
    YL = ylim;
    set(gca, 'ytick', [YL(1) 0 YL(2)])
end

subplot(2,3,1)
title('Control (6 Month)')
ylabel({'Population Spiking';'(Spike Count)'})
axes_h = gca;
legend(axes_h.Children([3 1]), {'Pre' 'Post'}, 'location', 'best')
subplot(2,3,2)
title('Bicu. (6 Month)')
subplot(2,3,3)
title('Bicu. (8 Month)')

subplot(2,3,4)
set(gca, 'xtick', [0 1 2])
xlabel('Time (s)')
ylabel({'LFP'; '(uV)'})
axes_h = gca;
legend(axes_h.Children([3 1]), {'Pre' 'Post'}, 'location', 'best')
subplot(2,3,5)
set(gca, 'xtick', [0 1 2])
subplot(2,3,6)
set(gca, 'xtick', [0 1 2])
%%
% plot significance
for i=1:6
    for tt=1:3000
        [hh(tt) pv(tt)] = ttest2(ctr_data{i}(tt,:),drug_data{i}(tt,:));
    end
    subplot(2,3,i)
    YL = ylim;
    hold on
    plot((find(pv<0.01)-500)/fs, ones(1,length(find(pv<0.01)))*YL(2), 'k.')
    hold off
end
%%
subplot(2,2,1); ylim([0 31])
subplot(2,2,2); ylim([0 31])
subplot(2,2,3); ylim([-1.5 1.1])
subplot(2,2,4); ylim([-3 4.2])

%% for whatever reason this plot takes forever to make, probably the shading
nice_figure(gcf, [fig_folder '4_pharm'],[6 4])
%% first pharm exp (25/26)
num_events = [cellfun(@length,kerTimes{25}); cellfun(@length,kerTimes{26})];
figure
subplot(1,2,1)
h1 = plot(num_events(:,[5 6 9 10]), 'o-', 'color', 'k', 'markerfacecolor', 'k');
hold on
h2 = plot(num_events(:,[7 11]), 'o-', 'color', colors(3,:), 'markerfacecolor', colors(3,:));
h3 = plot(num_events(:,[8 12]), 'o-', 'color', colors(2,:), 'markerfacecolor', colors(2,:));
hold off
xlim([0.8 2.2])
%legend([h1(1) h2(1) h3(1)], {'Ctrl', 'CNQX+AP5', 'Bicuculine'}, 'Location', 'northwest')
set(gca, 'ytick', [0 10 20])
ylabel('Network Events')
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel', {'Pre', 'Post'})
title('6 Months')
%nice_figure(gcf, [fig_folder '4_pharm_events1'],[2.25 4])

% second pharm exp (34/35)
num_events = [cellfun(@length,kerTimes{34}); cellfun(@length,kerTimes{35})];
%figure
subplot(1,2,2)
h1 = plot(num_events(:,[5 9]), 'o-', 'color', 'k', 'markerfacecolor', 'k');
hold on
h2 = plot(num_events(:,[6 10]), 'o-', 'color', colors(1,:), 'markerfacecolor', colors(1,:));
h3 = plot(num_events(:,[7 11]), 'o-', 'color', colors(3,:), 'markerfacecolor', colors(3,:));
h4 = plot(num_events(:,[8 12]), 'o-', 'color', colors(2,:), 'markerfacecolor', colors(2,:));
hold off
xlim([0.8 2.2])
legend([h1(1) h2(1) h3(1) h4(1)], {'Ctrl', 'Baclofen', 'CNQX+AP5', 'Bicuculline'}, 'Location', 'best', 'fontsize', 10)
set(gca, 'ytick', [0 30 60])
%ylabel('Network Events')
set(gca, 'xtick', [1 2])
set(gca, 'xticklabel', {'Pre', 'Post'})
title('8 Months')
%nice_figure(gcf, [fig_folder '4_pharm_events2'],[2.25 4])
nice_figure(gcf, [fig_folder '4_pharm_events'],[5 3.5])

%% PAC & MI plots
load('aggregate.mat', 'MI_all')
figure
plot([1 1.5], squeeze(mean(MI_all,1)), '-ok', 'markerfacecolor', 'k', 'linewidth', 1)
hold on
plot(1.5, squeeze(mean(MI_all(:,2,:),1)), 'o','color',colors(2,:), 'markerfacecolor', colors(2,:))
hold off
xlim([0.8 1.7])
set(gca,'xtick',[1 1.5])
set(gca,'xticklabel',{'Event' 'Non-Event'})
ylabel('Modulation Index')
%nice_figure(gcf, [fig_folder '4_MI'],[2 3])

%% plot subset of spike waveforms
load('organoids_processed/CTC_101816/LFP_Sp.mat', 'spike_shape')
figure
well = 10;
plot_win = (-10:10)/12500*1000;
for chan=1:9
    subplot(3,3,chan)
    plot(plot_win, spike_shape{well,chan}(41:61,:), 'color', [0 0 0 0.1], 'linewidth', 1)
    hold on
    idx = spike_shape{well,chan}(51,:)>0;
    plot(plot_win, mean(spike_shape{well,chan}(41:61,idx==1),2), 'color', colors(1,:), 'linewidth', 2);
    plot(plot_win, mean(spike_shape{well,chan}(41:61,idx==0),2), 'color', colors(2,:), 'linewidth', 2)
    hold off
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    xlim(plot_win([1 end]))
end
subplot(3,3,7)
set(gca,'xtick',[plot_win(1) 0 plot_win(end)])
xlabel('Time (ms)')
ylabel('Voltage')
nice_figure(gcf, [fig_folder 'SI_spkwaveform'],[8,8])

%% OLD CODE %%%%%%%%%%%%%%%%%%%%%%%%
%
%
% OLD CODE
%
%
%
% first recording, 7/31 is week 8
load names.mat
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
dayVec = (1:length(recs))+7; %use weeks instead
fs = 12500;
fs_ds = 1000;
dF = 0.5;
Fa = [0:dF:59 61:dF:200]; %frequency points
Xa = Fa/dF+1;
%% load data for plotting
load aggregate.mat
load('organoids_processed/CTC_101816/LFP_Sp.mat', 'LFP', 'spikes', 'spike_shape', 't_s', 't_ds')
%% schematic of data
% plot example raster
% recording Oct 18, 2016, corresponds to index 15 in data matrix
pkts = [6.319 43.02];
date = find(strcmp({dates.name}, 'CTC_101816'));
well = 12;
chan = 53;
XL = [0 60];
figure
subplot(2,1,1)
spike_raster(t_s,spikes,well)
xlabel('')
set(gca, 'xtick', [])
set(gca, 'ytick',[0,60])
xlim(XL)
title('Spike Raster')
hold on
for pk = 1:2
    xb = pkts(pk) +[-0.5 2.5];
    yb = ylim;
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
    set(fbox, 'facealpha', 0.1)
end
hold off

% plot network spike
subplot(4,1,3)
plot(t_ds, nws_smo{date}(:,well), 'k', 'linewidth', 1)
set(gca, 'xtick', [])
ylabel('Spikes')
xlim(XL)
ylim([0 20])
box off
title('Population Spiking')
hold on
for pk = 1:2
    xb = pkts(pk) +[-0.5 2.5];
    yb = ylim;
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
    set(fbox, 'facealpha', 0.1)
end
hold off

% plot example LFP
subplot(4,1,4)
plot(t_ds, LFP{well}(:,chan), 'linewidth', 1)
xlim(XL)
set(gca, 'ytick', 0)
set(gca, 'yticklabel', '0')
xlabel('Time (s)')
ylabel('Voltage')
box off
title('Local Field Potential (LFP)')
hold on
for pk = 1:2
    xb = pkts(pk) +[-0.5 2.5];
    yb = ylim;
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
    set(fbox, 'facealpha', 0.1)
end
hold off
nice_figure(gcf, [fig_folder '3a'],[6 6])

%% big time trace of pop spiking and LFP
% plot example of spiking kernel
well = 12;
chan = 53;
date = find(strcmp({dates.name}, 'CTC_101816'));
t = (-500:2500)/1000;

figure
subplot(2,2,1)
plot(t,kernels{date}{well},'color', [0 0 0 0.5], 'linewidth', 1)
ylabel('Spikes')
xlabel('Time (s)')
set(gca, 'xtick', [0 2])
set(gca, 'ytick', [0 20])
xlim([-0.5, 2.5])
box off
title('Pop. Spiking Kernel')

% plot example of LFP kernel
subplot(2,2,3)
plot(t,LFP_ker_smo{date}{well}(:,:,chan),'color', [colors(1,:) 0.5], 'linewidth', 1)
xlabel('Time (s)')
ylabel('Voltage')
set(gca, 'xtick', [0 2])
set(gca, 'ytick', 0)
set(gca, 'yticklabel', '0')
xlim([-0.5, 2.5])
box off
title('LFP Kernel')

% sample time trace from different wells on same day
% plot kernels from different wells
wells = [5 6 7 8];
subp_idx = [3 4 7 8];
for well = 1:4   
    subplot(4,4,subp_idx(well))
    plot(t,kernels{date}{wells(well)},'color', [0 0 0 0.5], 'linewidth', 1)
    box off
    xlim([-0.5, 2.5])
    ylim([0 25])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
end

% plot LFP kernels from different wells
chans = [9 7 3 35];
subp_idx = [11 12 15 16];
for well = 1:4 
    subplot(4,4,subp_idx(well))
    plot(t,LFP_ker_smo{date}{wells(well)}(:,:,chans(well)),'color', [colors(1,:) 0.5], 'linewidth', 1)
    box off
    xlim([-0.5, 2.5])    
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
end
nice_figure(gcf, [fig_folder '3b'],[6 6])

%% sample time trace from different days
% spiking kernel from different days
t = (-500:2500)/1000;
days = [2 8 15 24];
well = 5;
chan = 13;
figure
for day=1:4
    subplot(2,4,day)
    plot(t,kernels{recs(days(day))}{well},'color', [0 0 0 0.3], 'linewidth', 1)
    xlim([-0.5, 2.5])
    ylim([0 35])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    title([num2str(day*2), ' Months'])
    box off
end
subplot(2,4,1)
set(gca, 'ytick', [0 35])
ylabel('Pop. Spike')

% LFP kernel from different days
for day=1:4
    subplot(2,4,day+4)
    plot(t,LFP_ker_smo{recs(days(day))}{well}(:,:,chan),'color', [colors(1,:) 0.3], 'linewidth', 1)
    xlim([-0.5, 2.5])    
    ylim([-2.5e-5, 4e-5])
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    box off    
end
subplot(2,4,5)
set(gca, 'xtick', [0 2])
set(gca, 'ytick', 0)
set(gca, 'yticklabel', '0')
xlabel('Time (s)')
ylabel('Voltage')
nice_figure(gcf, [fig_folder '3c'],[8 3.5])

%% aggregate stats plot
%network inter-event latency, surrogate of event count
X = dayVec;
Y = nanmean(ker_latM(recs,:),2);
Y_s = nanstd(ker_latM(recs,:),1,2)/sqrt(8);
figure
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
set(gca,'ytick',0:50:150)
ylabel('Inter-Event Interval(s)')
nice_figure(gcf, [fig_folder '3d'],[3 3])

%total spiking under the event
Y = nanmean(ker_totalM(recs,:),2);
Y_s = nanstd(ker_totalM(recs,:),1,2)/sqrt(8);
figure
plot_filled(dayVec,Y,Y_s, 'k')
xlabel('Weeks')
set(gca,'ytick', 0:10000:30000)
set(gca,'yticklabel', num2str((0:3)'))
ylabel('Spikes In Event (x10^5)')
nice_figure(gcf, [fig_folder '3e'],[3 3])

%P1 amplitude
figure
Y = nanmean(ker_peakM(recs,:),2);
Y_s = nanstd(ker_peakM(recs,:),1,2)/sqrt(8);
subplot(2,1,1)
plot_filled(X,Y,Y_s,'k')
set(gca,'xticklabel',[])
ylabel('Peak 1 Amp.')

%P1 fullwidth half maximum
subplot(2,1,2)
Y = nanmean(cellfun(@mean,FWHM(recs,:)),2);
Y_s = nanstd(cellfun(@mean,FWHM(recs,:)),1,2)/sqrt(8);
plot_filled(X,Y,Y_s, 'k')
xlabel('Weeks')
ylabel('Peak 1 Width')
nice_figure(gcf, [fig_folder '3f'],[3 3])

%% aggregate subpeak plots
% plot number of subpeaks over time
Y = nanmean(nsubpks(recs,:),2);
Y_s = nanstd(nsubpks(recs,:)')'/sqrt(8);
figure
plot_filled(X, Y, Y_s, 'k')
set(gca,'ytick',0:1:6)
xlabel('Weeks')
ylabel('# of Subpeaks')
nice_figure(gcf, [fig_folder '3g'],[3 3])

%subpeak amplitude
Y = squeeze(nanmean(p2P(recs,:,2:4),2));
figure
subplot(2,1,1)
plot(X,Y)
set(gca,'xticklabel',[])
xlim([dayVec([1 end])])
ylim([0.2 0.8])
ylabel('Rel. Subpeak Amp.')
legend({'P2', 'P3', 'P4'}, 'location', 'northwest')

%subpeak delay
subplot(2,1,2)
Y = squeeze(nanmean(diff(p2T(recs,:,1:4),1,3),2))';
plot(X,Y)
xlim([dayVec([1 end])])
ylim([200 800])
xlabel('Weeks')
ylabel('Subpeak Latency (ms)')
legend({'P1-P2', 'P2-P3', 'P3-P4'}, 'location', 'northwest')
nice_figure(gcf, [fig_folder '3h'],[4 4])

%% event kernel stats
%event correlation
Y = nanmean(self_cor(recs,:),2);
Y_s = nanstd(self_cor(recs,:),1,2)/sqrt(8);
figure
plot_filled(dayVec,Y,Y_s,'k')
xlabel('Weeks')
ylabel('Pop. Spike Event Corr.')
nice_figure(gcf, [fig_folder '3i'],[3 3])

%LFP mean correlation
figure
subplot(2,1,1)
flat_corr = reshape(self_cor_LFP,size(self_cor_LFP,1),[]);
Y = nanmean(flat_corr(recs,:),2);
Y_s = nanstd(flat_corr(recs,:),1,2)/sqrt(8);
plot_filled(X, Y, Y_s, colors(1,:));
ylabel('LFP Mean Corr.')
set(gca,'xtick', [])

%LFP max correlation
subplot(2,1,2)
Y = nanmean(squeeze(max(self_cor_LFP(recs,:,:),[],3)),2);
Y_s = nanstd(squeeze(max(self_cor_LFP(recs,:,:),[],3)),1,2)/sqrt(8);
plot_filled(X, Y, Y_s, colors(1,:));
ylim([0.6 1])
ylabel('LFP Max Corr.')
xlabel('Weeks')
nice_figure(gcf, [fig_folder '3j'],[3 3])

% LFP spatial correlation
figure
Y = nanmean(LFP_spatialcorr(recs,:),2);
Y_s = nanstd(LFP_spatialcorr(recs,:),1,2)/sqrt(8);
plot_filled(X, Y, Y_s, colors(1,:))
xlabel('Weeks')
ylabel('LFP Spatial Correlation')
nice_figure(gcf, [fig_folder '3k'],[3 3])

%% schematic of interpretation
figure
days = [6, 15, 24];
well = 12;
chan = 22;
labels = {'Early', 'Middle', 'Late'};
for i=1:3
    subplot(3,1,i)
    plot(LFP_ker_smo{days(i)}{well}(:,:,chan), 'color', [colors(1,:), 0.5], 'linewidth',1)
    xlim([0 2500])
    axis off
    title(labels(i))
end
nice_figure(gcf, [fig_folder '3_interp'],[1.2 5])

%% figure 4 - the good shit --------------------------------
ctc = load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/organoids_processed/CTC_022417/LFP_Sp.mat');
fetal = load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/fetal_122016/LFP_Sp.mat');
eeg = load('/Users/rgao/Documents/Data/Kahana/MerkKaha15-selected/data_unpacked/UP043_19May14_1032.mat');
eeg = load('/Users/rgao/Documents/Data/Muotri/Pri_Corticoids/EEG/EEG.mat');

%% time trace
figure
y_pos = [0.79 0.55 0.3 0.075];
data_all = {ctc.LFP{5}(:,15), fetal.LFP{7}(:,2), sin(0:0.001:60*2*pi), eeg.data(eeg.fs*100:end,4)};
t_all = {ctc.t_ds, fetal.t_ds, (0:length(data_all{3})-1)/eeg.fs, (0:length(data_all{4})-1)/eeg.fs};
for i=1:4
   data_all{i} = butterpass(data_all{i}, 1./mean(diff(t_all{i})), [0.01 60]);
end
zoom_win = {[40 45], [13 18], [40 45], [38 43]};
labels = {'Organoid', 'Fetal Brain Tissue', 'Preterm Infant EEG', 'Adult EEG'};
for i=1:4
    subplot('position',[0.05,y_pos(i),0.7, 0.16])   
    plot(t_all{i}, data_all{i}, 'linewidth',1)
    set(gca,'yticklabel',[])
    xlim([0 60])
    set(gca,'xtick',[])
    box off
    title(labels{i})

    hold on    
    xb = zoom_win{i};
    yb = ylim;
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
    set(fbox, 'facealpha', 0.1)
    hold off
    
    if i==4
        set(gca,'xtick', 0:20:60)
        xlabel('Time (s)')
        ylabel('Voltage')
    end
    
    subplot('position',[0.8,y_pos(i),0.16, 0.16])
    plot(t_all{i}, data_all{i}, 'linewidth',1)
    set(gca,'yticklabel',[])
    xlim(zoom_win{i})
    set(gca,'xtick',[])
    box off
    if i==4
        set(gca,'xtick', zoom_win{i})
        xlabel('Time (s)')
        set(gca,'xtickLabel', [0 5])
    end
end
%nice_figure(gcf, [fig_folder '4_timeseries'],[6 6])

%% spectrogram & PSD
figure
y_pos = [0.79 0.55 0.3 0.075];
data_all = {ctc.LFP{5}(:,15), fetal.LFP{7}(:,2), sin(0:0.001:60*2*pi), eeg.data(eeg.fs*100:end,4)};
t_all = {ctc.t_ds, fetal.t_ds, (0:length(data_all{3})-1)/eeg.fs, (0:length(data_all{4})-1)/eeg.fs};
for i=1:4
   data_all{i} = butterpass(data_all{i}, 1./mean(diff(t_all{i})), [0.01 60]);
end
zoom_win = {[40 45], [13 18], [40 45], [38 43]};
labels = {'Organoid', 'Fetal Brain Tissue', 'Preterm Infant EEG', 'Adult EEG'};
FS_all = {1000, 1000, 1000, eeg.fs};

PSD_all = {PSDw{18}{12}, ...
    pwelch(fetal.LFP{8},fetal.fs_ds*2,fetal.fs_ds/2,fetal.fs_ds*2,fetal.fs_ds),...
    pwelch(eeg.data,eeg.fs*2,eeg.fs/2,eeg.fs*2,eeg.fs),...
    pwelch(eeg.data,eeg.fs*2,eeg.fs/2,eeg.fs*2,eeg.fs)};

freq_all = {0:0.5:55, 0:0.5:55, 0:0.5:55, 0:0.5:55};
osc_range = {[2.5 4.5], [3 5], [6 10], [6 10]};

for i=1:4
    subplot('position',[0.07,y_pos(i),0.65, 0.16])   
    [S, F, T] = spectrogram(data_all{i}, FS_all{i}*2, FS_all{i}*1.875, FS_all{i}*2, FS_all{i});
    imagesc(T,F,abs(S))    
    xlim([T(1) 60])
    ylim([0 10])
    set(gca,'xtick',[])
    box off
    title(labels{i})

    hold on    
    xb = zoom_win{i};
    yb = ylim;
    fbox = fill([xb(1) xb(2) xb(2) xb(1)], [yb(1) yb(1) yb(2) yb(2)], 'w');
    set(fbox, 'facealpha', 0.1)
    set(fbox, 'edgecolor', 'w')
    hold off
    set(gca,'ytick',0:10:10)
    
    if i==4
        set(gca,'xtick', 0:20:60)        
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
    end
    
    subplot('position',[0.8,y_pos(i),0.16, 0.16])
    P = log10(PSD_all{i}(1:length(freq_all{i}),:));
    plot_filled(log10(freq_all{i}+1e-10),mean(P,2), std(P,1,2), colors(1,:))   
    hold on
    xb = osc_range{i};
    yb = ylim;
%     fbox = fill(log10([xb(1) xb(2) xb(2) xb(1)]), [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
%     set(fbox, 'facealpha', 0.1)
    xlim([0 log10(55)])
    set(gca,'ytick',[])
    set(gca,'xtick',[])    
    box off
    hold off
    axis 'auto y'
    if i==4
        set(gca,'xtick', log10([1 10]))
        set(gca,'xticklabel', {'1' '10'})
        ylabel('Power (au)')
        xlabel('Frequency (Hz)')
    end
end
%nice_figure(gcf, [fig_folder '4_spgpsd'],[6 6])

%% PSD comparison

figure
PSD_all = {PSDw{18}{12}, ...
    pwelch(fetal.LFP{8},fetal.fs_ds*2,fetal.fs_ds/2,fetal.fs_ds*2,fetal.fs_ds),...
    pwelch(eeg.data,eeg.fs*2,eeg.fs/2,eeg.fs*2,eeg.fs),...
    pwelch(eeg.data,eeg.fs*2,eeg.fs/2,eeg.fs*2,eeg.fs)};

freq_all = {0:0.5:55, 0:0.5:55, 0:0.5:55, 0:0.5:55};
osc_range = {[2.5 4.5], [3 5], [6 10], [6 10]};
for i=1:4
   subplot(2,2,i)
   %loglog(freq_all{i}, PSD_all{i}, 'color',[0 0 0 0.05],'linewidth',1);
   %hold on
   %loglog(freq_all{i}, mean(PSD_all{i},2), 'color',colors(1,:),'linewidth',2);
   %hold off
   %xlim([0 55])
   P = log10(PSD_all{i}(1:length(freq_all{i}),:));
   plot_filled(log10(freq_all{i}+1e-10),mean(P,2), std(P,1,2), colors(1,:))   
   hold on
   xb = osc_range{i};
   yb = ylim;
   fbox = fill(log10([xb(1) xb(2) xb(2) xb(1)]), [yb(1) yb(1) yb(2) yb(2)], colors(3,:));
   set(fbox, 'facealpha', 0.1)
   xlim([0 log10(55)])
   set(gca,'ytick',[])
   set(gca,'xtick',[])
   title(labels{i})
   box off
   hold off
end
subplot(2,2,3)
set(gca,'xtick',[0.1 1])
set(gca,'xticklabel',[1 10])
xlabel('Frequency (Hz)')
ylabel('Power')
nice_figure(gcf, [fig_folder '4b'],[6 6])

%% osc power over time
% relative oscillatory power via fitting 1/f
orgPSD_idx = [11, 22, 36];
figure

% pool channels: do we want to do this here?
% P_flat = reshape(osc_power_foof,size(osc_power_foof,1),[]);
% plot_filled(dayVec, mean(P_flat(recs,:),2), std(P_flat(recs,:),1,2)/sqrt(size(P_flat,2)), colors(1,:))

plot_filled(dayVec, mean(mean(osc_power_foof(recs,:,:),3),2), std(mean(osc_power_foof(recs,:,:),3),1,2), colors(1,:))
hold on
for i=1:3
    plot(dayVec(find(recs==orgPSD_idx(i))), 0, '^', 'color', colors(i,:), 'markerfacecolor', colors(i,:), 'markersize', 8);
end
hold off
xlabel('Weeks')
ylabel('Oscillatory Power')
%nice_figure(gcf, [fig_folder '4c'],[6 3])

%%
% barplot of osc power
dF=0.5;
fit_freq_all = {1:0.5:10, 1:0.5:20, 0.5:0.5:20, 4:0.5:30};
osc_freq_all = {2.5:dF:4.5, 2.5:dF:4.5, 6:dF:10, 6:dF:9};
labels = {'Early', 'Mid.', 'Late', 'Fetal', 'Baby', 'Adult'};
figure
hold on
% first 3 bars are time points of organoids
for i=1:3
    %pool channels
    %OP = reshape(squeeze(osc_power_foof(orgPSD_idx(i),:,:)),1,[]);
    
    OP = mean(osc_power_foof(orgPSD_idx(i),:,:),3);
    bar(i, mean(OP), 'facecolor', colors(i,:));
    plot(i+(rand(size(OP))-0.5)*0.3,OP, '.k', 'markersize',8)
    errorbar(i,mean(OP), std(OP), 'k', 'linewidth', 2);
end

% last 3 bars are fetal, baby, and adult ECoG
for i=2:4
    foof_res = foof_power(PSD_all{i}, fit_freq_all{i}, dF, osc_freq_all{i});
    OP = foof_res(:,3);
    bar(i+2, mean(OP), 'facecolor', 'w');
    plot(i+2+(rand(size(OP))-0.5)*0.3,OP, '.k', 'markersize',5)
    errorbar(i+2,mean(OP), std(OP), 'k', 'linewidth', 2);
    
end
hold off
xlim([0 7])
set(gca,'xtick', 1:6)
set(gca,'xticklabel', labels)
ylabel('Oscillatory Power')
nice_figure(gcf, [fig_folder '4d'],[4 4])



%% -------------------------
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % 

%% probably supplementary ------------------------
% PSD over time
C = redblue(length(kernels));
figure
well=5;
for rec=recs
    loglog(Fa,median(PSDw{rec}{well}(Xa,:),2), 'color', [1-C(rec,:) 0.5], 'linewidth',1)
    hold on
end
hold off
xlim([1 200])
xlabel('Frequency (Hz)')
set(gca,'ytick',[])
title('Well median over time; B: earlier, R: later')
nice_figure(gcf, [fig_folder 'supp_PSDovertime'],[4 4])

%% LFP event correlation over time for each well
for well = 1:8
    figure
    for i = 1:length(recs)
        subplot(6,5,i)
        imagesc(reshape(squeeze(self_cor_LFP(recs(i),well,:)),8,8))
        set(gca,'xticklabel', [])
        set(gca,'yticklabel', [])
        caxis([0 1])
        title(dayVec(i))
    end
    nice_figure(gcf, [fig_folder 'lfpcorrs/well' num2str(well)],[8 8])
    close
end



%% 1/f slope over time at different freqs
figure
for slp=1:3
    subplot(1,3,slp)
    plot_filled(dayVec, mean(median(slopes(recs,:,:,slp),3),2), std(median(slopes(recs,:,:,slp),3),1,2), 'k')
    hold on
    plot(dayVec, mean(slopes(recs,:,:,slp),3), 'k.')
    hold off
    set(gca,'xtick',[])
    title(['[' num2str(slope_range{slp}([1 end])) ']'])
end
subplot(1,3,1)
set(gca,'xtick',10:5:35)
xlabel('Weeks')
ylabel('Slope')
nice_figure(gcf, [fig_folder 'supp_slopes'],[12 4])

%% plot all kernels
for well = 5:12
    figure
    for i=1:length(recs)
        subplot(5,6,i)
        if ~isempty(kernels{recs(i)}{well})
            norm = repmat(max(kernels{recs(i)}{well}),length(kernels{recs(i)}{well}),1);
            plot((-500:2500)/1000,kernels{recs(i)}{well}./norm, 'color', [0 0 0 0.2], 'linewidth', 1);
            line([-0.5 2.5], [0.01 0.01], 'linewidth',1)
            line([-0.5 2.5], [0.05 0.05], 'linewidth',1)
            set(gca, 'xtick', [])
            set(gca, 'ytick', [])
            xlim([-0.5 2.5])
            title(['Week ', num2str(dayVec(i))])
        end
    end
    subplot(5,6,25)
    set(gca, 'xtick', [0, 2])
    xlabel('Time (s)')
    nice_figure(gcf, [fig_folder 'popkernels/well' num2str(well)],[8 8])
    close
end

%% LFP kernel histogram plot
figure
A = reshape(self_cor_LFP,length(kernels),[]);
H = hist(A',0:0.05:1);
imagesc(dayVec, 0.05:0.05:0.95, H(2:end-1,recs))
set(gca, 'ydir', 'normal')
hold on
plot(dayVec, nanmean(A(recs,:),2), 'w', 'linewidth',3)
hold off
colorbar
%xlabel('Days')
xlabel('Weeks')
ylabel('Similarity')
nice_figure(gcf, [fig_folder '5_longitudinal_LFPKerCorrDist'],[4 3])

%%
% relative oscillatory power via baseline subtraction
figure
plot_filled(dayVec, mean(mean(osc_power(recs,:,:),3),2), std(mean(osc_power(recs,:,:),3),1,2), colors(1,:))
xlabel('Weeks')
ylabel('Oscillatory Power')
nice_figure(gcf, [fig_folder 'supp_oscPW_longitudinal'],[6 3])


%% consistency plot for Fig 2
data_2d3d = {FR2D(:,2:13), FR3D(1:13,2:end)};
CC = {};
PV = {};
CC_val = {};
figure
hold on
for set = 1:2
    CC{set} = zeros(size(data_2d3d{set},2));
    PV{set} = zeros(size(data_2d3d{set},2));
    for i=1:size(data_2d3d{set},2)
        for j=(i+1):size(data_2d3d{set},2)
            %plot(data_2d3d{set}(:,i),data_2d3d{set}(:,j),'o')
            X = data_2d3d{set}(:,i);
            Y = data_2d3d{set}(:,j);
            good_idx = find((isnan(X)+isnan(Y))==0);
            [CC{set}(i,j), PV{set}(i,j)] = corr(X(good_idx), Y(good_idx));
            
        end
    end    
    CC_val{set} = CC{set}(find(CC{set}));

    plot(set, CC_val{set}, 'ok')
end
xlim([0 3])
set(gca, 'xtick', [0, 2])
set(gca,'xticks',[1 2])
%set(gca,'xticklabel', {'2D' , '3D'})