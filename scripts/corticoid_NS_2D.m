%% neurosphere data
cd('/Users/rgao/Documents/data/Muotri/Pri_Corticoids/neurosphere/')
load LFP_Sp.mat
%% 
fs = 12500;
bsp = binarize_spikes(ceil(t_ds), fs, spikes,fs_ds);
nws = (squeeze(sum(bsp,2))');
nws_smo = tsmovavg(nws','s',50)';
plot_tight(nws_smo(:,1:end), [3 4], [], [0, 10])
%% SI figure - example of neurosphere data
well=5;
chan = 6;
XL = [0 60];
figure
subplot(2,1,1)
spike_raster(t_s,spikes,well)
xlabel('')
set(gca, 'xtick', [])
set(gca, 'ytick',[0,60])
xlim(XL)
title('Spike Raster')

% plot network spike
subplot(4,1,3)
plot(t_ds, nws_smo(1:length(t_ds),well), 'k', 'linewidth', 1)
set(gca, 'xtick', [])
ylabel('Spikes')
xlim(XL)
ylim([0 10])
box off
title('Population Spiking')

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
nice_figure(gcf, [fig_folder 'SI_NS_raster'],[6 6])
%% SI figure - example of neurosphere kernel
[ker, ker_ind] = kernel_findpeak(nws_smo(:,well), [-500, 1000], 0.75, fs_ds);
%make smoothed LFP
smo_mask = gausswin(100)./sum(gausswin(100));
LFP_smo = conv(LFP{well}(:,chan), smo_mask,'same');
LFP_smo = LFP{well}(:,chan);
% get LFP kernel
lfp_ker = collect_spikes(LFP_smo,[],ker_ind,[-500 1000]);

figure
subplot(2,1,1)
plot((-500:1000)/1000, ker, 'color', [0 0 0 0.2], 'linewidth', 1)
ylabel('Spikes')
subplot(2,1,2)
plot((-500:1000)/1000, lfp_ker, 'color', [colors(1,:) 0.2], 'linewidth', 1)
xlabel('Time (s)')
ylabel('Voltage')
nice_figure(gcf, [fig_folder 'SI_NS_kernel'],[3 6])

%% 2D 
cd('/Users/rgao/Documents/data/Muotri/Pri_Corticoids/2D')
load LFP_Sp.mat 
%% 
fs = 12500;
bsp = binarize_spikes(ceil(t_ds(end)), fs,spikes,fs_ds);
nws = (squeeze(sum(bsp,2))');
nws_smo = tsmovavg(nws','s',50)';
plot_tight(nws_smo(:,1:end), [3 4], [], [0, 10])
%%
figure
plot(t_ds, LFP{7}(:,6)*1e5)
hold on
plot((0:length(nws_smo)-1)/fs_ds,nws_smo(:,7), 'k')
hold off
spike_raster(t_s,spikes,7)