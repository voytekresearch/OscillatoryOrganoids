function MEA_summary(agg_file, output_folder, wells, print_report)
% MEA_summary(agg_file, output_folder, wells, print_report)
%   load a LFP_Sp.mat file and produce a print out of network spike vector
%   and PSDs, and save out summary file to disk

disp('Loading file.')
load(agg_file)
fs = mean(1./diff(t_s));

% binarize spikes and get network spike
disp('Computing.')
bsp = squeeze(binarize_spikes(ceil(t_s(end)), fs, spikes, 1000));
nws = squeeze(sum(bsp,2));
nws_smo = zeros(size(nws));

for well=wells
    PSDm{well} = mPSD(LFP{well}, fs_ds, fs_ds*2, fs_ds/2, 500);
    PSDw{well} = pwelch(LFP{well},fs_ds*2,fs_ds,fs_ds*2,fs_ds);
    nws_smo(well,:) = conv(nws(well,:),gausswin(50),'same');    
end
ac = autocorr(nws', 5000);

% plotting
disp('Plotting.')
figure(1)
figure(2)
for well=wells
    figure(1)
    subplot(3,4,well)
    loglog(0:0.5:fs_ds/2, PSDm{well})
    xlim([1,fs_ds/2])
    axis tight
    figure(2)
    subplot(3,4,well)
    loglog(0:0.5:fs_ds/2, PSDw{well})
    xlim([1,fs_ds/2])
    axis tight
end
figure(1)
title('Median PSD')
figure(2)
title('Welch PSD')

plot_tight(nws_smo',[3,4],[],[])
plot_tight(ac, [3,4],[],[])


% saving plots
disp('Saving out.')
if ~exist(output_folder,'file')
    mkdir(output_folder);
end
if print_report    
    figure(1)
    nice_figure(gcf, [output_folder 'medianpsd'], [12 9])
    figure(2)
    nice_figure(gcf, [output_folder 'welchpsd'], [12 9])
    figure(3)
    nice_figure(gcf, [output_folder 'networkspikes'], [24 9])
    figure(4)
    nice_figure(gcf, [output_folder 'nws_ac'], [12 9])
end
% close all
% save([output_folder, 'summary.mat'], 'PSDm', 'PSDw', 'nws', 'nws_smo')

