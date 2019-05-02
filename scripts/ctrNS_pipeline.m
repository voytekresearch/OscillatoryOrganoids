do_convert = 1;
wells2process = 1:12;
batch_folder = '~/Documents/data/Muotri/Neurosphere/';
cd(batch_folder)
if do_convert
    p1_convertpreprocess(batch_folder, wells2process)
end 
%% collect
D = dir(batch_folder);
D = D(3:end); % ignore ./ and ../
D = D([D.isdir]);
psd = {};
nws_smo = {};
recday = zeros(1,length(D));
for f = 1:length(D)
    % get date string
    date_str = D(f).name(end-10:end-5);
    recday(f)=datenum(date_str,'mmddyy');
    cur_data = load([D(f).name, '/summary.mat']);
    psd{f} = cur_data.PSDw(wells2process);
    nws_smo{f} = cur_data.nws_smo(wells2process,:)';
end
dayVec = (recday-min(recday))/7+10;

%%
% computing features
fs = 1000;
fp_params.MPH = 1;
fp_params.MPH_ratio = 0.8;
fp_params.MPD = 200;
fp_params.toss_thresh_ratio = 0.001;
psd_freq = 0:0.5:500;
[event_feats, pk_times] = compute_peakfeats(nws_smo, fs, fp_params);
ctrlFC_EMAfeatures = compute_nEEGfeats(psd_freq, psd, event_feats);
save([batch_folder, 'ctrlFC_EMAfeatures.mat'], 'dayVec', 'ctrlFC_EMAfeatures')