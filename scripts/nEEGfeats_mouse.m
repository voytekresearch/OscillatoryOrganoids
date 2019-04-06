%% script for computing mouse primary culture (mpc) nEEG features
warning('off')
datapath = '/Users/rdgao/Documents/data/Muotri/Spencer/';
fs = 1000;
fp_params.MPH = 1;
fp_params.MPH_ratio = 0.8;
fp_params.MPD = 200;
fp_params.toss_thresh_ratio = 0.001;
psd_freq = 0:0.5:500;

files = {'39127_summary.mat', '41111_summary.mat'};
wt_inds = {[1 2 4 5 6 9], [1 2 5 6 7 8 9 10]};

%% 
for f = 1:length(files)    
    disp('Loading...')
    load([datapath,files{f}])
    disp([datapath,files{f}])
    
    DIV = [agg_summary.date_num] - agg_summary(1).date_num+3;
    dayVec = DIV/7; % this is fucking named dayVec even though it's weeks. Why do I do this to myself.
    for i=1:length(agg_summary)
        nws_smo{i} = agg_summary(i).nws_smo(wt_inds{f},:)';
        psd{i} = agg_summary(i).PSDw(wt_inds{f});
    end
    
    % computing features
    [event_feats, pk_times] = compute_peakfeats(nws_smo, fs, fp_params);
    mpc_EMAfeatures = compute_nEEGfeats(psd_freq, psd, event_feats);
    
    save([datapath, 'mpc', int2str(f), '_EMAfeatures.mat'], 'dayVec', 'mpc_EMAfeatures')
    clear nws_smo psd
end


