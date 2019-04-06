load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/aggregate.mat', 'PSDw', 'recs');
load('/Users/rdgao/Documents/data/Muotri/Pri_Corticoids/event_timing.mat');

num_feat = 23;
psd_freq = 0:0.5:500;

% stack all features in 3D array
ctc_EMAfeatures = zeros(length(recs), length(wells), num_feat);

% 7: events per hour
ctc_EMAfeatures(:,:,7) = EPH;

% 8-11: RMS, 50%, 5%, 95%, event (SAT) duration
ctc_EMAfeatures(:,:,8) = ED_rms;
ctc_EMAfeatures(:,:,9) = ED_50;
ctc_EMAfeatures(:,:,10) = ED_05;
ctc_EMAfeatures(:,:,11) = ED_95;

% 12-15: RMS, 50%, 5%, and 95% inter-event duration
ctc_EMAfeatures(:,:,12) = IEI_rms;
ctc_EMAfeatures(:,:,13) = IEI_50;
ctc_EMAfeatures(:,:,14) = IEI_05;
ctc_EMAfeatures(:,:,15) = IEI_95;

% loopdiloop for spectral features
for r_idx = 1:length(recs)
    rec = recs(r_idx);    
    for well = wells               
        % 20-23: relative spectral power
        % delta (0-3), theta(3-8), alpha(8-15), beta(15-30)
        % sum to 1 from 0-30Hz
        psd = PSDw{rec}{well}(1:61,:);
        psd_norm = psd./repmat(sum(psd,1),61,1);
        ctc_EMAfeatures(r_idx, well-4, 20) = mean(sum(psd_norm(1:7,:),1));
        ctc_EMAfeatures(r_idx, well-4, 21) = mean(sum(psd_norm(8:17,:),1));
        ctc_EMAfeatures(r_idx, well-4, 22) = mean(sum(psd_norm(18:31,:),1));
        ctc_EMAfeatures(r_idx, well-4, 23) = mean(sum(psd_norm(32:end,:),1));                
    end
end

