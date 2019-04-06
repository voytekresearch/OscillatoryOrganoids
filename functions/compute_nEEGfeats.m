function EMAfeatures = compute_nEEGfeats(psd_freq, psd_cell, event_feats)

num_feat = 23;
num_wells = size(event_feats.ED_05,2);
num_recs = length(psd_cell);

% stack all features in 3D array
EMAfeatures = zeros(num_recs, num_wells, num_feat);

% 7: events per hour
EMAfeatures(:,:,7) = event_feats.EPH;

% 8-11: RMS, 50%, 5%, 95%, event (SAT) duration
EMAfeatures(:,:,8) = event_feats.ED_rms;
EMAfeatures(:,:,9) = event_feats.ED_50;
EMAfeatures(:,:,10) = event_feats.ED_05;
EMAfeatures(:,:,11) = event_feats.ED_95;

% 12-15: RMS, 50%, 5%, and 95% inter-event duration
EMAfeatures(:,:,12) = event_feats.IEI_rms;
EMAfeatures(:,:,13) = event_feats.IEI_50;
EMAfeatures(:,:,14) = event_feats.IEI_05;
EMAfeatures(:,:,15) = event_feats.IEI_95;

% loopdiloop for spectral features
% delta (0-3), theta(3-8), alpha(8-15), beta(15-30)
band_bounds = [0, 3, 8, 15, 30];
for r_idx = 1:num_recs
    for well = 1:num_wells
        % 20-23: relative spectral power        
        % sum to 1 from 0-30Hz
        psd = psd_cell{r_idx}{well}(find(psd_freq<=band_bounds(end)),:);
        psd_norm = psd./repmat(sum(psd,1),length(find(band_bounds(end))),1);
        for bb = 1:length(band_bounds)-1
            bb_inds = find(psd_freq>=band_bounds(bb) & psd_freq<band_bounds(bb+1));
            EMAfeatures(r_idx, well, 20+bb-1) = mean(sum(psd_norm(bb_inds,:),1));
        end
    end
end
end

