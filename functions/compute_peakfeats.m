function [event_feats, pktimes] =  compute_peakfeats(nws_smo, fs, fp_params)
num_recs = length(nws_smo);
num_wells = size(nws_smo{1},2);

disp('Peak finding...')
% find peaks and get peak features
pktimes = cell(1,num_recs);
for r_idx = 1:num_recs   
    % disp(dates(r_idx).name)
    for well = 1:num_wells        
        pktimes{r_idx}{well} = kernel_findpeak(nws_smo{r_idx}(:,well), fp_params.MPH, fp_params.MPH_ratio, fp_params.MPD, fp_params.toss_thresh_ratio);
    end
end

disp('Computing features...')
% compute features
EPH = zeros(num_recs, num_wells);
%recDays_IEI_flat = cell(1,num_recs);
%recDays_ED_flat = cell(1,num_recs);
IEI_flat = cell(1,num_recs);
ED_flat = cell(1,num_recs);

for r_idx = 1:num_recs    
    for well=1:num_wells
        n_events = size(pktimes{r_idx}{well},1);
        % ------- EVENTS PER HOUR ------ % 
        EPH(r_idx, well) = (n_events/(size(nws_smo{1},1)/fs))*3600;
        if n_events>0
            % ----------- IEI FEATURES -----------------
            % time vector for regression and stuff
            % recDays_IEI_flat{r_idx} = [recDays_IEI_flat{r_idx} recDays(r_idx).*ones(1,n_events-1)];
            
            % inter event interval
            iei = diff(pktimes{r_idx}{well}(:,2))/fs;
            IEI_flat{r_idx} = [IEI_flat{r_idx} iei'];
            IEI_mean(r_idx,well) = mean(iei);
            IEI_std(r_idx,well) = std(iei);
            
            % quantile features (for baby EEG)
            IEI_rms(r_idx,well) = mean(iei.^2).^0.5;
            IEI_50(r_idx, well) = quantile(iei, 0.5);
            IEI_05(r_idx, well) = quantile(iei, 0.05);
            IEI_95(r_idx, well) = quantile(iei, 0.95);
            
            % ----------- EVENT DURATION FEATURES -----------------
            % recDays_ED_flat{r_idx} = [recDays_ED_flat{r_idx} recDays(r_idx).*ones(1,n_events)];
            
            % event duration
            ev_dur = (pktimes{r_idx}{well}(:,3)-pktimes{r_idx}{well}(:,1))/fs;
            ED_flat{r_idx} = [ED_flat{r_idx} ev_dur'];
            ED_mean(r_idx,well) = mean(ev_dur);
            ED_std(r_idx,well) = std(ev_dur);
            
            ED_rms(r_idx,well) = mean(ev_dur.^2).^0.5;
            ED_50(r_idx, well) = quantile(ev_dur, 0.5);
            ED_05(r_idx, well) = quantile(ev_dur, 0.05);
            ED_95(r_idx, well) = quantile(ev_dur, 0.95);            
        end
    end
end

% look i'm not proud of this either but it will have to do for now
event_feats = struct();
event_feats.EPH = EPH;
event_feats.ED_rms = ED_rms;
event_feats.ED_50 = ED_50;
event_feats.ED_05 = ED_05;
event_feats.ED_95 = ED_95;
event_feats.IEI_rms = IEI_rms;
event_feats.IEI_50 = IEI_50;
event_feats.IEI_05 = IEI_05;
event_feats.IEI_95 = IEI_95;
