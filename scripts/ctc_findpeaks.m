%% find peaks and get peak features
for r_idx = 1:length(nws_smo)
    disp(dates(r_idx).name)
    for well = wells
        pktimes{r_idx}{well} = kernel_findpeak(nws_smo{r_idx}(:,well), MPH, MPH_ratio, MPD, toss_thresh_ratio);
    end
end

%%
recDays_IEI_flat = cell(1,length(recs));
recDays_ED_flat = cell(1,length(recs));
IEI_flat = cell(1,length(recs));
ED_flat = cell(1,length(recs));
for r_idx = 1:length(recs)
    rec = recs(r_idx);
    for well=wells
        n_events = size(pktimes{rec}{well},1);
        % ------- EVENTS PER HOUR ------ %         
        EPH(r_idx, well-4) = n_events/T{rec}*3600;
        
        
        if n_events>0
            % ----------- IEI FEATURES -----------------
            % time vector for regression and stuff
            recDays_IEI_flat{r_idx} = [recDays_IEI_flat{r_idx} recDays(r_idx).*ones(1,n_events-1)];
            
            % inter event interval
            iei = diff(pktimes{rec}{well}(:,2))/fs;
            IEI_flat{r_idx} = [IEI_flat{r_idx} iei'];
            IEI_mean(r_idx,well-4) = mean(iei);
            IEI_std(r_idx,well-4) = std(iei);
            
            % quantile features (for baby EEG)
            IEI_rms(r_idx,well-4) = mean(iei.^2).^0.5;
            IEI_50(r_idx, well-4) = quantile(iei, 0.5);
            IEI_05(r_idx, well-4) = quantile(iei, 0.05);
            IEI_95(r_idx, well-4) = quantile(iei, 0.95);
            
            % ----------- EVENT DURATION FEATURES -----------------
            recDays_ED_flat{r_idx} = [recDays_ED_flat{r_idx} recDays(r_idx).*ones(1,n_events)];
            
            % event duration
            ev_dur = (pktimes{rec}{well}(:,3)-pktimes{rec}{well}(:,1))/fs;
            ED_flat{r_idx} = [ED_flat{r_idx} ev_dur'];
            ED_mean(r_idx,well-4) = mean(ev_dur);
            ED_std(r_idx,well-4) = std(ev_dur);
            
            ED_rms(r_idx,well-4) = mean(ev_dur.^2).^0.5;
            ED_50(r_idx, well-4) = quantile(ev_dur, 0.5);
            ED_05(r_idx, well-4) = quantile(ev_dur, 0.05);
            ED_95(r_idx, well-4) = quantile(ev_dur, 0.95);            
        end
    end
end