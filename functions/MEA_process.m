function MEA_process(data_folder, wells2process)
%   MEA_process(data_folder, wells2process)
%       data_folder: file path of folder where raw .mat well data is stored
%       wells2process: wells to process further
%   MEA_process takes groups of raw time series and does the following:
%       - downsample to get LFP at 1000Hz
%       - detect and grab spikes
%   really this could be applied to high sampling rate data in general, but
%   for now we keep it to MEA data

%analysis parameters-------------------------------------------------------
% currently only works for 12 well plates
num_well = 12;
num_chan = 64;
fs = 12500;
fs_ds = 1000;
sp_freq = [300,3000];
auto_spike = 1;
std_thr = 5.5; %5.5 std for spike detection
spike_len = 12;
wells = wells2process;

%LFP, spike detection & MUA-----------------------------------------------------
spikes = cell(num_well,num_chan);
spike_cnt = zeros(num_well,num_chan);
spike_shape = cell(num_well,num_chan);
spike_avg = zeros(num_well,num_chan,spike_len*2+1);
LFP = cell(num_well,1);

%filtering
for well=wells    
    f=[data_folder, '/well_' num2str(well) '.mat'];
    disp(['Loading... ',f])
    if exist(f, 'file')        
        load(f)
        
        % downsampling LFP to 1000Hz
        disp('Downsampling... ')
        LFP{well} = resample(MEA,2,25);
        
        % spike detection
        % first apply median well filter
        % find all 0 channels, exclude those in median calc
        mask = 1-all(MEA==0);
        well_med = median(MEA(:,find(mask)),2);
        % re-reference only the channels that have data
        MEA = MEA - well_med*mask;
        
        disp('Filtering... ')
        % filtering and do spike detection
        filtered = butterpass(MEA,fs,sp_freq,3);
        try
            [spikes(well,:), spike_cnt(well,:)] = spike_detect_abs(filtered,fs,std_thr);
        catch
            keyboard
        end
        
        % grab spike waveforms and compute average spike shape
        for chan = 1:num_chan
            spike_shape{well,chan} = collect_spikes(filtered,[],spikes{well,chan},spike_len);
            %.*repmat(-sign(spike_shape{well,chan}(spike_len+1,:)),spike_len*2+1,1)
            
            if ~isempty(spike_shape{well,chan})
                % get average spike waveform
                % spike_avg(well,chan,:) = mean(spike_shape{well,chan},2);
                
                % little matrix trick to get the average: taking the inner
                % product between spike shape and its polarity flips all
                % the spikes to the same direction
                spike_avg(well,chan,:) = -sign(spike_shape{well,chan}(spike_len+1,:))*spike_shape{well,chan}';
            end
        end   
        
    else
        disp('Non-existent well, skip.')
    end
end

t_s = t;
t_ds = (t(1):(1/fs_ds):t(end))';
final_output = [data_folder, '/LFP_Sp.mat'];
save(final_output, 'LFP', 't_ds', 'fs_ds', 'spikes', 'spike_shape', 'spike_avg', 'spike_cnt', 't_s', '-v7.3')
