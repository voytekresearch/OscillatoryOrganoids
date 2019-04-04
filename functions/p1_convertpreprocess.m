function p1_convertpreprocess(batch_folder)
% Conversion and Preprocessing
% ** requires the Axion functions also included **
% 
% Each Axion .raw file in this folder is decoded and saved into a folder,
% where each well is saved as its own .mat file.
%
% For preprocessing, MEA_process does the following:
%   - downsamples raw data to 1000Hz LFP
%   - detect MUA/spikes based on adapted standard deviation estimation
%   from Quiroga et al., Neural Comp (2004)
%   - save LFP, spike timing, and spike waveforms to LFP_sp.mat
%
% For summary, MEA_summary computes and saves (& optionally prints out):
%   - median and welch's PSD
%   - network spike vector
%   - network spiking autocorrelation
% for all the wells, 

% define raw data folder to run conversion
%raw_data_folder = '/Users/rdgao/Documents/data/Lipton/MEA/';
cd(batch_folder)
raw_data_folder = batch_folder;
raw_files = dir('*.raw');
wells2process=1:12; % which of the 12 wells in 12-well plate to analyze

for f=1:length(raw_files)
    output_folder = [raw_data_folder raw_files(f).name(1:end-4)];
    processed_file = [output_folder '/LFP_Sp.mat'];
    disp(raw_files(f).name)
    % do the conversion
    MEA_convert(raw_files(f).name, output_folder, wells2process);
    % do preprocessing
    MEA_process(output_folder, wells2process);
    % compute and saveout summary info
    MEA_summary(processed_file, [output_folder '/'], wells2process, 1)
end