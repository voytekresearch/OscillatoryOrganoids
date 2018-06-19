  %
  %  calc_spike_rate(aInputFile, aOutputFile, aIntervalSeconds, aStartSeconds, aEndSeconds)
  %   
  %   Required Parameters
  %     aInputFile        path and name of the input spike file
  %     aOutputFile       path and name of the spike rate file to write (in CSV format)
  %
  %   Optional Parameters
  %     aStartSeconds     time you want to start from in the file, in seconds (default 0)
  %     aEndSeconds       time you want to stop, in seconds (default Inf, end of file)
  %     aIntervalSeconds  interval that will be used between aStartSeconds and aEndSeconds (default 60)
  %
  %   Note: Parameters must be specified in order.  If optional parameters are omitted,
  %         the default is used, as follows:
  %
  %     calc_spike_rate(aInputFile, aOutputFile)
  %         => uses defaults for aIntervalSeconds, aStartSeconds, aEndSeconds
  %
  %     calc_spike_rate(aInputFile, aOutputFile, aIntervalSeconds)
  %         => uses defaults for StartSeconds, aEndSeconds
  %
  %     calc_spike_rate(aInputFile, aOutputFile, aIntervalSeconds, aStartSeconds)
  %         => uses defaults for aEndSeconds
  %
  %
  %  The script will divide the time between aStartSeconds and aEndSeconds into N intervals,
  %  where N = aEndSeconds - aStartSeconds) / aIntervalSeconds, rounded up, and output
  %  spike counts for the intervals:
  %
  %  [aStartSeconds,                           aStartSeconds +   aIntervalSeconds)
  %  [aStartSeconds +       aIntervalSeconds,  aStartSeconds + 2*aIntervalSeconds)
  %  [aStartSeconds +     2*aIntervalSeconds,  aStartSeconds + 3*aIntervalSeconds)
  %  ....
  %  [aStartSeconds + (N-1)*aIntervalSeconds, aStartSeconds +  N*aIntervalSeconds)
  %

function calc_spike_rate(aInputFile, aOutputFile, aIntervalSeconds, aStartSeconds, aEndSeconds)

    ELECTRODE_ROWS = 8;
    ELECTRODE_COLS = 8;
    NUM_ELECTRODES = ELECTRODE_ROWS * ELECTRODE_COLS;

    % Handle default arguments:
    % aIntervalSeconds = 60
    % aStartSeconds = 0
    % aEndSeconds = Inf
    if nargin < 3
        aIntervalSeconds = 60;
    end
    
    if nargin < 4
        aStartSeconds = 0;
    end
    
    if nargin < 5
        aEndSeconds = Inf;
    end
    
    if nargin > 5
        error('calc_spike_rate:tooManyParameters', ['calc_spike_rate: Too many input parameters']);
    end

    if ischar(aInputFile)
        disp(['Loading spikes from file ' aInputFile]);
        
        fSpikeData = load_AxIS_file(aInputFile, [aStartSeconds aEndSeconds]);
    else
        error('calc_spike_rate:invalidFileName', ['calc_spike_rate: Invalid file name: ' aInputFile]);
    end
    
    if ~strcmp(fSpikeData.fileType, 'spike') 
        if strcmp(fSpikeData.fileType, 'spike-beta')
            % This is a beta-format spike file.  Not currently supported by calc_spike_rate.
            error('calc_spike_rate:betaFormat', ['Beta-format spike file not supported: ' aInputFile]);
        else
            % File was not a spike file
            error('calc_spike_rate:invalidSpikeFile', ['Invalid spike file format in file ' aInputFile]);
        end
    end
    
    % If aEndSeconds is infinite, we want to read until the end of the file, but we can't leave
    % it as Inf, because then we'll try to create a table with an infinite number of rows.  Instead,
    % base it on the last spike in the file.
    if aEndSeconds == Inf
        aEndSeconds = (fSpikeData.spikes(end).startingSample + fSpikeData.spikes(end).triggerSampleOffset) / ...
                      fSpikeData.samplingFrequency;  
    end
    
    fNumIntervals = ceil((aEndSeconds - aStartSeconds) / aIntervalSeconds);
    
    fSpikeCounts = zeros(fNumIntervals, ELECTRODE_ROWS, ELECTRODE_COLS);
    
    for i=1:length(fSpikeData.spikes)
        % Iterate over each loaded spike and determine whether or where it falls in the spike
        % rate table.
        fSpikeTime = (fSpikeData.spikes(i).startingSample + fSpikeData.spikes(i).triggerSampleOffset) / ...
                     fSpikeData.samplingFrequency;
        fSpikeInterval = floor((fSpikeTime - aStartSeconds) / aIntervalSeconds);
        
        if fSpikeInterval < 0
            % This spike was too early to be counted
            continue;
        elseif fSpikeInterval >= fNumIntervals
            % This spike was too late to be counted.  Since the file is ordered by spike time,
            % we can stop reading.
            break;
        else
            % we have the channel number ranging from 11 to 88 (this is more a label than a number)
            % Decode it:
            fChannelY = floor(fSpikeData.spikes(i).channel / 10);
            fChannelX = mod(fSpikeData.spikes(i).channel, 10);
            %fChannelIndex = fChannelX + (fChannelY - 1) * ELECTRODE_COLS;

            fSpikeCounts(fSpikeInterval + 1, fChannelY, fChannelX) =  ...
                    fSpikeCounts(fSpikeInterval + 1, fChannelY, fChannelX) + 1;
        end
    end
    
    % Make sure file doesn't exist so we don't overwrite it
    if exist(aOutputFile, 'file')
        error('calc_spike_rate:outputFileExists', ['Output file already exists: ' aOutputFile]);
    end
    
    % Write spike rates in CSV format
    fid = fopen(aOutputFile, 'w');
    if fid == -1
        % error opening output file
        error('calc_spike_rate:errorOpeningOutputFile', ['Couldn''t write to output file ' aOutputFile]);
    end
    
    fprintf(fid, 'Interval Start (S),Interval End (S)');
    for i=1:ELECTRODE_ROWS
        for j=1:ELECTRODE_COLS
            fprintf(fid, ',Channel %d%d', i, j);
        end
    end
    
    fprintf(fid, '\n');
    
    for fInterval=1:fNumIntervals
        fprintf(fid, '%f,%f', aStartSeconds + (fInterval - 1) * aIntervalSeconds, ...
                              aStartSeconds + fInterval * aIntervalSeconds);
                              
        for i=1:ELECTRODE_ROWS
            for j=1:ELECTRODE_COLS
                fprintf(fid, ',%d', fSpikeCounts(fInterval, i, j));
            end
        end
        
        fprintf(fid, '\n');
    end
    
    fclose(fid);
    
    disp(['Finished writing ' aOutputFile]);

