classdef BlockVectorLegacyLoader < BlockVectorData
    %BLOCKVECTORLEGACYLOADER Static class with resources to load file data
    %   to legacy style structures.
    
    properties(GetAccess = private, Constant = true)
        % Fixed by file format - there's only 1 byte for each
        MAX_CHANNEL_ARTICHOKE = 256;
        MAX_CHANNEL_INDEX     = 256;
    end
    
    methods (Static = true)
        function aData = Legacy_Load_Raw_v1(aSourceData, aHeader, aTimeRange)
            
            % Check to make sure the file format makes sense for a raw file
            if aHeader.header.numSamplesPerBlock ~= 1
                error('load_AxIS_raw:argNumSamplesPerBlock', ...
                    'Invalid header for RAW file: incorrect samples per block');
            end
            
            if aHeader.header.blockHeaderSize ~= 0
                error('load_AxIS_raw:argBlockHeaderSize', ...
                    'Invalid header for RAW file: incorrect block header size');
            end
            
            fseek(aSourceData.FileID, int64(aSourceData.Start), 'bof');
            
            % Copy input structure to output; then add information to output
            aData = aHeader;
            aData.numChannels = aData.header.numChannelsPerBlock;
            if aData.numChannels ~= aData.channelArray.numChannels
                error('load_AxIS_raw:mismatchNumChannels', ...
                    'Invalid RAW file: mismatched number of channels');
            end
            
            if isempty(aHeader.loadedChannels)
                % No channels requested - don't read any data
                return;
            end
            
            % Read data
            fSampleFreq =  int64(aData.samplingFrequency);
            fChannelCount = int64(aData.numChannels);
            fBytesPerSecond = fSampleFreq * fChannelCount * int64(2);
            fMaxTime = int64(aSourceData.EntryRecord.Length / fBytesPerSecond);
            if ~strcmp(aTimeRange, 'all')
                fStart = aTimeRange(1);
				fEnd =aTimeRange(2);
                
                if(fStart >= fEnd)
                    warning('Invalid timespan argument: end time < start time. No valid waveform can be returned.');
                    return;
                end
                if(fStart > fMaxTime)
                    warning( 'DataSet only contains %d Seconds of data (%d Seconds was the requested start time). No waveform will be returned.', fMaxTime, fStart);
                    return;
                end
                if(fEnd > fMaxTime)
                   fEnd = fMaxTime; 
                   warning('DataSet only contains %d Seconds of data. Returned waveforms will be shorter than requested (%d Seconds).', fMaxTime, fEnd - fStart);
                end
                
                fSkipInitialSamples = fStart * fSampleFreq; 
                fSkipInitialBytes = fSkipInitialSamples * fChannelCount * int64(2);
                
                fNumSamples  = int64((fEnd - fStart) * fSampleFreq);
                
				fseek(aSourceData.FileID, fSkipInitialBytes, 'cof');				
            else
                fNumSamples = int64(fMaxTime * fSampleFreq);
            end
            
            if length(aHeader.loadedChannels) == 1
                % We're only reading one channel. For efficiency, take advantage of fread's
                % 'skip' argument.
                
                % First, read and throw away other channels at the beginning of the data
                fread(aSourceData.FileID, aHeader.loadedChannels - 1, 'int16=>int16');
                
                % Read the data for the given channel, skipping from one to the next
                aData.channelData = fread(aSourceData.FileID, fNumSamples, 'int16=>int16', 2*(aData.numChannels - 1));
            else
                % We're reading multiple channels.  There's no easy way to translate this into a
                % single fread, so we read all of the channels in the sample range in question,
                % then remap, throwing away what we didn't want.
                
                % First, read everything within the time range
                fNumSamples      = fNumSamples * int64(aData.numChannels);
                fTempChannelData = fread(aSourceData.FileID, fNumSamples, 'int16=>int16');
                
                % This test (which can fail only when we read the the end of the file)
                % makes sure that we didn't get a number of samples that's not divisible
                % by the number of channels.
                fRemainderCount = mod(length(fTempChannelData), aData.numChannels);
                if fRemainderCount ~= 0
                    warning('load_AXiS_raw:remainderCheck', ...
                        'File %s has wrong number of samples for %u channels', ...
                        aSourceData.FileID, aData.numChannels);
                end
                
                % Convert the 1D array to a 2D array, with channel as the second dimension
                % (starts as first dimension and then is transposed)
                fTempChannelData = reshape(fTempChannelData(:), aData.numChannels, []);
                
                % Remap
                aData.channelData = fTempChannelData(aData.loadedChannels, :);
                aData.channelData = aData.channelData';
            end
        end
        
        function aData = Legacy_Load_spike_v1(aSourceData, aHeader, aTimeRange)
            % Check to make sure the file format makes sense for a spike file
            if aHeader.header.numChannelsPerBlock ~= 1
                error('Load_spike_v1:argNumChannelsPerBlock', ...
                    'Invalid header for SPIKE file: incorrect channels per block');
            end
            
            if aHeader.header.blockHeaderSize < Spike_v1.LOADED_HEADER_SIZE
                error('Load_spike_v1:argBlockHeaderSize', ...
                    'Invalid header for SPIKE file: block header size too small');
            end
            
            if aHeader.header.numSamplesPerBlock < 1
                error('load_AxIS_spike:argNumSamplesPerBlock', ...
                    'Invalid header for SPIKE file: number of samples per block < 1');
            end
            
            % Copy input structure to output; then add information to output
            aData = aHeader;
            
            if isempty(aHeader.loadedChannels)
                % No channels requested - don't read any data
                return;
            end
            
            if ischar(aTimeRange) && strcmp(aTimeRange, 'all')
                fFirstSample = 0;
                fLastSample  = Inf;
            else
                fFirstSample = aTimeRange(1) * aData.samplingFrequency;
                fLastSample  = aTimeRange(2) * aData.samplingFrequency;
            end
            
            fSparseChannelIndex = sparse(...
                double([aData.channelArray.channel.channelAchk] + 1) , ...
                double([aData.channelArray.channel.channelIndex] + 1), ...
                double(1:aData.channelArray.numChannels), ...
                BlockVectorLegacyLoader.MAX_CHANNEL_ARTICHOKE,...
                BlockVectorLegacyLoader.MAX_CHANNEL_INDEX);
            
            spikeIndex = 1;
            
            % Manage the spike array to avoid "growing" it
            fSpikeArrayCapacity = 0;
            
            fseek(aSourceData.FileID, int64(aSourceData.Start), 'bof');
            fStop = aSourceData.Start + aSourceData.EntryRecord.Length;
            
            while ftell(aSourceData.FileID) < fStop
                % Load the block header for this spike
                fCurrentSpike = [];
                fCurrentSpike.startingSample      = fread(aSourceData.FileID, 1, 'int64'); %Changed uint64 to int64
                
                if feof(aSourceData.FileID)
                    % we were at the end of the file
                    break;
                end
                
                fHardwareChannelIndex              = fread(aSourceData.FileID, 1, 'uint8');
                fHardwareChannelAchk               = fread(aSourceData.FileID, 1, 'uint8');
                fCurrentSpike.triggerSampleOffset  = fread(aSourceData.FileID, 1, 'uint32');
                fCurrentSpike.stDev                = fread(aSourceData.FileID, 1, 'double');
                fCurrentSpike.threshold            = fread(aSourceData.FileID, 1, 'double');
                
                % Get well and electrode row and column
                fCurrentSpike.channelArrayIndex = full(fSparseChannelIndex(fHardwareChannelAchk + 1, fHardwareChannelIndex + 1));
                if fCurrentSpike.channelArrayIndex == 0
                    error('load_AxIS_spike:invalidSpikeChannel', ...
                        ['Spike has invalid Artichoke/Channel number ' num2str(fHardwareChannelAchk) ...
                        ' / ' num2str(fHardwareChannelIndex)]);
                end
                
                fCurrentSpike.channelInfo = aData.channelArray.channel(fCurrentSpike.channelArrayIndex);
                
                % Historical channel number -- XY, where X = electrodeColumn and Y = electrodeRow
                % This field is deprecated, since it doesn't contain multiwell information.
                fCurrentSpike.channel = 10*fCurrentSpike.channelInfo.electrodeColumn + ...
                    fCurrentSpike.channelInfo.electrodeRow;
                
                % Seek forward to the end of the header.
                % It might seem cleaner to do this by saving the current position
                % before we read the header and seeking forward from there.
                % However, that method will not handle files larger than 4GB if fseek
                % and ftell are not 64-bit clean.
                fseek(aSourceData.FileID, aData.header.blockHeaderSize - Spike_v1.LOADED_HEADER_SIZE, 'cof');
                
                if ftell(aSourceData.FileID) >= fStop
                    % file ended before the start of the data -- looks like a truncated file
                    warning('load_AxIS_spike:headerTruncated', ...
                        'Spike header truncated at end of file');
                    break;
                end
                
                % Read the data for this spike
                fCurrentSpike.waveform = fread(aSourceData.FileID, aData.header.numSamplesPerBlock, 'int16=>int16');
                
                if length(fCurrentSpike.waveform) < aData.header.numSamplesPerBlock
                    % file ended in the middle of the data
                    warning('load_AxIS_spike:dataTruncated', ...
                        'Spike data truncated at end of file');
                    break;
                end
                
                % Does this spike match the channel filter?  If not, skip it.
                if ~ismember(fCurrentSpike.channelArrayIndex, aData.loadedChannels)
                    continue;
                end
                
                % Does this spike match the time filter?
                if fCurrentSpike.startingSample < fFirstSample
                    % Too early
                    continue;
                elseif fCurrentSpike.startingSample >= fLastSample
                    % Too late.  Since spikes are in chronological order, we
                    % can stop reading now.
                    break;
                end
                
                
                % Check for capacity in the spike array
                if spikeIndex > fSpikeArrayCapacity
                    if fSpikeArrayCapacity == 0
                        % Let it grow on its own
                        fSpikeArrayCapacity = 1;
                    else
                        % Double the capacity
                        fSpikeArrayCapacity = fSpikeArrayCapacity * 2;
                        aData.spikes(fSpikeArrayCapacity).startingSample = -1;
                    end
                end
                
                % save the spike
                aData.spikes(spikeIndex) = fCurrentSpike;
                spikeIndex = spikeIndex + 1;
            end
            
            % trim down the spike array
            if spikeIndex > 1
                aData.spikes = aData.spikes(1:spikeIndex-1);
            else
                aData.spikes = [];
            end
        end
    end
    
end

