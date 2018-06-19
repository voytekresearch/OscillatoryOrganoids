classdef BlockVectorData < Entry
    %BlockVectorData contains instructions for loading the Data types 
    %   (See BlockVectorDataType.m) from the data portions of the file
    %   listed in the header.
    
    properties(GetAccess = protected, SetAccess = private)
        FileID % File Handle to be used for file reading operations.
    end
    
    methods
        function this = BlockVectorData(aEntryRecord, aFileID)
            %BlockVectorData: Constructs a new BlockVectorData corresponding  
            % to an Entry Record and the file handle it came from
            this = this@Entry(aEntryRecord, int64(ftell(aFileID)));
            
            this.FileID = aFileID;
            fseek(this.FileID, double(this.EntryRecord.Length), 'cof');
            
            if ~(ftell(this.FileID) == (this.Start + this.EntryRecord.Length) || isinf(this.EntryRecord.Length))
                error('Unexpected BlockVectorHeader length')
            end
            
        end
        
        function Waveforms = GetRawV1ContinuousWaveforms(...
                this, ...
                aSourceSet, ...
                aChannelsToLoad, ...
                aTimeRange, ...
                aDimensions)
            
            if aSourceSet.Header.NumSamplesPerBlock ~= 1
                error('load_AxIS_raw:argNumSamplesPerBlock', ...
                    'Invalid header for RAW file: incorrect samples per block');
            end
            
            if aSourceSet.Header.BlockHeaderSize ~= 0
                error('load_AxIS_raw:argBlockHeaderSize', ...
                    'Invalid header for RAW file: incorrect block header size');
            end
            fseek(this.FileID, int64(this.Start), 'bof');
            fStart = 0;
            
            % Read data
            fSampleFreq =  int64(aSourceSet.Header.SamplingFrequency);
            fChannelCount = int64(aSourceSet.Header.NumChannelsPerBlock);
            fBytesPerSecond = fSampleFreq * fChannelCount * int64(2);
            fMaxTime = double(this.EntryRecord.Length) / double(fBytesPerSecond);
            
            if ~strcmp(aTimeRange, 'all')
                fStart = aTimeRange(1);
                fEnd =aTimeRange(2);
                
                if(fStart >= fEnd)
                    warning('Invalid timespan argument: end time < start time. No valid waveform can be returned.');
                    Waveforms = Waveform.empty;
                    return;
                end      
                
                fSkipInitialSamples = fStart * fSampleFreq;                
                fSkipInitialBytes = fSkipInitialSamples * fChannelCount * int64(2);
          
                if(fStart > fMaxTime)
                    warning( 'DataSet only contains %d Seconds of data (%d Seconds was the requested start time). No waveform will be returned.', fMaxTime, fStart);
                    Waveforms = Waveform.empty;
                    return;
                end
                
                if(fEnd > fMaxTime)
                   fEnd = fMaxTime; 
                   warning('DataSet only contains %d Seconds of data. Returned waveforms will be shorter than requested (%d Seconds).', fMaxTime, fEnd - fStart);
                end
                
                fNumSamples  = int64((fEnd - fStart) * aSourceSet.Header.SamplingFrequency);
                
                % skip past samples that are before the current time range
                fseek(this.FileID, fSkipInitialBytes, 'cof');
            else
                fNumSamples =  int64((fMaxTime) * aSourceSet.Header.SamplingFrequency);
            end
            
            % Read the data for the given channel, skipping from one to the next
            fNumChannels = length(aSourceSet.ChannelArray.Channels);
            
            fMaxExtents = PlateTypes.GetElectrodeDimensions(aSourceSet.ChannelArray.PlateType);
            
            if (isempty(fMaxExtents))
                fMaxExtents = [max([aSourceSet.ChannelArray.Channels(:).WellRow]), ...
                    max([aSourceSet.ChannelArray.Channels(:).WellColumn]), ...
                    max([aSourceSet.ChannelArray.Channels(:).ElectrodeColumn]), ...
                    max([aSourceSet.ChannelArray.Channels(:).ElectrodeRow])];
            end
                
            switch aDimensions
                    case {LoadArgs.ByPlateDimensions}
                        Waveforms = [];
                    case {LoadArgs.ByWellDimensions}
                        Waveforms = cell(double(fMaxExtents(1:2)));
                    case {LoadArgs.ByElectrodeDimensions}
                        Waveforms = cell(double(fMaxExtents));
            end
            
            if length(aChannelsToLoad) == 1
                % We're only reading one channel. For efficiency, take advantage of fread's
                % 'skip' argument.
                
                fread(this.FileID, aChannelsToLoad - 1, 'int16=>int16');
                fWaveform = fread(this.FileID, fNumSamples, 'int16=>int16', 2*(fNumChannels - 1));
                fChannelMapping = aSourceSet.ChannelArray.Channels(aChannelsToLoad);
                
                fWaveform = Waveform(fChannelMapping, fStart, fWaveform, aSourceSet);
               
                fOuputIndex = double([...
                        fChannelMapping.WellRow, ...
                        fChannelMapping.WellColumn, ...
                        fChannelMapping.ElectrodeColumn, ...
                        fChannelMapping.ElectrodeRow]);
                    
                switch aDimensions
                        
                        case {LoadArgs.ByPlateDimensions}
                            Waveforms = fWaveform;
                                                       
                        case {LoadArgs.ByWellDimensions}
                            Waveforms{fOuputIndex(1),fOuputIndex(2)} = fWaveform;
                            
                        case {LoadArgs.ByElectrodeDimensions}
                            Waveforms{fOuputIndex(1),fOuputIndex(2), fOuputIndex(3),fOuputIndex(4)} = fWaveform;
                   
                end
            else
                
                fNumSamples      = fNumSamples * fNumChannels;
                fTempChannelData = fread(this.FileID, fNumSamples, 'int16=>int16');
                
                % This test (which can fail only when we read the the end of the file)
                % makes sure that we didn't get a number of samples that's not divisible
                % by the number of channels.
                fRemainderCount = mod(length(fTempChannelData), fNumChannels);
                if fRemainderCount ~= 0
                    warning('load_AXiS_raw:remainderCheck', ...
                        'This Data has the wrong number of samples for %u channels, File may be corrupt', ...
                        fNumChannels);
                    fNumSamples = int64((length(fTempChannelData)/ fNumChannels) - 1) * fNumChannels;
                    fTempChannelData = fTempChannelData(1:fNumSamples);
                end
                
                % Convert the 1D array to a 2D array, with channel as the second dimension
                % (starts as first dimension and then is transposed)
                fTempChannelData = reshape(fTempChannelData(:), fNumChannels, []);
                
                for fChannelIndex = aChannelsToLoad
                    
                    fChannelMapping = aSourceSet.ChannelArray.Channels(fChannelIndex);
                    fWaveform = fTempChannelData(fChannelIndex,:)';
                    fWaveform = Waveform(fChannelMapping, fStart, fWaveform, aSourceSet);
                    
                    fOuputIndex = double([...
                        fChannelMapping.WellRow, ...
                        fChannelMapping.WellColumn, ...
                        fChannelMapping.ElectrodeColumn, ...
                        fChannelMapping.ElectrodeRow]);
                    
                    switch aDimensions
                        
                        case {LoadArgs.ByPlateDimensions}
                            if(isempty(Waveforms))
                                Waveforms = fWaveform;
                            else
                                Waveforms(length(Waveforms) + 1) = fWaveform;
                            end
                            
                        case {LoadArgs.ByWellDimensions}
                            Waveforms{fOuputIndex(1),fOuputIndex(2)}(...
                                length(Waveforms{fOuputIndex(1),fOuputIndex(2)}) + 1) = fWaveform;
                            
                        case {LoadArgs.ByElectrodeDimensions}
                            Waveforms{fOuputIndex(1),fOuputIndex(2), fOuputIndex(3),fOuputIndex(4)} = fWaveform; %We only expect one Waveform per channel
                   
                    end
                end
            end
        end
        
        function Waveforms = GetSpikeV1Waveforms(...
                this, ...
                aSourceSet, ...
                aChannelsToLoad, ...
                aTimeRange, ...
                aDimensions)
            
            fStorageType = 0;
            
            switch aDimensions
                case {LoadArgs.ByWellDimensions}
                    fStorageType = 1;
                case {LoadArgs.ByElectrodeDimensions}
                    fStorageType = 2;
            end
            
            if aSourceSet.Header.NumChannelsPerBlock ~= 1
                error('Load_spike_v1:argNumChannelsPerBlock', ...
                    'Invalid header for SPIKE file: incorrect channels per block');
            end
            
            if aSourceSet.Header.BlockHeaderSize < Spike_v1.LOADED_HEADER_SIZE
                error('Load_spike_v1:argBlockHeaderSize', ...
                    'Invalid header for SPIKE file: block header size too small');
            end
            
            if aSourceSet.Header.NumSamplesPerBlock < 1
                error('load_AxIS_spike:argNumSamplesPerBlock', ...
                    'Invalid header for SPIKE file: number of samples per block < 1');
            end
            
            fHeader = aSourceSet.Header;
            
            if ischar(aTimeRange) && strcmp(aTimeRange, 'all')
                fFirstSample = 0;
                fLastSample  = Inf;
            else
                fFirstSample = aTimeRange(1) * fHeader.SamplingFrequency;
                fLastSample  = aTimeRange(2) * fHeader.SamplingFrequency;
            end
            
            fMaxExtents = PlateTypes.GetElectrodeDimensions(aSourceSet.ChannelArray.PlateType);
            
            if (isempty(fMaxExtents))
                fMaxExtents = [max([aSourceSet.ChannelArray.Channels(:).WellRow]), ...
                    max([aSourceSet.ChannelArray.Channels(:).WellColumn]), ...
                    max([aSourceSet.ChannelArray.Channels(:).ElectrodeColumn]), ...
                    max([aSourceSet.ChannelArray.Channels(:).ElectrodeRow])];
            end
            
            if (fStorageType == 0)
                Waveforms = [];
            elseif (fStorageType == 1)
                Waveforms = cell(double(fMaxExtents(1:2)));
            elseif fStorageType == 2
                Waveforms = cell(double(fMaxExtents));
            end
            
            fDesiredChannelsLut = zeros(length(aSourceSet.ChannelArray.Channels),1);
            
            fDesiredChannelsLut(aChannelsToLoad) = 1;
            
            fWaveformBytesSize = (fHeader.NumChannelsPerBlock * fHeader.NumSamplesPerBlock * 2) + fHeader.BlockHeaderSize;
            fseek(this.FileID, int64(this.Start), 'bof');
            fData = fread(this.FileID, this.EntryRecord.Length, 'int8=>int8');
            fData = reshape(fData, fWaveformBytesSize, []);
            fNumWaves = size(fData,2);
            fNumWaves = fNumWaves(1);
            fChannelArray = aSourceSet.ChannelArray;
            for iWave = 1 : fNumWaves
                fStartingSample = typecast(fData(1:8,iWave), 'int64');
                
                if fStartingSample < fFirstSample
                    % Too early
                    continue;
                elseif fStartingSample >= fLastSample
                    % Too late.  Since spikes are in chronological order, we
                    % can stop reading now.
                    break;
                end
                
                fHardwareChannelIndex = typecast(fData(9,iWave), 'uint8');
                fHardwareChannelAchk = typecast(fData(10,iWave), 'uint8');
                fChannelIndex = fChannelArray.LookupChannel(fHardwareChannelAchk, fHardwareChannelIndex);
                fChannelMapping = fChannelArray.Channels(fChannelIndex);
                
                fOuputIndex = double([...
                    fChannelMapping.WellRow, ...
                    fChannelMapping.WellColumn, ...
                    fChannelMapping.ElectrodeColumn, ...
                    fChannelMapping.ElectrodeRow]);
                
                if fDesiredChannelsLut(fChannelIndex) ~= 0
                    
                    fWaveform = Spike_v1( ...
                        fChannelMapping, ...
                        double(fStartingSample) / fHeader.SamplingFrequency, ...
                        typecast(fData(fHeader.BlockHeaderSize+1:fWaveformBytesSize,iWave)','int16')', ...
                        aSourceSet, ...
                        typecast(fData(11:14,iWave), 'uint32'), ...
                        typecast(fData(15:22,iWave), 'double'), ...
                        typecast(fData(23:30,iWave), 'double'));
                    
                    if (fStorageType == 2)
                        Waveforms{fOuputIndex(1),fOuputIndex(2), fOuputIndex(3),fOuputIndex(4)}(...
                            length(Waveforms{fOuputIndex(1),fOuputIndex(2), fOuputIndex(3),fOuputIndex(4)}) + 1) = fWaveform;

                    elseif (fStorageType == 0)
                        if(isempty(Waveforms))
                            Waveforms = fWaveform;
                        else
                            Waveforms(length(Waveforms) + 1) = fWaveform;
                        end
                            
                    elseif (fStorageType == 1)
                        Waveforms{fOuputIndex(1),fOuputIndex(2)}(...
                            length(Waveforms{fOuputIndex(1),fOuputIndex(2)}) + 1) = fWaveform;
                   end
                end
            end
        end
        
        function [aElectrodes, aTimes] = GetAllSpikeTimes(this, aSourceSet)
            if aSourceSet.Header.NumChannelsPerBlock ~= 1
                error('Load_spike_v1:argNumChannelsPerBlock', ...
                    'Invalid header for SPIKE file: incorrect channels per block');
            end
            
            if aSourceSet.Header.BlockHeaderSize < Spike_v1.LOADED_HEADER_SIZE
                error('Load_spike_v1:argBlockHeaderSize', ...
                    'Invalid header for SPIKE file: block header size too small');
            end
            
            if aSourceSet.Header.NumSamplesPerBlock < 1
                error('load_AxIS_spike:argNumSamplesPerBlock', ...
                    'Invalid header for SPIKE file: number of samples per block < 1');
            end
            
            fHeader = aSourceSet.Header;
            
            fWaveformBytesSize = (fHeader.NumChannelsPerBlock * fHeader.NumSamplesPerBlock * 2) + fHeader.BlockHeaderSize;
            fseek(this.FileID, int64(this.Start), 'bof');
            fData = fread(this.FileID, this.EntryRecord.Length, 'int8=>int8');
            fData = reshape(fData, fWaveformBytesSize, []);
            fNumWaves = size(fData,2);
            fNumWaves = fNumWaves(1);
            
            %pre-allocate data
            aTimes(fNumWaves) = 0;
            
            %copy chip and channel data
            aElectrodes.Achk = fData(10,:);
            aElectrodes.Channel = fData(9,:);
            
            %calculate spike times
            for iWave = 1 : fNumWaves
                fStartingSample = double(typecast(fData(1:8,iWave), 'int64'));
                fSampleOffset = double(typecast(fData(11:14,iWave), 'uint32'));
                
                aTimes(iWave) = (fStartingSample + fSampleOffset) / fHeader.SamplingFrequency;
            end
        end
    end
end
