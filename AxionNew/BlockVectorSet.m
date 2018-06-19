classdef BlockVectorSet < handle
    %BLOCKVECTORSET is a grouping of data and metadata for a series of data contained within an AxisFile.
    %   this class is composed of 4 major parts:
    %
    %   ChannelArray:    The channel array (See ChannelArray.m) is a listing of
    %                    all of the channels that were recorded in this
    %                    loaded set.
    %
    %   Header:          The header of the data (See BlockVectorHeader.m) contains the basic infomation
    %                    that is used in loading and using the data in this
    %                    set (e.g. Sampling Frequency, Voltage Scale, etc...)
    %
    %   HeaderExtension: The header extension (See BlockVectorHeaderExtension.m)
    %                    contains metadata about the data capture / Reprocessing.
    %
    %   Data:            The data in this file (See BlockVectorData.m)
    %                    contains the methods for loading the sample data
    %                    from this set.
    
    properties (SetAccess = private, GetAccess = private)
        sourceFile
    end
    
    properties(GetAccess = public, SetAccess = private)
        ChannelArray
        Header
        HeaderExtension
        Data
    end
    
    methods
        
        function this = BlockVectorSet(varargin)
            this.setval(varargin{:});
        end
        
        
        function clone = Clone(this, varargin)
            clone = BlockVectorSet();
            if(isa(this.ChannelArray,'ChannelArray'))
                clone.ChannelArray = this.ChannelArray;
            end
            if(isa(this.Header,'BlockVectorHeader'))
                clone.Header = this.Header;
            end
            if(isa(this.HeaderExtension,'BlockVectorHeaderExtension'))
                clone.HeaderExtension = this.HeaderExtension;
            end
            if(isa(this.Data,'BlockVectorData'))
                clone.Data = this.Data;
            end
            if(isa(this.sourceFile,'AxisFile'))
                clone.sourceFile = this.sourceFile;
            end
            clone.setval(varargin{:});
        end
        
    end
    
    methods(Access = private)
        
        function setval(this, varargin)
            for i = 1:length(varargin)
                
                arg = varargin{i};
                
                if(isa(arg,'ChannelArray'))
                    this.ChannelArray = arg;
                elseif(isa(arg,'BlockVectorHeader'))
                    this.Header = arg;
                elseif(isa(arg,'BlockVectorHeaderExtension'))
                    this.HeaderExtension = arg;
                elseif(isa(arg,'BlockVectorData'))
                    this.Data = arg;
                elseif(isa(arg,'AxisFile'))
                    this.sourceFile = arg;
                else
                    error('Unknown member type');
                end
                
            end
        end
        
    end
    
    methods (Access = public)
        
        
        function aData = LoadData(this, varargin)
            % LoadData loads a Dataset, creaing a data structure similar
            % to the one created by load_AxIS_file (Deprecated)
            %
            %  Legal forms:
            %     data = LoadData();
            %     data = LoadData(well);
            %     data = LoadData(electrode);
            %     data = LoadData(well, electrode);
            %     data = LoadData(timespan);
            %     data = LoadData(well, timespan);
            %     data = LoadData(electrode, timespan);
            %     data = LoadData(well, electrode, timespan);
            %     data = LoadData(dimensions);
            %     data = LoadData(well, dimensions);
            %     data = LoadData(electrode, dimensions);
            %     data = LoadData(well, electrode, dimensions);
            %     data = LoadData(timespan, dimensions);
            %     data = LoadData(well, timespan, dimensions);
            %     data = LoadData(electrode, timespan, dimensions);
            %     data = LoadData(well, electrode, timespan, dimensions);
            %
            %  Optional arguments:
            %    well        String listing which wells (in a multiwell file) to load.
            %                Format is a comma-delimited string with whitespace ignored, e.g.
            %                'A1, B2,C3' limits the data loaded to wells A1, B2, and C3.
            %                Also acceptable: 'all' to load all wells.
            %                If this parameter is omitted, all wells are loaded.
            %                For a single-well file, this parameter is ignored.
            %
            %    electrode   Which electrodes to load.  Format is either a comma-delimited string
            %                with whitespace ignored (e.g. '11, 22,33') or a single channel number;
            %                that is, a number, not part of a string.
            %                Also acceptable: 'all' to load all channels and 'none', '-1', or -1
            %                to load no data (returns only header information).
            %                If this parameter is omitted, all channels are loaded.
            %
            %    timespan    Span of time, in seconds, over which to load data.  Format is a two-element
            %                array, [t0 t1], where t0 is the start time and t1 is the end time and both
            %                are in seconds after the first sample in the file.  Samples returned are ones
            %                that were taken at time >= t0 and <= t1.  The beginning of the file
            %                is at 0 seconds.
            %                If this parameter is omitted, the data is not filtered based on time.
            %
            %    dimensions  Preferred number of dimensions to report the waveforms in.
            %                Value must be a whole number scalar: 1, 3, or 5 (Other values are ignored):
            %
            %                dimensions = 1 -> ByPlate: returns a vector of Waveform objects, 1 Waveform
            %                                  per signal in the plate
            %
            %                dimensions = 3 -> ByWell: Cell Array of vectors of waveform 1 Waveform per signal
            %                                  in the electrode with size (well Rows) x (well Columns)
            %
            %                dimensions = 5 -> ByElectrode: Cell Array of vectors of waveform 1 Waveform per .
            %                                  signal in the electrode with size (well Rows) x (well Columns) x
            %                                  (electrode Columns) x (electrode Rows)
            %
            %                NOTE: The default loading dimensions for
            %                continous raw data is 5 and the default for
            %                spike data is 3.
            
            fLoadArgs = LoadArgs(varargin);
            fTargetWell = fLoadArgs.Well;
            fTargetElectrode = fLoadArgs.Electrode;
            
            fChannelsToLoad = BlockVectorSet.get_channels_to_load(this.ChannelArray, fTargetWell, fTargetElectrode);
            
            if(isempty(this.HeaderExtension))
                fDataType = this.sourceFile.PrimaryDataType;
            else
                fDataType = this.HeaderExtension.DataType;
            end
            
            if(fDataType == BlockVectorDataType.Raw_v1)
                fDimensions = fLoadArgs.Dimensions;
                if(isempty(fDimensions))
                    fDimensions = LoadArgs.ByElectrodeDimensions;
                end
                aData = this.Data.GetRawV1ContinuousWaveforms( ...
                    this, ...
                    fChannelsToLoad, ...
                    fLoadArgs.Timespan, ...
                    fDimensions);
                
            elseif(fDataType == BlockVectorDataType.Spike_v1)
                fDimensions = fLoadArgs.Dimensions;
                if(isempty(fDimensions))
                    fDimensions = LoadArgs.ByElectrodeDimensions;
                end
                aData = this.Data.GetSpikeV1Waveforms( ...
                    this, ...
                    fChannelsToLoad, ...
                    fLoadArgs.Timespan, ...
                    fDimensions);
                
            end
            
        end
        
        function [aElectrodes, aTimes] = LoadAllSpikes(this)
            % LoadAllSpikes attempts to load all spikes from file
            % 
            % Function returns an array of electrodes whre spikes occur
            % and array of times when spikes occur.  There is one entry for
            % every spike.
            %
            % aElectrodes : a structure with two elements (arrays of int8):  
            %    aElectrodes.Achk    - Artichoke chip where spike occurred
            %    aElectrodes.Channel - Artichoke's channel where spike
            %    occurred
            %  
            % aTimes : array of times (in seconds since file start) when spike occured
            %
            
            if(isempty(this.HeaderExtension))
                fDataType = this.sourceFile.PrimaryDataType;
            else
                fDataType = this.HeaderExtension.DataType;
            end
            
            if(fDataType == BlockVectorDataType.Spike_v1)
                [aElectrodes, aTimes] = this.Data.GetAllSpikeTimes(this);
            else
                error('File is not a spike file');
            end
        end
        
        function aData = load_as_legacy_struct(this, varargin)
            % load_AsLegacyStruct loads a dataset, creating a data structure similar
            % to the one created by load_AxIS_file (Deprecated)
            %
            %  Legal forms:
            %     data = load_as_legacy_struct();
            %     data = load_as_legacy_struct(well);
            %     data = load_as_legacy_struct(electrode);
            %     data = load_as_legacy_struct(well, electrode);
            %     data = load_as_legacy_struct(timespan);
            %     data = load_as_legacy_struct(well, timespan);
            %     data = load_as_legacy_struct(electrode, timespan);
            %     data = load_as_legacy_struct(well, electrode, timespan);
            %
            %
            %  Optional arguments:
            %    well        String listing which wells (in a multiwell file) to load.
            %                Format is a comma-delimited string with whitespace ignored, e.g.
            %                'A1, B2,C3' limits the data loaded to wells A1, B2, and C3.
            %                Also acceptable: 'all' to load all wells.
            %                If this parameter is omitted, all wells are loaded.
            %                For a single-well file, this parameter is ignored.
            %
            %    electrode   Which electrodes to load.  Format is either a comma-delimited string
            %                with whitespace ignored (e.g. '11, 22,33') or a single channel number;
            %                that is, a number, not part of a string.
            %                Also acceptable: 'all' to load all channels and 'none', '-1', or -1
            %                to load no data (returns only header information).
            %                If this parameter is omitted, all channels are loaded.
            %
            %    timespan    Span of time, in seconds, over which to load data.  Format is a two-element
            %                array, [t0 t1], where t0 is the start time and t1 is the end time and both
            %                are in seconds after the first sample in the file.  Samples returned are ones
            %                that were taken at time >= t0 and <= t1.  The beginning of the file
            %                is at 0 seconds.
            %                If this parameter is omitted, the data is not filtered based on time.
            %
            
            fLoadArgs = LoadArgs(varargin{:});
            fTargetWell = fLoadArgs.Well;
            fTargetElectrode = fLoadArgs.Electrode;
            fTimeRange = fLoadArgs.Timespan;
            
            aData = [];
            aData.header = [];
            
            %Check for Required ource file
            fSourceFile = this.sourceFile;
            if(~isa(fSourceFile,'AxisFile'))
                error('BlockVectorSet: Source file not specified');
            end
            
            %Check for optional BlockVectorHeaderExtension and Notes
            if(~isempty(fSourceFile.Notes))
                [~,idx]=max([fSourceFile.Notes.Revision]);

                fNotes = fSourceFile.Notes(idx);
                aData.notes = [];
                aData.notes.investigator = fNotes.Investigator;
                aData.notes.experimentId = fNotes.ExperimentID;
                aData.notes.notes = fNotes.Description;
            end
            
            aData.header.fileTypeNumber     = fSourceFile.PrimaryDataType; %Default incase we don't have a header extension
            aData.header.headerVersionMajor = fSourceFile.HeaderVersionMajor;
            aData.header.headerVersionMinor = fSourceFile.HeaderVersionMinor;

            %Check for Required Entries
            fChannelArray = this.ChannelArray;
            if(~isa(fChannelArray,'ChannelArray'))
                error('BlockVectorSet: Channel array was not specified for this file.');
            end
            
            fHeader = this.Header;
            if(~isa(fHeader,'BlockVectorHeader'))
                error('BlockVectorSet: BlockVectorHeader was not specified for this file.');
            end
            
            fData = this.Data;
            if(~isa(fData,'BlockVectorData'))
                error('BlockVectorSet: No Data found');
            end
            
            %Check for optional BlockVectorHeaderExtension
            fHeaderExtension = this.HeaderExtension;
            if(isa(fHeaderExtension,'BlockVectorHeaderExtension'))
                aData.header.fileTypeNumber = fHeaderExtension.DataType;
            end
            
            %Handle Channel Array Data
            aData.channelArray = [];
            aData.channelArray.plateType = fChannelArray.PlateType;
            aData.channelArray.numChannels = length(fChannelArray.Channels);

            aData.channelArray.channel = [];
            for i = 1:aData.channelArray.numChannels
                aData.channelArray.channel(i).wellColumn      = fChannelArray.Channels(i).WellColumn;
                aData.channelArray.channel(i).wellRow         = fChannelArray.Channels(i).WellRow;
                aData.channelArray.channel(i).electrodeColumn = fChannelArray.Channels(i).ElectrodeColumn;
                aData.channelArray.channel(i).electrodeRow    = fChannelArray.Channels(i).ElectrodeRow;
                aData.channelArray.channel(i).channelAchk     = fChannelArray.Channels(i).ChannelAchk;
                aData.channelArray.channel(i).channelIndex    = fChannelArray.Channels(i).ChannelIndex;
                aData.channelArray.channel(i).auxData         = fChannelArray.Channels(i).AuxData;
            end

            %Handle BlockVectorHeaders Data
            
            aData.samplingFrequency = fHeader.SamplingFrequency;
            aData.voltageScale = fHeader.VoltageScale;

            aData.header.numChannelsPerBlock = fHeader.NumChannelsPerBlock;
            aData.header.numSamplesPerBlock  = fHeader.NumSamplesPerBlock;
            aData.header.blockHeaderSize     = fHeader.BlockHeaderSize;

            if(~isempty(fHeader.FileStartTime))

                aData.fileStartTime = [];

                aData.fileStartTime.year        = fHeader.FileStartTime.Year;
                aData.fileStartTime.month       = fHeader.FileStartTime.Month;
                aData.fileStartTime.day         = fHeader.FileStartTime.Day;
                aData.fileStartTime.hour        = fHeader.FileStartTime.Hour;
                aData.fileStartTime.minute      = fHeader.FileStartTime.Minute;
                aData.fileStartTime.second      = fHeader.FileStartTime.Second;
                aData.fileStartTime.millisecond = fHeader.FileStartTime.Millisecond;

            end

            if(~isempty(fHeader.ExperimentStartTime))

                aData.experimentStartTime = [];
                aData.experimentStartTime.year        = fHeader.ExperimentStartTime.Year;
                aData.experimentStartTime.month       = fHeader.ExperimentStartTime.Month;
                aData.experimentStartTime.day         = fHeader.ExperimentStartTime.Day;
                aData.experimentStartTime.hour        = fHeader.ExperimentStartTime.Hour;
                aData.experimentStartTime.minute      = fHeader.ExperimentStartTime.Minute;
                aData.experimentStartTime.second      = fHeader.ExperimentStartTime.Second;
                aData.experimentStartTime.millisecond = fHeader.ExperimentStartTime.Millisecond;

            end

            fChannelsToLoad = BlockVectorSet.get_channels_to_load(fChannelArray, fTargetWell, fTargetElectrode);
            aData.loadedChannels = fChannelsToLoad;
            
            %Handle Data 
            
            if aData.header.fileTypeNumber == BlockVectorDataType.Raw_v1
                aData.fileType = 'raw';
                aData = BlockVectorLegacyLoader.Legacy_Load_Raw_v1(fData, aData, fTimeRange);

            elseif aData.header.fileTypeNumber == BlockVectorDataType.Spike_v1
                aData.fileType = 'spike';
                aData = BlockVectorLegacyLoader.Legacy_Load_spike_v1(fData, aData, fTimeRange);

            else
                error('BlockVectorSet: Unsupported file type number %u', aData.header.fileTypeNumber);
            end
            
            aData = rmfield(aData, 'header');
        end
    end
    
    methods(Static, Access = private)
        
        function ChannelListOut = get_channels_to_load(aChannelArray, aTargetWells, aTargetElectrodes)
            
            % Decode the aTargetWells string
            if strcmp(aTargetWells, 'all')
                % User has requested all wells - figure out what those
                % are from the channel array
                fTargetWells = BlockVectorSet.all_wells_electrodes([aChannelArray.Channels.WellColumn], ...
                    [aChannelArray.Channels.WellRow]);
            else
                fTargetWells = aTargetWells;
            end
            
            % Decode the aTargetElectrodes string
            if strcmp(aTargetElectrodes, 'all')
                % User has requested all electrodes - figure out what those
                % are from the channel array
                fTargetElectrodes = BlockVectorSet.all_wells_electrodes([aChannelArray.Channels.ElectrodeColumn], ...
                    [aChannelArray.Channels.ElectrodeRow]);
            elseif strcmp(aTargetElectrodes, 'none')
                % User has requested no electrodes
                fTargetElectrodes = [];
            else
                fTargetElectrodes = aTargetElectrodes;
            end
            
            ChannelListOut = zeros(1, size(fTargetWells, 1) * size(fTargetElectrodes, 1));
            if ~isempty(ChannelListOut)
                
                for fChannelArrayIndex = 1:length(aChannelArray.Channels)
                    fCurrentChannel = aChannelArray.Channels(fChannelArrayIndex);
                    
                    [fFoundWell, fIdxWell] = ismember( [fCurrentChannel.WellColumn  fCurrentChannel.WellRow], ...
                        fTargetWells, 'rows');
                    if ~any(fFoundWell)
                        continue;
                    end
                    
                    [fFoundElectrode, fIdxElectrode ] = ismember( [fCurrentChannel.ElectrodeColumn  fCurrentChannel.ElectrodeRow], ...
                        fTargetElectrodes, 'rows');
                    
                    if ~any(fFoundElectrode)
                        continue;
                    end
                    
                    ChannelListOut( (fIdxWell - 1) * size(fTargetElectrodes, 1) + fIdxElectrode ) = fChannelArrayIndex;
                end
                
                % Notify the user of any requested channels that weren't found in the channel array.
                % This is not necessarily an error; for example, if a whole well is requested, and
                % some channels in that well weren't recorded, we should return the well without
                % the "missing" channel.
                fChannelIdxZeros = find(ChannelListOut == 0);
                for i=1:length(fChannelIdxZeros)
                    fIdxNotFound = fChannelIdxZeros(i);
                    fMissingWell = floor((fIdxNotFound-1) / size(fTargetElectrodes, 1)) + 1;
                    fMissingElectrode = mod(fIdxNotFound-1, size(fTargetElectrodes, 1)) + 1;
                    warning('get_channels_to_load:invalidWellElectrode', ...
                        sprintf('Well/electrode %d %d / %d %d not recorded in file', ...
                        fTargetWells(fMissingWell, 1), fTargetWells(fMissingWell, 2), ...
                        fTargetElectrodes(fMissingElectrode, 1), fTargetElectrodes(fMissingElectrode, 2)));
                end
                
                % Strip out any zeros from aChannelListOut, because these correspond to channels that weren't in
                % the loaded channel array, and therefore won't be loaded.
                ChannelListOut = ChannelListOut( ChannelListOut ~= 0 );
            end
            
        end % end function
        
        % Subfunction to expand an 'all' well or electrode list
        function fOutput = all_wells_electrodes(aColumns, aRows)
            
            fOutput  = [];
            aColumns = unique(aColumns); % sort ascending and dedup
            aRows    = unique(aRows);
            
            for fiRow = 1:length(aRows)
                for fiCol = 1:length(aColumns)
                    fOutput = [ fOutput ; aColumns(fiCol) aRows(fiRow) ];
                end
            end
            
        end
        
        % Subfunction to help with channel array search
        function aMatch = match_well_electrode(aChannelStruct, aWellElectrode)
            
            if aChannelStruct.wellColumn == aWellElectrode(1) && ...
                    aChannelStruct.wellRow    == aWellElectrode(2) && ...
                    aChannelStruct.electrodeColumn == aWellElectrode(3) && ...
                    aChannelStruct.electrodeRow    == aWellElectrode(4)
                aMatch = 1;
            else
                aMatch = 0;
            end
            
        end
    end
    
end

