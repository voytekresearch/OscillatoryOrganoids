classdef AxisFile < handle
    %AXISFILE: Class that holds all Data for a loaded Axion 1.X File
    %
    % Calling A = AxisFile('filename.raw') or AxisFile('filename.spk') will
    % create an AxisFile object which contains the data set corresponding
    % to file 'y', as well as the notes and file name.
    %
    %   A = AxisFile('filename.raw') - loads file header information
    %
    % From here, the command 'z = A.DataSets(1)' stores a copy of the handle
    % to the first "BlockVectorSet" object in the file, which acts as a data
    % map for the user.
    %
    % y = z.LoadData() of z, loads specific data regions of the dataset
    % to y without reloading the entire file each time. For example:
    %
    %   a1data = z.LoadData('A1', 11)
    %
    % will load all A1_11 information in the file and store it in a1Data.
    %
    %   y = x.DataSets.LoadData(well, electrode, timespan, dimensions) - loads
    %      channel data into array
    %
    %      dimension 1 ->   by plate, vector with 1 waveform per signal in
    %                       the plate
    %      dimension 3 ->   by well, cell array of vectors 1 waveform per
    %                       signal with size (well rows, well columns)
    %      dimension 5 ->   by electrode, reference electrodes by y{well row,
    %                       well column, electrode column, electrode row}
    %
    % The LoadData() function loads data into instances of the class
    % 'Waveform', where the voltage and time vector data can be accessed via
    % the methods GetTimeVector() and GetVoltageVector(). See Also: Waveform.m
    %
    %   y{wr, wc, ec, er}.GetTimeVector - returns time vector based on start
    %                       time, length of data, and Sampling Frequency
    %                       in units of Seconds
    %   y{wr, wc, ec, er}.GetVoltageVector - returns the voltage vector
    %                       with units of Volts
    %
    % Properties of AxisFile include:
    %
    %   FileName:           Path to file that this instance points to
    %
    %   PrimaryDataType:    The type of the Original / most relevant data
    %                       in this file. See Also: BlockVectorDataType.m
    %
    %   HeaderVersionMajor: Major version number of the loaded file.
    %                       (modern files should be 1)
    %
    %   HeaderVersionMinor: Minor version number of the loaded file.
    %
    %   Notes:              User Notes on this file See Also: Note.m
    %
    %   DataSets            Objects with access to the data and meta-data
    %                       contained by this file See Also: BlockVectorSet.m
    %
    % See Also: BlockVectorSet, Waveform, BlockVectorSet,
    % BlockVectorDataType, Note
    
    
    
    properties (Constant = true, GetAccess = private)
        %The following are contants used in the opening of Axis Files
        MAGIC_WORD      = 'AxionBio';           % Preface to all modern Axis files
        MAGIC_BETA_FILE = 64;                   % Preface to Some legacy Axis files
        EXPECTED_NOTES_LENGTH_FIELD = 600;      % Number used as a validity check in Axis 1.0 file headers
        
        %Header CRC32 calculation constants
        CRC_POLYNOMIAL = hex2dec('edb88320');
        CRC_SEED = hex2dec('ffffffff');
        
        %Header Size Constants, see documentation for details
        PRIMARY_HEADER_CRCSIZE = 1018;
        SUBHEADER_CRCSIZE = 1016;
        PRIMARY_HEADER_MAXENTRIES = 123;
        SUBHEADER_MAXENTRIES = 126;
    end
    
    properties (Constant = true, GetAccess = public)
        % Version of AxIS this script is released with
        AXIS_VERSION='2.1.2.9';
    end
    
    properties (SetAccess = private, GetAccess = private)
        FileID;                % File Handle for file Access (fread, ftell, etc...)
        NotesStart;            % Location in file (in Number of bytes from the beginning) of the primary notes field
        EntriesStart;          % Staring byte of file entries
    end
    
    properties (SetAccess = private, GetAccess = public)
        %Basic File Data
        FileName;
        PrimaryDataType;
        
        %Version Number: Current version is (1.0)
        HeaderVersionMajor;
        HeaderVersionMinor;
        
        %Contained File Data
        Notes;
        DataSets;
        Annotations;
        PlateMap;
        StimulationEvents;
    end
    
    methods
        
        function this = AxisFile(varargin)
            %AxisFile Opens a new handle to an Axis File
            %  Required arguments:
            %    filename    Pathname of the file to load
            %                Note: Calling with mulitple file names results
            %                in vector of AxisFile objects corresponding to
            %                the input argument file names
            
            
            if(nargin == 0)
                this.FileID = [];
                return;
            end
            if(nargin > 1)
                this(nargin) = AxisFile;
                for i = 1: nargin
                    this(i) = AxisFile(varargin{i});
                end
                return;
            end
            fFilename = varargin{1};            
            this.FileName = fFilename;
            this.FileID = fopen(fFilename,'r');
            
            fSetMap = containers.Map('KeyType', 'int64', 'ValueType', 'any');
            this.Notes = Note.empty(0,0);
            
            if (this.FileID <= 0)
                error(['AxisFile: ' this.FileName ' not found.']);
            end
            % Make sure that this is a format that we understand
            versionOk = false;
            versionWarn = false;
            
            % Check for the "magic word" sequence
            fMagicRead = fread(this.FileID, length(AxisFile.MAGIC_WORD), '*char').';
            if ~strcmp(AxisFile.MAGIC_WORD, fMagicRead)
                
                % Magic phrase not found -- check to see if this is an old-style file
                if ~isempty(fMagicRead) && uint8(fMagicRead(1)) == AxisFile.MAGIC_BETA_FILE
                    % This looks like a deprecated beta file
                    warning('AxisFile:versionCheck', ['File ' fFilename ' looks like a deprecated AxIS v0.0 file format, Please Re-record it in Axis to update the header data']);
                    
                    [fType, fData, fChannelMapping, fHeader] = LegacySupport.GenerateRolstonEntries(this.FileID, fFilename((end-3):end));
                    
                    this.HeaderVersionMajor = 0;
                    this.HeaderVersionMinor = 0;
                    this.PrimaryDataType = fType;
                    this.DataSets = BlockVectorSet(this, fData, fChannelMapping, fHeader);
                    
                    return;
                else
                    fclose(this.FileID);
                    error('File format not recognized: %s', fFilename);
                end
                
            else
                
                this.PrimaryDataType         = fread(this.FileID, 1, 'uint16=>uint16');
                this.HeaderVersionMajor      = fread(this.FileID, 1, 'uint16=>uint16');
                this.HeaderVersionMinor      = fread(this.FileID, 1, 'uint16=>uint16');
                this.NotesStart              = fread(this.FileID, 1, 'uint64=>uint64');
                fNotesLength                 = fread(this.FileID, 1, 'uint32=>uint32');
                
                if(fNotesLength ~= AxisFile.EXPECTED_NOTES_LENGTH_FIELD)
                    error('Incorrect legacy notes length field');
                end
                
                if this.HeaderVersionMajor == 0
                    if this.HeaderVersionMinor == 1
                        versionOk = true;
                    elseif this.HeaderVersionMinor == 2
                        versionOk = true;
                    end
                    
                    this.EntriesStart = int64(this.NotesStart);
                    fEntryRecords = LegacySupport.GenerateEntries(this.FileID, this.EntriesStart);
                    
                elseif this.HeaderVersionMajor == 1
                    versionOk = true;
                    
                    this.EntriesStart        = fread(this.FileID, 1, 'int64=>int64');
                    fEntrySlots = fread(this.FileID, AxisFile.PRIMARY_HEADER_MAXENTRIES, 'uint64=>uint64');
                    fEntryRecords = EntryRecord.FromUint64(fEntrySlots);
                    
                    % Check CRC
                    fseek(this.FileID, 0, 'bof');
                    fCRCBytes = fread(this.FileID, AxisFile.PRIMARY_HEADER_CRCSIZE, 'uint8');
                    fReadCRC = fread(this.FileID, 1, 'uint32');
                    fCalcCRC = CRC32(AxisFile.CRC_POLYNOMIAL, AxisFile.CRC_SEED).Compute(fCRCBytes);
                    
                    if(fReadCRC ~= fCalcCRC)
                        error('File header checksum was incorrect: %s', fFilename);
                    end
                    
                    if this.HeaderVersionMinor > 0
                        versionWarn = true;
                    end
                end
                
            end
            
            if ~versionOk
                error('Unsupported file version %u.%u', ...
                    this.HeaderVersionMajor, ...
                    this.HeaderVersionMinor);
            end
            
            % Start Reading Entries
            fseek(this.FileID, this.EntriesStart, 'bof');
            
            fTerminated = false;
            
            fTagEntries = TagEntry.empty(0);
            
            % Load file entries from the header
            while(~fTerminated)
                
                for entryRecord = fEntryRecords
                    switch(entryRecord.Type)
                        
                        case EntryRecordID.Terminate
                            fTerminated = true;
                            break
                            
                        case EntryRecordID.ChannelArray
                            fChannelArray = ChannelArray(entryRecord, this.FileID);
                            if(~isa(fCurrentBlockVectorSet.ChannelArray , 'ChannelArray'))
                                fCurrentBlockVectorSet = fCurrentBlockVectorSet.Clone(fChannelArray);
                                fSetMap(int64(fCurrentHeader.FirstBlock)) = fCurrentBlockVectorSet;
                            else
                                error('AxisFile: Only one ChannelArray per BlockVectorSet');
                            end
                            
                        case EntryRecordID.BlockVectorHeader
                            fCurrentHeader = BlockVectorHeader(entryRecord, this.FileID);
                            
                            fCurrentBlockVectorSet = BlockVectorSet(this, fCurrentHeader);
                            
                            fKey = int64(fCurrentHeader.FirstBlock);
                            fSetMap(fKey) = fCurrentBlockVectorSet;
                            
                        case EntryRecordID.BlockVectorHeaderExtension
                            if(~isempty(fCurrentBlockVectorSet.HeaderExtension) || ...
                                    isa(fCurrentBlockVectorSet.HeaderExtension, 'BlockVectorHeaderExtension'))
                                error('AxisFile: Only one BlockVectorHeaderExtension per BlockVectorSet');
                            end
                            fCurrentBlockVectorSet = fCurrentBlockVectorSet.Clone(BlockVectorHeaderExtension(entryRecord, this.FileID));
                            fSetMap(fCurrentBlockVectorSet.Header.FirstBlock) = fCurrentBlockVectorSet;
                            
                        case EntryRecordID.BlockVectorData
                            fData = BlockVectorData(entryRecord, this.FileID);
                            if(~isempty(fCurrentBlockVectorSet.Data) || ...
                                    isa(fCurrentBlockVectorSet.Data, 'BlockVectorData'))
                                error('AxisFile: Only one BlockVectorData per BlockVectorSet');
                            end
                            fTargetSet = fSetMap(int64(fData.Start));
                            if (~isa(fTargetSet, 'BlockVectorSet'))
                                error('AxisFile: No header to match to data');
                            end
                            fTargetSet = fTargetSet.Clone(fData);
                            fSetMap(fData.Start) = fTargetSet;
                            
                        case EntryRecordID.NotesArray
                            this.Notes = [this.Notes ; Note.ParseArray(entryRecord, this.FileID)];
                            
                        case EntryRecordID.Tag
                            fTagEntries(end+1) = TagEntry(entryRecord, this.FileID);
                            
                        otherwise
                            fSkipSpace = double(entryRecord.Length);
                            if(0 ~= fseek(this.FileID, fSkipSpace, 'cof'))
                                error(ferror(this.FileID));
                            end
                            
                    end
                end
                
                if(~fTerminated)
                    
                    %Check Magic Bytes
                    fMagicRead = fread(this.FileID, length(AxisFile.MAGIC_WORD), '*char').';
                    if ~strcmp(AxisFile.MAGIC_WORD, fMagicRead)
                        error('Bad sub header magic numbers: %s', fFilename);
                    end
                    
                    %Read Entry Records
                    fEntrySlots = fread(this.FileID, AxisFile.SUBHEADER_MAXENTRIES, 'uint64=>uint64');
                    fEntryRecords = EntryRecord.FromUint64(fEntrySlots);
                    
                    %Check CRC of subheader
                    fseek(this.FileID,( -1 * length(AxisFile.MAGIC_WORD)) - (8 * AxisFile.SUBHEADER_MAXENTRIES),'cof');
                    fCRCBytes = fread(this.FileID, AxisFile.SUBHEADER_CRCSIZE, 'uint8');
                    fReadCRC = fread(this.FileID, 1, 'uint32');
                    fCalcCRC = CRC32(AxisFile.CRC_POLYNOMIAL, AxisFile.CRC_SEED).Compute(fCRCBytes);
                    if(fReadCRC ~= fCalcCRC)
                        error('Bad sub header checksum : %s', fFilename);
                    end
                    
                    %skip 4 reserved bytes
                    fseek(this.FileID, 4,'cof');
                end
                
            end
            
            fValueSet = fSetMap.values;
            
            %Record Final Data Sets
            this.DataSets = BlockVectorSet.empty(0,length(fSetMap));
            for i = 1 : length(fValueSet)
                this.DataSets(i) = fValueSet{i};
            end
            
            %Sort Notes
            [~,idx]=sort([this.Notes.Revision]);
            this.Notes = this.Notes(idx);
            
            %Collect Tags
            fTagMap = containers.Map();
            for fEntryNum = 1:length(fTagEntries)
                fEntry = fTagEntries(fEntryNum);
                fGuid = fEntry.TagGuid;
                if fTagMap.isKey(fGuid)
                    fTag = fTagMap(fGuid);
                else
                    fTag = Tag(fGuid);
                    fTagMap(fGuid) = fTag;
                end
                fTag.AddNode(fEntry);
            end
            this.Annotations = Annotation.empty(0);
            this.PlateMap = WellInformation.empty(0);
            this.StimulationEvents = StimulationEvent.empty(0);
            
            for fKey = fTagMap.keys;
                ffKey = fKey{1};
                fTag = fTagMap(ffKey).Promote(this.FileID);
                if isa(fTag, 'Annotation')
                    this.Annotations(end+1) = fTag;
                elseif isa(fTag, 'WellInformation')
                    this.PlateMap(end+1) = fTag;                    
                elseif isa(fTag, 'StimulationEvent')
                    this.StimulationEvents(end+1) = fTag;
                end
                fTagMap(ffKey) = fTag;
            end
            
           
            for fStimEvent = this.StimulationEvents
                fStimEvent.Link(fTagMap);
            end
            
            this.Annotations = this.Annotations';
            this.PlateMap = this.PlateMap';
            this.StimulationEvents = this.StimulationEvents';
                        
        end
        
        
        function delete(this)
            %DELETE is the destructor for the class, ensures that the file
            %stream is closed as the file reference is cleared from the
            %workspace
            
            if ~isempty(this.FileID)
                fclose(this.FileID);
            end
        end
        
    end
    
end



