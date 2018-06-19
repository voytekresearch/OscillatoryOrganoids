classdef Tag < handle
    %TAG Base class of user genreated metadata for a file
    %   Currently, AxionFileTags include Annotation (from the play bar), 
    %   StimulationEvent (when stimulations occur), and WellInformation
    %   (from the Plate Map Editor).
    
    properties(GetAccess = private, SetAccess = private)
        EntryNodes;
    end
    
    properties(GetAccess = public, SetAccess = private)
        % TagGuid: Unique GUID for this tag and its revision history. 
        % When a new tag is added to a TagCollection, if its GUID matches 
        % that of an existing tag, that tag is updated with the new tag as a revision.
		TagGuid;
        
        % HeadRevisionNumber: Each Revision of a tag by axis has a new, 
        % iterated revison number. This is the most recent tag's revision
        % number
        HeadRevisionNumber;    
        
        %Type: The TypeID attributed to this tag
        Type;
    end
    
    methods
        function this = Tag(aGuid)
           this.TagGuid = aGuid;
           this.HeadRevisionNumber = -1;
           this.EntryNodes = TagEntry.empty(0,1);
        end
        
        function new = Promote(this, aFileId)
            %%% Promote:
            % Converts a base tag to the type that is dictated by its Type
            % property. Note that the returned object is a new instance
           [~,idx] = sort(arrayfun(@(a)(a.RevisionNumber), this.EntryNodes));
           fEntryNodes = this.EntryNodes(idx);
           fHead = fEntryNodes(end);
           switch fHead.Type
               case TagType.UserAnnotation
                  new = Annotation(aFileId, fHead);
               case TagType.SystemAnnotation    
                  new = Annotation(aFileId, fHead);
               case TagType.WellTreatment
                  new = WellInformation(aFileId, fHead);
               case TagType.StimulationEvent
                  new = StimulationEvent(aFileId, fHead);
               case TagType.StimulationChannelGroup
                  new = StimulationChannels(aFileId, fHead);
               case TagType.StimulationWaveform
                  new = StimulationWaveform(aFileId, fHead);
               case TagType.CalibrationTag
                  %For Calibration Tags, No Additonal Parsing supported
                  new = Tag(this.TagGuid);                   
               otherwise
                  new = Tag(this.TagGuid);
                  if(this.Type ~= TagType.Deleted)
                     warning('Unknown Tag Type found. Is this loader out of date?');
                  end
           end
           new.EntryNodes = fEntryNodes;
           new.HeadRevisionNumber = fHead.RevisionNumber;
           new.Type = fHead.Type;
        end
        
        function AddNode(this, aNode)
           %%% AddNode
           % Adds a TagEntry to the revision history of this Tag series
           this.EntryNodes = [this.EntryNodes; aNode];
           %sort by revison number
           [this.HeadRevisionNumber, idx] = max(arrayfun(@(a)(a.RevisionNumber), this.EntryNodes));
           this.Type = this.EntryNodes(idx).Type;
        end
    end
    
end

