classdef Waveform
    %WAVEFORM Container for single dimensional recorded sample data
    %
    %   Channel:    Source Location of the waveform
    %
    %   Start:      Time(In Seconds) of Recording Start
    %
    %   Data:       Sample data
    %
    %   Source:     BlockVectorDataSet that contains this Waveform
    
    properties(GetAccess = public, SetAccess = protected)
        Channel;
        Start;
        Data;
        Source;
    end
    
    methods
        
        function this = Waveform(aChannel, aStart, aData, aSource)
            
            if(nargin == 0)
                return;
            end
            
            if(~isa(aChannel,'ChannelMapping'))
                error(['Waveform: Unexpected Argument for aChannel: ' aChannel]);
            end
            
            if(~isa(aSource,'BlockVectorSet'))
                error(['Waveform: Unexpected Argument for aSource: ' aSource]);
            end
            
            this.Channel = aChannel;
            this.Start = aStart;
            this.Data = aData;
            this.Source = aSource;
        end
        
        function [timeData, voltageData] = GetTimeVoltageVector(this)
            %GetTimeVoltageVector: Returns a vector for time and voltage
            % for this waveform in a single call
            timeData = this.GetTimeVector();
            voltageData = this.GetVoltageVector();
        end
        
        function voltageData = GetVoltageVector(this)
            % GetTimeVector: returns a voltage vector for this waveform based
            % on the uncasted sample data (Stored as int16) and the source
            % header's specified voltage scale 
            %
            % If this Method is called on an array of waveforms, the
            % lengths of the waveforms MUST agree
            fData = double([this(:).Data]);
            fSource = [this(:).Source];
            fHeader = [fSource(:).Header];
            fVoltageScale = [fHeader(:).VoltageScale];
            voltageData = fData * diag(fVoltageScale);
        end
        
        function timeData = GetTimeVector(this)
            % GetTimeVector: returns a time vector for this waveform based
            % on the Start time, Length of the data, and the Sampling
            % Frequency of the source header
            %
            % If this Method is called on an array of waveforms, the
            % lengths of the waveforms MUST agree
            fSource = [this(:).Source];
            fHeader = [fSource(:).Header];
            fSamplingPeriod = 1./[fHeader(:).SamplingFrequency];
            
            timeData = repmat((0 : (length(this(1).Data) - 1))', 1,length(this));
            timeData = timeData * diag(fSamplingPeriod);
            
            fStart = ones(size(timeData));
            fStart = fStart * diag([this(:).Start]);
            
            timeData = timeData + fStart;
        end
        
    end
    
end

