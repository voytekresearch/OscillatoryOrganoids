classdef PlateTypes
    %PlateTypes finds dimentions for known plate types
    %   This class has to be modified every time a new plate type is added 
    %   to AxIS
    
    properties(Constant, Access=public)
        % All the channels of a Muse, in artichoke order 
        LinearSingleWell64 = uint32(hex2dec('0400000')); 
        
        % Muse single-well plate
        P200D30S = uint32(hex2dec('0400001')); 
        
        % All the channels of a Maestro, in artichoke order 
        LinearTwelveWell = uint32(hex2dec('3000000'));

		% Standard Maestro 12 well plate
		TwelveWell = uint32(hex2dec('3000001'));

		% Opaque Maestro 48 well plate
		FortyEightWell = uint32(hex2dec('3000002'));

		% Standard Maestro 96 well plate
		NinetySixWell =	uint32(hex2dec('3000003'));

		% Transparent Maestro 48 well plate
		FortyEightWellTransparent =	uint32(hex2dec('3000004'));
		
		% Standard Maestro 384 well plate
		ThreeEightyFourWell =	uint32(hex2dec('3000005'));
    end
    
    properties (Constant, Access=private)
        
        MUSE_MASK = uint32(hex2dec('0400000'));
        MAESTRO_MASK = uint32(hex2dec('3000000'));
        
        MuseElectrodeMap = [1, 1, 8, 8; ...  % LinearSingleWell64
                            1, 1, 8, 8];     % P200D30S
                        
        MaestroElectrodeMap = [ 3, 4, 8, 8; ...  %LinearTwelveWell
                                3, 4, 8, 8; ...  %TwelveWell
                                6, 8, 4, 4; ...  %FortyEightWell
                                8, 12, 3, 3; ... %NinetySixWell
                                6, 8, 4, 4;  ... %FortyEightWellTransparent
								16, 24, 2, 1];...%ThreeEightyFourWell%
    end
    
    methods(Static)
        function fPlateDimentions = GetWellDimensions(aPlateType)
            % GetWellDimensions returns a 2-element array of plate
            % dimensions.
            %
            % First element is the number of well rows, second element
            % is the number of well columns.
            %
            
            offset = bitand(aPlateType, 15);
            
            if (bitand(aPlateType, PlateTypes.MUSE_MASK) == PlateTypes.MUSE_MASK)
               fPlateDimentions = PlateTypes.MuseElectrodeMap(offset + 1, (1:2));
            elseif (bitand(aPlateType, PlateTypes.MAESTRO_MASK) == PlateTypes.MAESTRO_MASK)
               fPlateDimentions = PlateTypes.MaestroElectrodeMap(offset + 1, (1:2));
            else
                warning('File has and unknown plate type. These Matlab Scripts may be out of date.');
                fPlateDimentions = [];
            end
        end
        
        function fElectrodeDimentions = GetElectrodeDimensions(aPlateType)
            % GetElectrodeDimensions returns a 4-element array of plate
            % dimensions (wells and electrodes within wells).
            %
            % Format is [well rows, well columns, electrode rows, electrode
            % columns].
            %
            % NOTE:  wells of a 96-well plates have 3 electrode rows an 3
            % electrode columns.  However, the second row contains only 2 
            % valid electrodes.
            %
            
            offset = bitand(aPlateType, 15);
            
            if (bitand(aPlateType, PlateTypes.MUSE_MASK) == PlateTypes.MUSE_MASK)
               fElectrodeDimentions = PlateTypes.MuseElectrodeMap(offset + 1, :);
            elseif (bitand(aPlateType, PlateTypes.MAESTRO_MASK) == PlateTypes.MAESTRO_MASK)
               fElectrodeDimentions = PlateTypes.MaestroElectrodeMap(offset + 1, :);
            else
                warning('File has and unknown plate type. These Matlab Scripts may be out of date.');
                fElectrodeDimentions = [];
            end
        end
    end
    
end

