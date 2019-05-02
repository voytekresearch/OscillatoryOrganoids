function MEA_convert(raw_file, output_folder, wells2process)
%   MEA_convert(raw_file, output_folder)
%       raw_file: file path of .raw file to be unpacked
%       output_folder: file path of output folder to be created
%       wells2process: subset of wells to unpack
warning('off')
% set up intake of raw data
if ~exist(output_folder)
    % if the folder doesn't exist, create and do conversion
    mkdir(output_folder);
end
% load the axion file into matlab
disp('Loading .raw file...')
file = AxisFile(raw_file);

% compute plate dimension from file header
num_wells = length(unique([file.DataSets.ChannelArray.Channels.WellRow])) * ...
        length(unique([file.DataSets.ChannelArray.Channels.WellColumn]));    

if num_wells~=12 && num_wells~=48
    % if you get something weird, default to 12 wells
    num_wells=12;
end
%parameters & loading data into MATLAB
if num_wells == 12
    %dimensions of the plate (# wells per row,col)
    plate_dim=[3,4];
    %number of channels per well
    num_chan=64;
    row = {'A','B','C'};
    col = {'1' '2' '3' '4'};
elseif num_wells==48
    %dimensions of the plate (# wells per row,col)
    plate_dim=[6,8]; %[6,8] or [3,4]
    %number of channels per well
    num_chan=16; %16 or 64
    row = {'A','B','C' 'D','E','F'};
    col = {'1' '2' '3' '4' '5' '6' '7' '8'};
end

%gather wells & save into matfiles
for k=1:plate_dim(1)
    for l=1:plate_dim(2)
        
        % extract wells
        disp(sprintf('Well: (%i,%i)', k,l));        
        outname = sprintf([output_folder '/well_%i.mat'], (k-1)*plate_dim(2)+l);
        
        if any(wells2process==((k-1)*plate_dim(2)+l))                        
            
            % if we want to convert this well
            if ~exist(outname, 'file')
                
                %dataset too large, load well by well
                %load well
                data = file.DataSets.LoadData([row{k} col{l}]);                                
                if all(all(squeeze([cellfun(@isempty,data(k,l,:,:))])))
                    %there is nothing in the well, skip to next iteration
                    %of fore loop
                    disp(sprintf('Whole well (%i, %i) empty, skipped.',k,l));
                    continue
                    
                elseif any(any(squeeze([cellfun(@isempty,data(k,l,:,:))])))
                    %some channels are empty, have to go through each channel
                    %fill empty channels with 0s
                    disp('Empty channels encountered. Sorting...')
                    
                    %find first non-empty channel
                    [r,c] = find(1-squeeze(cellfun(@isempty,data(k,l,:,:))),1);
                    
                    % get time vector and data length from here
                    t = data{k,l,r,c}.GetTimeVector;
                    data_len = length(t);
                    
                    % pre initialize MEA array
                    MEA = zeros(data_len, num_chan);
                    chan_dims = size(squeeze(data(k,l,:,:)));                    
                    for rr = 1:chan_dims(1)
                        for cc = 1:chan_dims(2)
                            if ~isempty(data{k,l,rr,cc})
                                %get data of non-empty channel
                                MEA(:,(cc-1)*chan_dims(1)+rr) = data{k,l,rr,cc}.GetVoltageVector;
                            else
                                disp(sprintf('(%i, %i) empty.',rr,cc));
                            end
                        end
                    end                                                   
                else
                    %all channels have data, just convert and save
                    temp = [data{k,l,:,:}];
                    MEA = temp.GetVoltageVector;
                    t = temp(1).GetTimeVector;
                end                
                save(sprintf([output_folder '/well_%i.mat'], (k-1)*plate_dim(2)+l), 't', 'MEA', '-v7.3')
            else
                disp('Well conversion exists.')
            end
        else
            disp('Well skipped by user.')
        end                
    end
end