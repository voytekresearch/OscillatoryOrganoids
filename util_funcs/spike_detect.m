function [spikes spike_shapes] = spike_detect(data,fs,auto_detect)
%function [spikes spike_shapes] = spike_detect(data,auto_detect)
% data: time X chan (filtered)
% auto_detect: whether to use auto threshold (5.5std) or manual
mpd = 0.002*fs; %min peak distance 2ms
winlen = 0.004*fs; %spike shape, 8ms spike-centered window

[len, numchan] = size(data);
spikes = cell(numchan,3);
spike_shapes = cell(numchan,3);
for chan = 1:numchan    
    thr1 = median(abs(data(:,chan))/0.6745)*5; %robust estimation by Q.Quiroga
    thr2 = median(abs(data(:,chan))/0.6745)*5.5;   
    
    spikes{chan,3}=zeros(1,2);
    spike_shapes{chan,3} = zeros(1,2);
    
    auto_thr = thr2; %5.5 std
    if auto_detect
        %auto-detect at 5.5 std        
        if max(data(:,chan))>auto_thr 
            %skip if channel does not go above auto_thr            
            [p spikes{chan,1}] = findpeaks(data(:,chan),'minpeakheight',auto_thr, 'minpeakdistance',mpd); %peak spikes
            spike_shapes{chan,1} = collect_spikes(data(:,chan),[],spikes{chan,1},winlen); %peak spike shape
        end
        if min(data(:,chan))<-auto_thr
            [p spikes{chan,2}] = findpeaks(-data(:,chan),'minpeakheight',auto_thr, 'minpeakdistance',mpd); %trough spikes
            spike_shapes{chan,2} = collect_spikes(data(:,chan),[],spikes{chan,2},winlen); %trough spike shape
        end
        spikes{chan,3} = [length(spikes{chan,1}) length(spikes{chan,2})]; %num of spikes                                   
        spike_shapes{chan,3} = spikes{chan,3}; %num spikes        
    else
        %manual-detect, thresholding with pre-cutoff at 5.5 std
        if max(data(:,chan))>thr2 | min(data(:,chan))<-thr2
            %only display if signal crosses threshold
            figure
            plot(data(:,chan));
            line([ones(1,4); len*ones(1,4)], [thr1 -thr1 thr2 -thr2; thr1 -thr1 thr2 -thr2]);
            title(chan)
            set(gcf, 'Position', get(0,'Screensize')); % Maximize figure
            [x, y] = ginput(2);
            if isempty(y)
                %no selected spikes
                spikes{chan,1}=[];
                spikes{chan,2}=[];
                spikes{chan,3}=zeros(1,2);
            else
                %first thr is peak, second is trough
                if y(2)>y(1); y = flipud(y); end;
                [p spikes{chan,1}] = findpeaks(data(:,chan),'minpeakheight',y(1), 'minpeakdistance',mpd);
                [p spikes{chan,2}] = findpeaks(-data(:,chan),'minpeakheight',-y(2), 'minpeakdistance',mpd);
                spikes{chan,3} = [length(spikes{chan,1}) length(spikes{chan,2})];
                
                spike_shapes{chan,1} = collect_spikes(data(:,chan),[],spikes{chan,1},winlen);
                spike_shapes{chan,2} = collect_spikes(data(:,chan),[],spikes{chan,2},winlen);
                for i=1:2
                    if spikes{chan,3}(i)>0
                        subplot(1,2,i)
                        plot(mean(spike_shapes{chan,i},2), 'linewidth',2)
                        hold on
                        plot(mean(spike_shapes{chan,i},2)+std(spike_shapes{chan,i},0,2), 'color', 'r')
                        hold off
                        title(spikes{chan,3}(i))
                    end
                end
                pause
            end
            spike_shapes{chan,3} = spikes{chan,3};
            close all
        else
            disp(sprintf('Skipped Channel: %i',chan));
        end
    end
end