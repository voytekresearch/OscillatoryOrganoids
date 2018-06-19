function plot_tight(data_array, dim, x_ax, y_ax)
%plot_tight(data_array, dim, x_ax, y_ax)
figure
nrow = dim(1);
ncol = dim(2);

h = 1/nrow - 0.01;
w = 1/ncol - 0.01;
if isempty(x_ax)
    x_ax = 1:size(data_array,1);
end


for row = 1:nrow
    for col = 1:ncol
        %2D dataset
        if length(size(data_array))==2
            if any(data_array(:, col+(row-1)*ncol))
                subplot('Position', [(col-1)/ncol+0.005, 1-(row)/nrow+0.005 w h]);
                plot(x_ax,data_array(:, col+(row-1)*ncol), 'k');
                %loglog(x_ax,data_array(:, col+(row-1)*ncol), 'k');
                xlim([x_ax(1),x_ax(end)])
                if ~isempty(y_ax)
                    ylim([y_ax(1),y_ax(end)])
                end
            end
            %Hleg=legend(num2str(col+(row-1)*ncol));
            %set(Hleg,'fontsize',10)
        end
        
        %3D dataset
        if length(size(data_array))==3
            if any(any(data_array(:,:, col+(row-1)*ncol)))
                subplot('Position', [(col-1)/ncol+0.005, 1-(row)/nrow+0.005 w h]);
                imagesc(x_ax, y_ax, data_array(:,:,col+(row-1)*ncol)');
            end
        end
        set(gca,'XTick',[]);set(gca,'YTick',[]); set(gca, 'ydir', 'normal');        
    end
end
% 
% 
% d1=dim(1);
% d2=dim(2);
% for i=1:d1
%     for j=1:d2
%         
%         
%         %2D dataset
%         if length(size(data_array))==2
%             if any(data_array(:,j+(i-1)*d2))
%                 subplot('Position', [mod(i-1,d1)/d1+0.005, 1-(mod(j-1,d2)/d2)-1/d2+0.005,  1/d1-0.008, 1/d2-0.008 ]);
%                 if isempty(x_ax)
%                     plot(data_array(:,j+(i-1)*d2));
%                 else
%                     plot(x_ax, data_array(:,j+(i-1)*d2));
%                     %bar(x_ax, data_array(:,j+(i-1)*d2));
%                 end
%                 %xlim([x_ax(1),x_ax(end)])
%                 %ylim([min(data_array(:,j+(i-1)*d2)), max(data_array(:,i+(j-1)*d2))]);
%                 %ylim([0 1])
%             end
%         end
%         
%         %3D dataset
%         if length(size(data_array))==3
%             if any(data_array(:,:,i+(j-1)*d2))
%                 subplot('Position', [mod(i-1,d1)/d1+0.005, 1-(mod(j-1,d2)/d2)-1/d2+0.005,  1/d1-0.008, 1/d2-0.008 ]);
%                 
%                 imagesc(x_ax, y_ax, data_array(:,:,i+(j-1)*d2)); %transpose because MI somehow got stored backwards, should fix
%                 
%                 %for SFC (wrapped phase)
%                 %imagesc(x_ax, y_ax, [data_array(:,:,i+(j-1)*d2) data_array(:,:,i+(j-1)*d2)]); %transpose because MI somehow got stored backwards, should fix
%     
%                 caxis([0 0.3])
%             end
%         end
%         set(gca,'XTick',[]);set(gca,'YTick',[]); set(gca, 'ydir', 'normal');        
%     end
% end
% end