function nice_figure(f_handle, filename, dim)
%nice_figure(f_handle, filename, dim)
% plots a nice figure at filename (full file name) with the specified
% dimensions (dim) in inches (I think)

set(f_handle,'Units','Inches');
pos = get(f_handle,'Position');
set(f_handle,'PaperPosition',[0.05 0.05 dim],'PaperUnits','Inches','PaperSize',dim)
%f_handle.PaperPosition = [0 0 dim];
print('-dpdf','-r300', filename)