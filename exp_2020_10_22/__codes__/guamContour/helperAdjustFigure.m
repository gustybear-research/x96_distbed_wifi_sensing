function [] = helperAdjustFigure(figure_in, figure_out_name)
  % ax = gca;
  % fig = gcf;

  % set(gcf,'color',[0.5 0.5 0.5]);
  %pbaspect([1.5 1 1])
  
  %Added
  fileLocation = '';
  
  if isgraphics(figure_in)  %If is a handle to a fig, don't need to load
      fig = figure_in;
      ax = fig.CurrentAxes;
  else
    fig = openfig(figure_in);
    ax = gca; %Need to get axis of current figure
  end
  
  %set(ax,'color',[0.94117647058 0.94117647058 0.94117647058]);
  set(ax, 'FontSize', 25);
  
  set(ax, 'FontName', 'Helvetica');
  set(fig,'color',[0.94117647058 0.94117647058 0.94117647058]);
  set(fig, 'Position', get(0, 'Screensize'));
  fig.PaperPosition = [0 0 10 10];
  fig.InvertHardcopy = 'off';
  set(ax,'Position',[0.175 0.12 0.8 0.8]);
  
  set(ax,'DataAspectRatio',[1 1 1]); 
  %Probably don't need this, however
  %this can be used to change the x/y ratios
  
  %figure_out_name = sprintf('%s',figure_out_name);
  
  saveas(gcf,figure_out_name,'fig')
  %saveas(gcf,figure_out_name,'eps')
  saveas(gcf,figure_out_name,'png')  %Saving a png straight from here
                                    %keeps font size, saving a png 
                                    %from fig doesn't
  %pbaspect([1 1.5 1])
  %set(gca,'Position',[0 0 0.95 0.95])

  %fig.PaperPositionMode = 'auto';
  %ax.Position = [0.1 0.1 0.7 0.7];
  %ax.OuterPosition =[0 0 0.8 0.8];
  %axes('OuterPosition',[0.13 0.58 0.77 0.34]);
end