function prettyPlot(fig, ax, figW, figH, hasLabel)
box on
grid on
fig.Units = 'centimeters';
ax.Units = 'centimeters';
fig.Color = [1.0, 1.0, 1.0];
% fig.Color = 'none';
% ax.Color = 'none';
fig.Position =  [3, 3, figW+2.5, figH+2.5];
ax.Position = [2, 1.5, figW, figH];

ax.LineWidth = 1;
ax.FontName = 'Times New Roman';
ax.FontSize = 12;
ax.TickLabelInterpreter = 'latex';
ax.Title.Interpreter = 'latex';

if hasLabel > 0
    ax.XLabel.Units = 'normalized';
    ax.XLabel.Interpreter = 'latex';

    ax.YLabel.Units = 'normalized';
    ax.YLabel.Interpreter = 'latex';
    if hasLabel == 3
        ax.ZLabel.Units = 'normalized';
        ax.ZLabel.Interpreter = 'latex';
    end
end
if ~isempty(ax.Legend)
    ax.Legend.Interpreter = 'latex';
end
% fig.MenuBar = 'none';
% fig.ToolBar = 'none';
end