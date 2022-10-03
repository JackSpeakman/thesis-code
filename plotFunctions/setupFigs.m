function [ax] = setupFigs(fig,forPub,useNiceFigPos,layout)

%% get positions
if useNiceFigPos || forPub
    figPos = getNiceFigPos(layout,forPub);
else
    % don't change from default/current layout
    figPos = zeros(5,4);
    for i = 1:5
        figPos(i,:) = get(fig(i),'Position');
    end
end

%% set positions and make axes
for i = 1:numel(fig)
    
    % change figure
    set(0,'CurrentFigure',fig(i))  
    
    % create axis
    ax(i) = axes(fig(i));
    
    % turn hold on
    hold(ax(i),'on')
    
    % set x-axis size to number of iterations
    xlabel(ax(i),'Iteration, k','Interpreter','latex')
    
    % run fixAxis to set-up consistent axis
    if forPub
        fixAxis(fig(i),ax(i),'linewidth',2.5,'fontsize',20)
    else
        fixAxis(fig(i),ax(i),'linewidth',1,'fontsize',12)
    end
    
    % pull axis to top
    set(ax(i),'Layer','Top')

    % set position
    set(fig(i),'Position',figPos(i,:))
end

end