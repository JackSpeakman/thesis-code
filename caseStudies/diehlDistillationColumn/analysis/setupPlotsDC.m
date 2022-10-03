function setupPlotsDC(ax,m1,m2)
% setupPlotsDC sets up the axes for the DC plots with the correct axis labels, lims, and
% legend details.
% 
% ------ INPUT VARIABLES ------
% ax        1-by-3          cell of axes for the figures
% m1                        cell of line specs (line 1)
% m2                        cell of line specs (line 2)
%
% ------ EXAMPLES ------
% See '../../../methods/distillationColumn/modifierAdaptationmakePlotsMA_DCshort.m'

% r plot
xlabel(ax(1),'Iteration, $k$','Interpreter','latex');
ylabel(ax(1),'Normalized set point, $\overline{\emph{\textbf{r}}}$','Interpreter','latex');
xlim(ax(1),[0,inf])
plot(ax(1),-1,-1,m1{:},'linewidth',2,'markersize',7)
plot(ax(1),-1,-1,m2{:},'linewidth',2,'markersize',7)
plot(ax(1),-1,-1,'k-','linewidth',2)
plot(ax(1),-1,-1,'k--','linewidth',2)

% J plot
xlabel(ax(2),'Iteration, $k$','Interpreter','latex');
ylabel(ax(2),'Plant profit, $J_p$','Interpreter','latex');
xlim(ax(2),[0,inf])
plot(ax(2),-1,-1,m1{:},'linewidth',2,'markersize',7)
plot(ax(2),-1,-1,'k-','linewidth',2)

% con plot
xlabel(ax(3),'Iteration, $k$','Interpreter','latex');
ylabel(ax(3),'Plant constraint, $\emph{\textbf{G}}_p$','Interpreter','latex');
xlim(ax(3),[0,inf])
ylim(ax(3),[-0.5,0.1])
plot(ax(3),-1,-1,m1{:},'linewidth',2,'markersize',7)
plot(ax(3),-1,-1,m2{:},'linewidth',2,'markersize',7)
plot(ax(3),-1,-1,'k-','linewidth',2)
plot(ax(3),-1,-1,'k--','linewidth',2)

end