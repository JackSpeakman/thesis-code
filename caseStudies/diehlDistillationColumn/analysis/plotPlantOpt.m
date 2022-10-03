function plotPlantOpt(ax,kmax,rpOpt,objpOpt,conpOpt)
% plots the plant optimum lines to axes
plot(ax(1),[0,kmax-1],[rpOpt(1),rpOpt(1)]/4-70/4,'k-','linewidth',1)
plot(ax(1),[0,kmax-1],[rpOpt(2),rpOpt(2)]/4-82/4,'k--','linewidth',1)
plot(ax(2),[0,kmax-1],-[objpOpt,objpOpt],'k-','linewidth',1)
plot(ax(3),[0,kmax-1],[conpOpt(1),conpOpt(1)],'k-','linewidth',1)
plot(ax(3),[0,kmax-1],[conpOpt(2),conpOpt(2)],'k--','linewidth',1)

end

