% makes the plots for standard MA on WO used in the literature review
clearvars -except MAfig

if exist('MAfig','var') && all(isvalid(MAfig))
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, MAfig));
    for i = 1:numel(MAfig)
        clf(MAfig(i))
    end
else
    for i = 1:5
        MAfig(i) = figure;
    end
end

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/williamsOttoCSTR/functions/')

% set up constraint
conFun = @(u,y)WOconFun(u,y);
con0 = conFun([],zeros(1,6));

% set up number of iterations
kmax = 21;

% set up figure
set(0,'CurrentFigure',MAfig(1))

for i = 1:numel(MAfig)
    set(0,'CurrentFigure',MAfig(i))
    ax(i) = axes(MAfig(i));
    hold(ax(i),'on')
    xlim(ax(i),[0,kmax-1]);
    
    fixAxis(MAfig(i),ax(i),'linewidth',2.5)
    set(ax(i),'Layer','Top')
    
    set(MAfig(i),'Position',[50+50*i,50+50*i,800,600])
    xlabel(ax(i),'Iteration, k','Interpreter','latex')
end

set(MAfig(1),'Position',[50+50*i,50+50*i,800,600])
set(MAfig(2),'Position',[50+50*i,50+50*i,800,600])
set(MAfig(3),'Position',[50+50*i,50+50*i,800*2/3,600*3/4])
set(MAfig(4),'Position',[50+50*i,50+50*i,800*2/3,600*3/4])
set(MAfig(5),'Position',[50+50*i,50+50*i,800*2/3,600*3/4])

ylabel(ax(1),'$J$','Interpreter','latex')
ylabel(ax(2),'$X_G$','Interpreter','latex')
ylabel(ax(3),'$F_A$','Interpreter','latex')
ylabel(ax(4),'$F_B$','Interpreter','latex')
ylabel(ax(5),'$T_R$','Interpreter','latex')

% find plant optimum
uGuess = [3.9,9.3,91];
umin = [3,6,80];
umax = [4.5,11,105];

plant = @(u)(WOplantFun(u));
objFun = @(u,y)WOobjFun(u,y);
fminconopts = optimoptions('fmincon','Display','off');

uOptp = fmincon(@(u)objFun(u,plant(u)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,plant(u)),[]),fminconopts);

yOptp = plant(uOptp);
objOptp = objFun(uOptp,yOptp);
conOptp = conFun(uOptp,yOptp);


% legend stuff
cp1 = [0.7,0.7,0.7];
mp1 = {'^-','Color',cp1,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp1};

cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',2.5};

for i = 1:numel(MAfig)
    plot(ax(i),-1,-1,mp1{:});
    plot(ax(i),-1,-1,mp{:});
    area(ax(i),-1,-1,ma{:});
end

% plant plot
ylim(ax(1),[100,220])
yticks(ax(1),100:20:220)
ylim(ax(2),[0.045,0.09])
yticks(ax(2),0.05:0.01:0.09)
ylim(ax(3),[3,4.5])
yticks(ax(3),3:0.5:4.5)
ylim(ax(4),[6,11])
yticks(ax(4),6:11)
ylim(ax(5),[80,105])
yticks(ax(5),80:5:105)

area(ax(2),[0,kmax-1],-[con0(1),con0(1)],'BaseValue',1,ma{:})

plot(ax(1),[0,kmax-1],-[objOptp,objOptp],mp{:})
plot(ax(2),[0,kmax-1],conOptp(1)-[con0(1),con0(1)],mp{:})
plot(ax(3),[0,kmax-1],[uOptp(1),uOptp(1)],mp{:})
plot(ax(4),[0,kmax-1],[uOptp(2),uOptp(2)],mp{:})
plot(ax(5),[0,kmax-1],[uOptp(3),uOptp(3)],mp{:})

drawnow

%% 1. Standard MA
% set-up functions
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
u0 = [];                                % initial condition
model = @(u)WOmodelFun([u],yGuess);     % model function
plant = @(u)WOplantFun([u],yGuess);     % plant function
n_u = 3;                                % number of inputs

% run
[ukMA,ykMA,conkMA,objkMA] = runMA_WO('filter',0.3,'kmax',kmax,'conFun',conFun,...
    'startingPoint',u0,'modelFun',model,'plantFun',plant,'num_inputs',n_u);

% plot
plot(ax(1),0:(kmax-1),-objkMA,mp1{:});
plot(ax(2),0:(kmax-1),conkMA(:,1)-con0(1),mp1{:});
plot(ax(3),0:(kmax-1),ukMA(:,1),mp1{:});
plot(ax(4),0:(kmax-1),ukMA(:,2),mp1{:});
plot(ax(5),0:(kmax-1),ukMA(:,3),mp1{:});

drawnow

%% 5. Save figure
% legend
% leg = {'Standard MA\quad','WCMA\quad','PMAi\quad','PAMj\quad','Plant Optimum\quad','Infeasible Region'};
% legend(ax(3,2),leg,'Location','southoutside','Orientation','horizontal','Interpreter','latex');

saveas(MAfig(1),'plots\MAobj_WO.eps','epsc')
saveas(MAfig(1),'plots\MAobj_WO.fig','fig')

saveas(MAfig(2),'plots\MAcon_WO.eps','epsc')
saveas(MAfig(2),'plots\MAcon_WO.fig','fig')

saveas(MAfig(3),'plots\MAu1_WO.eps','epsc')
saveas(MAfig(3),'plots\MAu1_WO.fig','fig')

saveas(MAfig(4),'plots\MAu2_WO.eps','epsc')
saveas(MAfig(4),'plots\MAu2_WO.fig','fig')

saveas(MAfig(5),'plots\MAu3_WO.eps','epsc')
saveas(MAfig(5),'plots\MAu3_WO.fig','fig')

