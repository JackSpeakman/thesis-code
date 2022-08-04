% script to run worst case MA on WO and plot results
clearvars -except PMAifig

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/williamsOttoCSTR/functions/')
addpath('../../modifierAdaptation/williamsOtto/')

% set up constraint
conFun = @(u,y)WOconFun2(u,y);
con0 = conFun([],zeros(1,6));

n_c = numel(con0);

% create figures
if exist('PMAifig','var') && all(isvalid(PMAifig)) && numel(PMAifig) == 4+n_c
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, PMAifig));
    for i = 1:numel(PMAifig)
        clf(PMAifig(i))
    end
else
    for i = 1:(4+n_c)
        PMAifig(i) = figure;
    end
end

% set up number of iterations
kmax = 21;

% set up figure
set(0,'CurrentFigure',PMAifig(1))

for i = 1:numel(PMAifig)
    set(0,'CurrentFigure',PMAifig(i))
    ax(i) = axes(PMAifig(i));
    hold(ax(i),'on')
    xlim(ax(i),[0,kmax-1]);
    
    fixAxis(PMAifig(i),ax(i),'linewidth',2.5)
    set(ax(i),'Layer','Top')
    
    set(PMAifig(i),'Position',[50+50*i,50+50*i,800,600])
    xlabel(ax(i),'Iteration, k','Interpreter','latex')

end

set(PMAifig(1),'Position',[50+50*i,-150+50*i,800,600])
set(PMAifig(2),'Position',[50+50*i,-150+50*i,800,600])
set(PMAifig(3),'Position',[50+50*i,-150+50*i,800,600])
set(PMAifig(n_c+2),'Position',[50+50*i,-150+50*i,800*2/3,600*3/4])
set(PMAifig(n_c+3),'Position',[50+50*i,-150+50*i,800*2/3,600*3/4])
set(PMAifig(n_c+4),'Position',[50+50*i,-150+50*i,800*2/3,600*3/4])

ylabel(ax(1),'$J$','Interpreter','latex')
ylabel(ax(2),'$X_G$','Interpreter','latex')
ylabel(ax(3),'$X_B-X_E$','Interpreter','latex')
ylabel(ax(n_c+2),'$F_A$','Interpreter','latex')
ylabel(ax(n_c+3),'$F_B$','Interpreter','latex')
ylabel(ax(n_c+4),'$T_R$','Interpreter','latex')

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

cp2 = [0,0,1];
mp2 = {'s-','Color',cp2,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp2};

cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',2.5};

for i = 1:numel(PMAifig)
    plot(ax(i),-1,-1,mp1{:});
    plot(ax(i),-1,-1,mp{:});
    area(ax(i),-1,-1,ma{:});
end

% plant plot
ylim(ax(1),[100,220])
yticks(ax(1),100:20:220)
ylim(ax(2),[0.045,0.09])
yticks(ax(2),0.05:0.01:0.09)
ylim(ax(3),[0.14,0.26])
yticks(ax(3),0.14:0.02:0.26)
ylim(ax(n_c+2),[3,4.5])
yticks(ax(n_c+2),3:0.5:4.5)
ylim(ax(n_c+3),[6,11])
yticks(ax(n_c+3),6:11)
ylim(ax(n_c+4),[80,105])
yticks(ax(n_c+4),80:5:105)

% obj
plot(ax(1),[0,kmax-1],-[objOptp,objOptp],mp{:})

% con
area(ax(2),[0,kmax-1],-[con0(1),con0(1)],'BaseValue',1,ma{:})
plot(ax(2),[0,kmax-1],conOptp(1)-[con0(1),con0(1)],mp{:})
if n_c == 2
    area(ax(3),[0,kmax-1],[con0(2),con0(2)],'BaseValue',0,ma{:})
    plot(ax(3),[0,kmax-1],[con0(2),con0(2)]-conOptp(2),mp{:})
end

% inputs
plot(ax(n_c+2),[0,kmax-1],[uOptp(1),uOptp(1)],mp{:})
plot(ax(n_c+3),[0,kmax-1],[uOptp(2),uOptp(2)],mp{:})
plot(ax(n_c+4),[0,kmax-1],[uOptp(3),uOptp(3)],mp{:})

drawnow

%% 1. Standard MA
% set-up functions
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];

% run
[ukMA,ykMA,conkMA,objkMA] = runMA_WO('filter',0.5,'kmax',kmax,'conFun',conFun,'num_inputs',3);

% plot
plot(ax(1),0:(kmax-1),-objkMA,mp1{:});
plot(ax(2),0:(kmax-1),conkMA(:,1)-con0(1),mp1{:});
plot(ax(3),0:(kmax-1),con0(n_c)-conkMA(:,n_c),mp1{:});
plot(ax(n_c+2),0:(kmax-1),ukMA(:,1),mp1{:});
plot(ax(n_c+3),0:(kmax-1),ukMA(:,2),mp1{:});
plot(ax(n_c+4),0:(kmax-1),ukMA(:,3),mp1{:});

drawnow

%% 2. Probabilistic MA
% set-up functions
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];

% parameters
th = [0,0; 70,160; -70,160; 70,-160; -70,-160]*3;

% run
[ukPMAi,ykPMAi,conkPMAi,objkPMAi] = runPMAi_WO('filter',1,'kmax',kmax,...
    'conFun',conFun,'th',th,'constraintChance',0.95);

% plot
plot(ax(1),0:(kmax-1),-objkPMAi,mp2{:});
plot(ax(2),0:(kmax-1),conkPMAi(:,1)-con0(1),mp2{:});
plot(ax(3),0:(kmax-1),con0(n_c)-conkPMAi(:,n_c),mp2{:});
plot(ax(n_c+2),0:(kmax-1),ukPMAi(:,1),mp2{:});
plot(ax(n_c+3),0:(kmax-1),ukPMAi(:,2),mp2{:});
plot(ax(n_c+4),0:(kmax-1),ukPMAi(:,3),mp2{:});

drawnow

%% 5. Save figures
saveas(PMAifig(1),'plots\PMAiobj_WO.eps','epsc')
saveas(PMAifig(1),'plots\PMAiobj_WO.fig','fig')

saveas(PMAifig(2),'plots\PMAicon_WO.eps','epsc')
saveas(PMAifig(2),'plots\PMAicon_WO.fig','fig')

if n_c == 2
    saveas(PMAifig(3),'plots\PMAicon2_WO.eps','epsc')
    saveas(PMAifig(3),'plots\PMAicon2_WO.fig','fig')
end

saveas(PMAifig(n_c+2),'plots\PMAiu1_WO.eps','epsc')
saveas(PMAifig(n_c+2),'plots\PMAiu1_WO.fig','fig')

saveas(PMAifig(n_c+2),'plots\PMAiu2_WO.eps','epsc')
saveas(PMAifig(n_c+2),'plots\PMAiu2_WO.fig','fig')

saveas(PMAifig(n_c+2),'plots\PMAiu3_WO.eps','epsc')
saveas(PMAifig(n_c+2),'plots\PMAiu3_WO.fig','fig')
