% makes the plots for DHMA - 3 variable used in the thesis
clearvars -except DHMAfig

if exist('DHMAfig','var') && all(isvalid(DHMAfig))
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, DHMAfig));
    for i = 1:numel(DHMAfig)
        clf(DHMAfig(i))
    end
else
    for i = 1:5
        DHMAfig(i) = figure;
    end
end

% add path
addpath('../Colours/')
addpath('../General Functions/');
addpath('../General Functions/Plot/');
addpath('../Case Studies/Willaims Otto Case/')
addpath('../Standard MA/')

% set up constraint
conFun = @(u,y)([y(:,6)-0.08]);
con0 = conFun([],zeros(1,6));

% set up number of iterations
kmax = 21;

% set up figures
set(0,'CurrentFigure',DHMAfig(1))

wasDocked = zeros(size(DHMAfig));
for i = 1:numel(DHMAfig)
    if ~strcmp(get(DHMAfig(i),'WindowStyle'),'normal')
        set(DHMAfig(i),'WindowStyle','normal')
        wasDocked(i) = 1;
    end
    set(DHMAfig(i),'Position',[50+50*i,50+50*i,800,600])
end

for i = 1:numel(DHMAfig)
    set(0,'CurrentFigure',DHMAfig(1))
    ax(i) = axes(DHMAfig(i));
    hold(ax(i),'on')
    xlim(ax(i),[0,kmax-1]);
    
    fixAxis(DHMAfig(i),ax(i),'linewidth',2.5)
    set(ax(i),'Layer','Top')
    
    set(DHMAfig(i),'Position',[50+50*i,50+50*i,800,600])
    xlabel(ax(i),'Iteration, k','Interpreter','latex')
end

for i = 1:numel(DHMAfig)
    set(DHMAfig(i),'Position',[50+50*i,50+50*i,800,500])
    if wasDocked(i)
        set(DHMAfig(i),'WindowStyle','docked')
    end
end

ylabel(ax(1),'$J$','Interpreter','latex')
ylabel(ax(2),'$X_G$','Interpreter','latex')
ylabel(ax(3),'$F_A$','Interpreter','latex')
ylabel(ax(4),'$F_B$','Interpreter','latex')
ylabel(ax(5),'$T_R$','Interpreter','latex')

% find plant optimum
uGuess = [3.9,9.3,91];

umin = [3,6,80];
umax = [4.5,11,105];

uNorm = @(u)((u-umin)./(umax-umin));
uRest = @(u)(umin+(umax-umin).*u);

plant = @(u)(WOplantFun(uRest(u)));
objFun = @(u,y)WOobjFun(uRest(u),y);
fminconopts = optimoptions('fmincon','Display','off');

uOptp = fmincon(@(u)objFun(u,plant(u)),uNorm(uGuess),[],[],[],[],[0,0,0],[1,1,1],...
    @(u)deal(conFun(u,plant(u)),[]),fminconopts);

uOptpR = uRest(uOptp);
yOptp = plant(uOptp);
objOptp = objFun(uOptp,yOptp);
conOptp = conFun(uOptp,yOptp);

% legend stuff
cp1 = [0.7,0.7,0.7];
cp2 = [1,0,0];
cp3 = [0,0,1];
cp4 = [0.1,0.7,0.3];
cp5 = brightBlue;

mp1 = {'^-','Color',cp1,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp1};
mp2 = {'.-','Color',cp5,'MarkerSize',20,'LineWidth',2.5,'MarkerFaceColor',cp5};

cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',2.5};

for i = 1:numel(DHMAfig)
    plot(ax(i),-1,-1,mp1{:});
    plot(ax(i),-1,-1,mp2{:});
    plot(ax(i),-1,-1,mp{:});
    area(ax(i),-1,-1,ma{:});
end

% plant plot
ylim(ax(1),[100,220])
yticks(ax(1),100:20:220)
ylim(ax(2),[0.05,0.09])
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
plot(ax(3),[0,kmax-1],[uOptpR(1),uOptpR(1)],mp{:})
plot(ax(4),[0,kmax-1],[uOptpR(2),uOptpR(2)],mp{:})
plot(ax(5),[0,kmax-1],[uOptpR(3),uOptpR(3)],mp{:})

drawnow

%% 1. Standard MA
% run
[ukMA,ykMA,conkMA,objkMA] = runWO_MA('K',0.5,'kmax',kmax);
ukMAr = ukMA;

% plot
plot(ax(1),0:(kmax-1),-objkMA,mp1{:});
plot(ax(2),0:(kmax-1),conkMA(:,1)-con0(1),mp1{:});
plot(ax(3),0:(kmax-1),ukMAr(:,1),mp1{:});
plot(ax(4),0:(kmax-1),ukMAr(:,2),mp1{:});
plot(ax(5),0:(kmax-1),ukMAr(:,3),mp1{:});

drawnow

%% 2. Directional-Hessian MA
% run
[ukDHMA,ykDHMA,conkDHMA,objkDHMA] = runWO3_DHMA('K',0.9,'kmax',kmax);
ukDHMAr = uRest(ukDHMA);

% plot
plot(ax(1),0:(kmax-1),-objkDHMA,mp2{:});
plot(ax(2),0:(kmax-1),conkDHMA(:,1)-con0(1),mp2{:});
plot(ax(3),0:(kmax-1),ukDHMAr(:,1),mp2{:});
plot(ax(4),0:(kmax-1),ukDHMAr(:,2),mp2{:});
plot(ax(5),0:(kmax-1),ukDHMAr(:,3),mp2{:});

drawnow

%% 3. Save figure
% legend
leg = {'Standard MA\quad','DHMA\quad','Plant Optimum\quad','Infeasible Region'};
legend(ax(2),leg,'Location','southeast','Interpreter','latex');

saveas(DHMAfig(1),'Plots\DHMAobj_3var.eps','epsc')
saveas(DHMAfig(1),'Plots\DHMAobj_3var.fig','fig')

saveas(DHMAfig(2),'Plots\DHMAcon_3var.eps','epsc')
saveas(DHMAfig(2),'Plots\DHMAcon_3var.fig','fig')

saveas(DHMAfig(3),'Plots\DHMAu1_3var.eps','epsc')
saveas(DHMAfig(3),'Plots\DHMAu1_3var.fig','fig')

saveas(DHMAfig(4),'Plots\DHMAu2_3var.eps','epsc')
saveas(DHMAfig(4),'Plots\DHMAu2_3var.fig','fig')

saveas(DHMAfig(5),'Plots\DHMAu3_3var.eps','epsc')
saveas(DHMAfig(5),'Plots\DHMAu3_3var.fig','fig')


