% script to run robust MA on WO and plot results
clearvars -except fig
if exist('fig','var') && isvalid(fig)
    set(0,'CurrentFigure',fig);
    clf
else
    fig = figure;
end
% add path
addpath('../General Functions/');
addpath('../General Functions/Plot/');

% set up constraint
conFun = @(u,y)([y(:,6)-0.08,y(:,5)-y(:,2)+0.16]);
con0 = conFun([],zeros(1,6));

% set up number of iterations
kmax = 21;

% set up figure
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');
for i = 1:3
    for j = 1:3
        ax(i,j) = nexttile;
        hold on
        xlim(ax(i,j),[0,kmax-1]);
        
        fixAxis(fig,ax(i,j),'linewidth',2.5)
        set(ax(i,j),'Layer','Top')
    end
end

set(fig,'Position',[-3800,300,1800,1000])

for j = 1:3
    xlabel(ax(3,j),'Iteration, k','Interpreter','latex')
end

title(ax(1,1),'\bf{(a) WCMA}','Interpreter','latex')
title(ax(1,2),'\bf{(b) PMAi}','Interpreter','latex')
title(ax(1,3),'\bf{(c) PMAj}','Interpreter','latex')

ylabel(ax(1,1),'$$J$$','Interpreter','latex')
ylabel(ax(2,1),'$$X_G$$','Interpreter','latex')
ylabel(ax(3,1),'$$X_B - X_E$$','Interpreter','latex')

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
cp2 = [1,0,0];
cp3 = [0,0,1];
cp4 = [0.1,0.7,0.3];

mp1 = {'^-','Color',cp1,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp1};
mp2 = {'o-','Color',cp2,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp2};
mp3 = {'s-','Color',cp3,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp3};
mp4 = {'d-','Color',cp4,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp4};

cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',2.5};

plot(ax(3,2),-1,-1,mp1{:});
plot(ax(3,2),-1,-1,mp2{:});
plot(ax(3,2),-1,-1,mp3{:});
plot(ax(3,2),-1,-1,mp4{:});
plot(ax(3,2),-1,-1,mp{:});
area(ax(3,2),-1,-1,ma{:});


% plant plot
for j = 1:3
    ylim(ax(1,j),[100,220])
    yticks(ax(1,j),100:20:220)
    ylim(ax(2,j),[0.045,0.09])
    yticks(ax(2,j),0.05:0.01:0.09)
    ylim(ax(3,j),[0.14,0.26])
    yticks(ax(3,j),0.14:0.02:0.26)
    
    area(ax(2,j),[0,kmax-1],-[con0(1),con0(1)],'BaseValue',1,ma{:})
    area(ax(3,j),[0,kmax-1],[con0(2),con0(2)],'BaseValue',-1,ma{:})

    plot(ax(1,j),[0,kmax-1],-[objOptp,objOptp],mp{:})
    plot(ax(2,j),[0,kmax-1],conOptp(1)-[con0(1),con0(1)],mp{:})
    plot(ax(3,j),[0,kmax-1],-conOptp(2)+[con0(2),con0(2)],mp{:})
end

drawnow

%% 1. Standard MA
% add path
addpath('../Standard MA/');
addpath('../Case Studies/Willaims Otto Case/');

kmax = 21;
% run
[ukMA,ykMA,conkMA,objkMA] = runWO_MA('K',0.5,'kmax',kmax,...
    'conFun',conFun);

% plot

for j = 1:3
    plot(ax(1,j),0:(kmax-1),-objkMA,mp1{:});
    plot(ax(2,j),0:(kmax-1),conkMA(:,1)-con0(1),mp1{:});
    plot(ax(3,j),0:(kmax-1),-conkMA(:,2)+con0(2),mp1{:});
end
drawnow

%% 2. Worst-Case MA
% set theta
th = [0,0;70,160;-70,160;70,-160;-70,-160]*3; 

% set WC filter
n_r = 21;
filterWC = @(a,b,c)newFilter(a,b,c,n_r,1);

% set quadratic
Qk = [0.000001,0.000001,0.000001];

% run
[ukWC,ykWC,conkWC,objkWC] = runWO_WCMA('K',filterWC,'kmax',kmax,...
    'conFun',conFun,'th',th,'Qk',Qk);

% plot
plot(ax(1,1),0:(kmax-1),-objkWC,mp2{:});
plot(ax(2,1),0:(kmax-1),conkWC(:,1)-con0(1),mp2{:});
plot(ax(3,1),0:(kmax-1),-conkWC(:,2)+con0(2),mp2{:});

drawnow

%% 3. PMA individual
% generate theta
rng(100)
n_th = 100;
mu = [0,0];
s = [70,160];
corr = 0.85;
sigma = [s(1)*s(1),corr*s(2)*s(1);corr*s(1)*s(2),s(2)*s(2)];
th = [0,0;mvnrnd(mu,sigma,n_th-1)];

% set PMA filter
n_r = 21;
pK = 0.9;
filterPMA = @(a,b,c)newFilter(a,b,c,n_r,pK);

% set probability
pj = 0.9;

% run
[ukPi,ykPi,conkPi,objkPi] = runWO_PMAi('K',filterPMA,'kmax',kmax,...
    'conFun',conFun,'th',th,'Qk',Qk,'p',pj);

% plot
plot(ax(1,2),0:(kmax-1),-objkPi,mp3{:});
plot(ax(2,2),0:(kmax-1),conkPi(:,1)-con0(1),mp3{:});
plot(ax(3,2),0:(kmax-1),-conkPi(:,2)+con0(2),mp3{:});

drawnow

%% 4. PMA joint
% probability
p = 0.9;

% run
[ukPj,ykPj,conkPj,objkPj] = runWO_PMAj('K',filterPMA,'kmax',kmax,...
    'conFun',conFun,'th',th,'Qk',Qk,'p',p);

% plot
plot(ax(1,3),0:(kmax-1),-objkPj,mp4{:});
plot(ax(2,3),0:(kmax-1),conkPj(:,1)-con0(1),mp4{:});
plot(ax(3,3),0:(kmax-1),-conkPj(:,2)+con0(2),mp4{:});

drawnow

%% 5. Save figure
% legend
leg = {'Standard MA\quad','WCMA\quad','PMAi\quad','PAMj\quad','Plant Optimum\quad','Infeasible Region'};
legend(ax(3,2),leg,'Location','southoutside','Orientation','horizontal','Interpreter','latex');

saveas(fig,'Plots\WO_full.eps','epsc')
saveas(fig,'Plots\WO_full.fig','fig')
