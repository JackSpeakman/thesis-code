% makePlotsPMAj_WO makes the plots for probaiblistic MA using joint
% consraints on WO used in Chapter 4
% 
% Firstly, the figures will be set-up, along with the model/plant
% functions. Then, standard MA is run at K=0.5, then the PMAj 
% method is run for fixed filter of K=1

% set forPub to 1 for save (and set the default sizes of font/figure/etc.)
forPub = 1;

%% 0. Set-up plots
% 6 plots in total: J, con1, con2, u1, u2, u3
clearvars -except PMAjfig forPub

% create new figures or clear figures if they already exist
% by not deleting already created figures and just clearing, the positions
% do not change which is convenient for running the code multiple times
if exist('PMAjfig','var') && all(isvalid(PMAjfig))
    % already exists -> clear
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, PMAjfig));
    for i = 1:numel(PMAjfig)
        clf(PMAjfig(i))
    end
else
    % don't exist -> create
    for i = 1:6
        PMAjfig(i) = figure;
    end
end

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/williamsOttoCSTR/functions/')
addpath('../')
addpath('../../modifierAdaptation/')

% set up number of iterations
kmax = 11;

% set-up figure layout on screen
useAllFigPos = 1;

if useAllFigPos || forPub
    if forPub
        % correct sizes
        allFigPos = getAllFigPosPub('WO2'); %use DC as 2 constrains
    else
        % nice layout (can see all)
        allFigPos = getAllFigPos('WO2');
    end
else
    % don't change from default/current layout
    allFigPos = zeros(5,4);
    for i = 1:5
        allFigPos(i,:) = get(PMAjfig(i),'Position');
    end
end

% set up figure sizes, etc.
for i = 1:numel(PMAjfig)
    % change figure
    set(0,'CurrentFigure',PMAjfig(i))  
    
    % create axis
    ax(i) = axes(PMAjfig(i));
    
    % turn hold on
    hold(ax(i),'on')
    
    % set x-axis size to number of iterations
    xlim(ax(i),[0,kmax-1]);
    xlabel(ax(i),'Iteration, k','Interpreter','latex')
    
    % run fixAxis to set-up consistent axis
    if forPub
        fixAxis(PMAjfig(i),ax(i),'linewidth',2.5)
    else
        fixAxis(PMAjfig(i),ax(i),'linewidth',1,'fontsize',12)
    end
    
    % pull axis to top
    set(ax(i),'Layer','Top')
    
    % set position
    set(PMAjfig(i),'Position',allFigPos(i,:))
end

% y-axis labels
ylabel(ax(1),'$J$','Interpreter','latex')
ylabel(ax(2),'$X_G$','Interpreter','latex')
ylabel(ax(3),'$X_B-X_E$','Interpreter','latex')
ylabel(ax(4),'$F_A$','Interpreter','latex')
ylabel(ax(5),'$F_B$','Interpreter','latex')
ylabel(ax(6),'$T_R$','Interpreter','latex')

% find plant optimum (for plotting)
uGuess = [3.9,9.3,91];
umin = [3,6,80];
umax = [4.5,11,105];

plant = @(u)(WOplantFun(u));
objFun = @(u,y)WOobjFun(u,y);
conFun = @(u,y)WOconFun2(u,y);
fminconopts = optimoptions('fmincon','Display','off');

uOptp = fmincon(@(u)objFun(u,plant(u)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,plant(u)),[]),fminconopts);

yOptp = plant(uOptp);
objOptp = objFun(uOptp,yOptp);
conOptp = conFun(uOptp,yOptp);% set up constraint zero (used for setting up plots)
con0 = conFun([],zeros(1,6));

% MA line
cp1 = [0.7,0.7,0.7];
mp1 = {'^-','Color',cp1,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp1};

% WCMA line
cp2 = brightGreen;
mp2 = {'d-','Color',cp2,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp2};

% infeasible region
cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',2.5};

for i = 1:numel(PMAjfig)
    % plot these off the figure first in the order of the desired legend
    plot(ax(i),-1,-1,mp1{:});
    plot(ax(i),-1,-1,mp{:});
    area(ax(i),-1,-1,ma{:});
end

% set-up yaxis for each figure
ylim(ax(1),[100,220])
yticks(ax(1),100:20:220)
ylim(ax(2),[0.045,0.09])
yticks(ax(2),0.05:0.01:0.09)
ylim(ax(3),[0.14,0.26])
yticks(ax(3),0.14:0.02:0.26)
ylim(ax(4),[3,4.5])
yticks(ax(4),3:0.5:4.5)
ylim(ax(5),[6,11])
yticks(ax(5),6:11)
ylim(ax(6),[80,105])
yticks(ax(6),80:5:105)

% plot plant optimum
plot(ax(1),[0,kmax-1],-[objOptp,objOptp],mp{:})
area(ax(2),[0,kmax-1],-[con0(1),con0(1)],'BaseValue',1,ma{:})
plot(ax(2),[0,kmax-1],conOptp(1)-[con0(1),con0(1)],mp{:})
area(ax(3),[0,kmax-1],[con0(2),con0(2)],'BaseValue',0,ma{:})
plot(ax(3),[0,kmax-1],-conOptp(2)+[con0(2),con0(2)],mp{:})
plot(ax(4),[0,kmax-1],[uOptp(1),uOptp(1)],mp{:})
plot(ax(5),[0,kmax-1],[uOptp(2),uOptp(2)],mp{:})
plot(ax(6),[0,kmax-1],[uOptp(3),uOptp(3)],mp{:})

drawnow

%% 1. Run MA
% set-up functions
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
u0 = [];                                % initial condition
model = @(u)WOmodelFun(u);              % model function
plant = @(u)WOplantFun(u,yGuess);       % plant function
n_u = 3;                                % number of inputs

% run
[ukMA,ykMA,conkMA,objkMA] = runMA(...
    'filter',0.5,...        % RTO input filter gain
    'kmax',kmax,...         % maximum number of iterates
    'conFun',conFun,...     % constraint function @(r,y)
    'objFun',objFun,...     % objective function @(r,y)
    'startingPoint',u0,...  % initial conditions for RTO
    'modelFun',model,...    % model function @(u)
    'plantFun',plant);      % plant function @(u)

% plot
plot(ax(1),0:(kmax-1),-objkMA,mp1{:});
plot(ax(2),0:(kmax-1),conkMA(:,1)-con0(1),mp1{:});
plot(ax(3),0:(kmax-1),-conkMA(:,2)+con0(2),mp1{:});
plot(ax(4),0:(kmax-1),ukMA(:,1),mp1{:});
plot(ax(5),0:(kmax-1),ukMA(:,2),mp1{:});
plot(ax(6),0:(kmax-1),ukMA(:,3),mp1{:});

drawnow

%% 2. RunPMAj
% set-up functions
model = @(u,th)WOmodelFun(u,yGuess,th); % model function

% run
[ukPMAj,ykPMAj,conkPMAj,objkPMAj] = runRobustMA(...
    'method','PMAj',...     % robust MA method
    'filter',1,...        % RTO input filter gain
    'kmax',kmax,...         % maximum number of iterates
    'conFun',conFun,...     % constraint function @(r,y)
    'objFun',objFun,...     % objective function @(r,y)
    'startingPoint',u0,...  % initial conditions for RTO
    'modelFun',model,...    % model function @(u)
    'plantFun',plant,...    % plant function @(u)
    'probChance',0.8);      % probability constraint chance

% plot
plot(ax(1),0:(kmax-1),-objkPMAj,mp2{:});
plot(ax(2),0:(kmax-1),conkPMAj(:,1)-con0(1),mp2{:});
plot(ax(3),0:(kmax-1),-conkPMAj(:,2)+con0(2),mp2{:});
plot(ax(4),0:(kmax-1),ukPMAj(:,1),mp2{:});
plot(ax(5),0:(kmax-1),ukPMAj(:,2),mp2{:});
plot(ax(6),0:(kmax-1),ukPMAj(:,3),mp2{:});

drawnow

%% 3. Save figures

if forPub
    saveas(PMAjfig(1),'plots\PMAjobj_WO.eps','epsc')
    saveas(PMAjfig(1),'plots\PMAjobj_WO.fig','fig')
    
    saveas(PMAjfig(2),'plots\PMAjcon1_WO.eps','epsc')
    saveas(PMAjfig(2),'plots\PMAjcon1_WO.fig','fig')
    
    saveas(PMAjfig(3),'plots\PMAjcon2_WO.eps','epsc')
    saveas(PMAjfig(3),'plots\PMAjcon2_WO.fig','fig')
    
    saveas(PMAjfig(4),'plots\PMAju1_WO.eps','epsc')
    saveas(PMAjfig(4),'plots\PMAju1_WO.fig','fig')
    
    saveas(PMAjfig(5),'plots\PMAju2_WO.eps','epsc')
    saveas(PMAjfig(5),'plots\PMAju2_WO.fig','fig')
    
    saveas(PMAjfig(6),'plots\PMAju3_WO.eps','epsc')
    saveas(PMAjfig(6),'plots\PMAju3_WO.fig','fig')
end

