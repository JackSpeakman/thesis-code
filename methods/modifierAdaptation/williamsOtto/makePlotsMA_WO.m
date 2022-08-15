% makePlotsMA_WO makes the plots for standard MA on WO used in the 
% literature review
%
% Firstly, the figures will be set-up, along with the model/plant
% functions. Then, the standard MA method is run for K=0.5, which does not
% converge.

% set forPub to 1 for save (and set the default sizes of font/figure/etc.)
forPub = 1;

%% 0. Set-up plots
% clear all vars (other than MAfig and forPub)
clearvars -except MAfig forPub

% create new figures or clear figures if they already exist
% by not deleting already created figures and just clearing, the positions
% do not change which is convenient for running the code multiple times
if exist('MAfig','var') && all(isvalid(MAfig))
    % already exists -> clear
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, MAfig));
    for i = 1:numel(MAfig)
        clf(MAfig(i))
    end
else
    % don't exist -> create
    for i = 1:5
        MAfig(i) = figure;
    end
end

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/williamsOttoCSTR/functions/')
addpath('../')

% set up number of iterations
kmax = 21;

% set-up figure layout on screen
useAllFigPos = 1;

if useAllFigPos || forPub
    if forPub
        % correct sizes
        allFigPos = getAllFigPosPub('WO');
    else
        % nice layout (can see all)
        allFigPos = getAllFigPos('WO');
    end
else
    % don't change from default/current layout
    allFigPos = zeros(5,4);
    for i = 1:5
        allFigPos(i,:) = get(MAfig(i),'Position');
    end
end

% set up figure sizes, etc.
for i = 1:numel(MAfig)
    % change figure
    set(0,'CurrentFigure',MAfig(i))  
    
    % create axis
    ax(i) = axes(MAfig(i));
    
    % turn hold on
    hold(ax(i),'on')
    
    % set x-axis size to number of iterations
    xlim(ax(i),[0,kmax-1]);
    xlabel(ax(i),'Iteration, k','Interpreter','latex')
    
    % run fixAxis to set-up consistent axis
    if forPub
        fixAxis(MAfig(i),ax(i),'linewidth',2.5)
    else
        fixAxis(MAfig(i),ax(i),'linewidth',1,'fontsize',12)
    end
    
    % pull axis to top
    set(ax(i),'Layer','Top')
    
    % set position
    set(MAfig(i),'Position',allFigPos(i,:))
end

% y-axis labels
ylabel(ax(1),'$J$','Interpreter','latex')
ylabel(ax(2),'$X_G$','Interpreter','latex')
ylabel(ax(3),'$F_A$','Interpreter','latex')
ylabel(ax(4),'$F_B$','Interpreter','latex')
ylabel(ax(5),'$T_R$','Interpreter','latex')

% find plant optimum (for plotting)
uGuess = [3.9,9.3,91];
umin = [3,6,80];
umax = [4.5,11,105];

plant = @(u)(WOplantFun(u));
objFun = @(u,y)WOobjFun(u,y);
conFun = @(u,y)WOconFun(u,y);
fminconopts = optimoptions('fmincon','Display','off');

uOptp = fmincon(@(u)objFun(u,plant(u)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,plant(u)),[]),fminconopts);

yOptp = plant(uOptp);
objOptp = objFun(uOptp,yOptp);
conOptp = conFun(uOptp,yOptp);

% set up constraint zero (used for setting up plots)
con0 = conFun([],zeros(1,6));

% legend stuff
cp1 = [0.7,0.7,0.7];
mp1 = {'^-','Color',cp1,'MarkerSize',5,'LineWidth',2.5,'MarkerFaceColor',cp1};

cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',2.5};

for i = 1:numel(MAfig)
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
ylim(ax(3),[3,4.5])
yticks(ax(3),3:0.5:4.5)
ylim(ax(4),[6,11])
yticks(ax(4),6:11)
ylim(ax(5),[80,105])
yticks(ax(5),80:5:105)

% plot plant optimum
area(ax(2),[0,kmax-1],-[con0(1),con0(1)],'BaseValue',1,ma{:})
plot(ax(1),[0,kmax-1],-[objOptp,objOptp],mp{:})
plot(ax(2),[0,kmax-1],conOptp(1)-[con0(1),con0(1)],mp{:})
plot(ax(3),[0,kmax-1],[uOptp(1),uOptp(1)],mp{:})
plot(ax(4),[0,kmax-1],[uOptp(2),uOptp(2)],mp{:})
plot(ax(5),[0,kmax-1],[uOptp(3),uOptp(3)],mp{:})

drawnow

%% 1. Run MA
% set-up functions
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
u0 = [];                                % initial condition
model = @(u)WOmodelFun(u);              % model function
plant = @(u)WOplantFun(u,yGuess);       % plant function
n_u = 3;                                % number of inputs

% run
[ukMA,ykMA,conkMA,objkMA] = runMA('filter',0.5,'kmax',kmax,'conFun',conFun,...
    'startingPoint',u0,'modelFun',model,'plantFun',plant);

% plot
plot(ax(1),0:(kmax-1),-objkMA,mp1{:});
plot(ax(2),0:(kmax-1),conkMA(:,1)-con0(1),mp1{:});
plot(ax(3),0:(kmax-1),ukMA(:,1),mp1{:});
plot(ax(4),0:(kmax-1),ukMA(:,2),mp1{:});
plot(ax(5),0:(kmax-1),ukMA(:,3),mp1{:});

drawnow

%% 2. Save figure

if forPub
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
end

