% Makes the plots for WCMA on DC used in the comparison chapter
% Firstly, the figures will be set-up, along with the model/plant
% functions. Then, the WCMA method is run for K=0.9

% set forPub to 1 for save (and set the default sizes of font/figure/etc.)
forPub = 1;


%% 0. Set-up plots
% clear all vars (other than MAfig and forPub)
clearvars -except fig

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/diehlDistillationColumn/functions/')
addpath('../../../caseStudies/diehlDistillationColumn/analysis/')
addpath('../')
addpath('../')

%% 0.a Set-up
r0 = [70.001,85];   % starting point
k_max = 11;         % number of iterations
Ns = 14;            % number of stripping trays
Nr = 7:11;          % number of rectifying trays
mu = 3;             % which rectifying tray is nominal

% lines
c_MA = [0.6,0.6,0.6];
m1_MA = {'-o','color',c_MA,'linewidth',3,'markersize',7,'MarkerFaceColor',c_MA};
m2_MA = {'--o','color',c_MA,'linewidth',3,'markersize',7,'MarkerFaceColor',c_MA};

%% 0.b Get plant optimum
fprintf('\nFinding plant optimum\n');
tic
[rpOpt,upOpt,ypOpt,objpOpt,conpOpt] = getPlantOpt;
t = toc;
fprintf('Plant optimum found [%4.2f s]\n\n',t);

%% 1 Set-up DC functions
% Set-up variables to be given to runMA. The decision variables are the 
% controller setpoints (r). four functions are required (model(r), 
% plant(r), objFun(r,y), conFun(r,y)). The cost/cons functions are defined 
% in terms of the input variables (u). PIcontoller relates r/y to u. The
% following code sets up the distillation column for all models (variable
% size of rectifier, but only uses the nominal (Nr = 9) for the standard
% MA run.

[model_ri,plant,objFun,conFun] = setupDCshortFunctions;

% replace model_ri(r,i) with model(r)
model = @(r)model_ri(r,mu);

%% 1.a Run standard MA (K=0.9)
% run at K=0.9;
[rk_MA1,yk_MA1,conk_MA1,objk_MA1] = runMA(...
    'startingPoint',r0,...  % initial conditions for RTO
    'kmax',k_max,...        % maximum number of iterates
    'filter',0.9,...        % RTO input filter gain
    'conFun',conFun,...     % constraint function @(r,y)
    'objFun',objFun,...     % objective function @(r,y)
    'modelFun',model,...    % model function @(r)
    'plantFun',plant,...    % plant function @(r)
    'umin',[68,82],...      % fmincon minimum limit
    'umax',[74,86]);        % fmincon maximum limit

%% 1.b Plot standard MA
if exist('fig','var')
    newFig = 0;
else
    newFig = 1;
end

for i = 1:3
    if newFig || ~isvalid(fig{i})
        fig{i} = figure;
    else
        clf(fig{i})
        ax{i} = [];
    end
    ax{i} = axes(fig{i});
    hold(ax{i},'on')
    fixAxis(fig{i},ax{i})
    p = fig{i}.Position;
    set(fig{i},'Position',[p(1:2),p(3:4)-100])
end

setupPlotsDC(ax,m1_MA,m2_MA)
plotPlantOpt(ax,k_max,rpOpt,objpOpt,conpOpt)

plot(ax{1},0:k_max-1,(rk_MA1(:,1)-70)/4,m1_MA{:})
plot(ax{1},0:k_max-1,(rk_MA1(:,2)-82)/4,m2_MA{:})
plot(ax{2},0:k_max-1,-objk_MA1,m1_MA{:})
plot(ax{3},0:k_max-1,conk_MA1(:,1)*100,m1_MA{:})
plot(ax{3},0:k_max-1,conk_MA1(:,2)*100,m2_MA{:})

setupLegendDC(ax)

saveas(fig{1},'plots/MAr1_DCshort.eps','epsc')
saveas(fig{2},'plots/MAobj1_DCshort.eps','epsc')
saveas(fig{3},'plots/MAcon1_DCshort.eps','epsc')

%% 1.c Run standard MA (K=0.6)
% run at K=0.6;
[rk_MA2,yk_MA2,conk_MA2,objk_MA2] = runMA(...
    'startingPoint',r0,...  % initial conditions for RTO
    'kmax',k_max,...        % maximum number of iterates
    'filter',0.6,...        % RTO input filter gain
    'conFun',conFun,...     % constraint function @(r,y)
    'objFun',objFun,...     % objective function @(r,y)
    'modelFun',model,...    % model function @(r)
    'plantFun',plant,...    % plant function @(r)
    'umin',[68,82],...      % fmincon minimum limit
    'umax',[74,86]);        % fmincon maximum limit

%% 1.d Plot standard MA (K=0.6)
if exist('fig','var')
    newFig = 0;
else
    newFig = 1;
end

for i = 1:3
    if newFig || ~isvalid(fig{i})
        fig{i} = figure;
    else
        clf(fig{i})
        ax{i} = [];
    end
    ax{i} = axes(fig{i});
    hold(ax{i},'on')
    fixAxis(fig{i},ax{i})
    p = fig{i}.Position;
    set(fig{i},'Position',[p(1:2),p(3:4)-100])
end

setupPlotsDC(ax,m1_MA,m2_MA)
plotPlantOpt(ax,k_max,rpOpt,objpOpt,conpOpt)

plot(ax{1},0:k_max-1,(rk_MA2(:,1)-70)/4,m1_MA{:})
plot(ax{1},0:k_max-1,(rk_MA2(:,2)-82)/4,m2_MA{:})
plot(ax{2},0:k_max-1,-objk_MA2,m1_MA{:})
plot(ax{3},0:k_max-1,conk_MA2(:,1)*100,m1_MA{:})
plot(ax{3},0:k_max-1,conk_MA2(:,2)*100,m2_MA{:})

setupLegendDC(ax)

saveas(fig{1},'plots/MAr2_DCshort.eps','epsc')
saveas(fig{2},'plots/MAobj2_DCshort.eps','epsc')
saveas(fig{3},'plots/MAcon2_DCshort.eps','epsc')

% plot plant stuff
plotPlantContour
saveas(figrr,'plots/rr_Plant.eps','epsc')

% Remove contour
for i = 1:numel(pObjContour)
    delete(pObjContour{i})
end

if isvalid(l{4})
    delete(l{4})
end

% Standard MA
pMA1 = plot(rk_MA1(:,1),rk_MA1(:,2),m1_MA{:},'color',[0.4,0.4,0.4],'MarkerFaceColor',[0.4,0.4,0.4]);
pMA2 = plot(rk_MA2(:,1),rk_MA2(:,2),m1_MA{:},'color',[0.7,0.7,0.7],'MarkerFaceColor',[0.7,0.7,0.7]);
plot(rpOpt(1),rpOpt(2),'x','LineWidth',3,'MarkerSize',15,'color',brightOrange)
plot(r0(1),r0(2),'x','LineWidth',3,'MarkerSize',15,'color','k')

set(l{7},'linestyle','-','marker','o',m1_MA{2:end},'color',[0.4,0.4,0.4],'MarkerFaceColor',[0.4,0.4,0.4])
set(l{8},'linestyle','-','marker','o',m1_MA{2:end},'color',[0.7,0.7,0.7],'MarkerFaceColor',[0.7,0.7,0.7])
leg = {'$G_{1,p}=0$','$G_{2,p}=0$','Infeasible','$\emph{\textbf{r}}_p^{\ast}$','$\emph{\textbf{r}}_0$','MA ($K=0.9$)','MA ($K=0.6$)'};
legend(axrr,leg{:},'Interpreter','latex','location','southeast')
saveas(figrr,'plots/MArr_DCshort.eps','epsc')

