% makePlotsMA_DCshort makes the plots for:  "Worst-case MA"
% on the case study:                        "Distillation Column" 
% used in:                                  "Chapter 6"
%
% 1. plant/model/obj/con functions are set-up,
% 2. figures are set-up,
% 3. axes are set-up,
% 4. WCMA scheme is run and plotted,
% 5. plots are saved.


%% 0. General set-up
% clear all vars (other than MAfig and above), allows for figures to be
% updated rather than made new for running multiple times (QOL thing)
clearvars -except MAfig

% variables for figure screen layout
forPub = 0;         % set to 1 for publication, and 0 to for screen use
useNiceFigPos = 1;  % set to 1 for nice layout (on my laptop)

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/diehlDistillationColumn/functions/')
addpath('../../../caseStudies/diehlDistillationColumn/analysis/')
addpath('../../modifierAdaptation/')
addpath('../')


%% 1. Set-up functions
% get plant optimum
fprintf('\nFinding plant optimum\n');
tic
[rpOpt,upOpt,ypOpt,objpOpt,conpOpt] = getPlantOpt;
t = toc;
fprintf('Plant optimum found [%4.2f s]\n\n',t);

% Set-up variables to be given to runMA. The decision variables are the 
% controller setpoints (r). four functions are required (model(r), 
% plant(r), objFun(r,y), conFun(r,y)). The cost/cons functions are defined 
% in terms of the input variables (u). PIcontoller relates r/y to u. The
% following code sets up the distillation column for all models (variable
% size of rectifier, but only uses the nominal (Nr = 9) for the standard
% MA run.
[model_ri,plant,objFun,conFun] = setupDCshortFunctions;

% parameters
r0 = [70.001,85];   % starting point
k_max = 11;         % number of iterations
Ns = 14;            % number of stripping trays
Nr = 7:11;          % number of rectifying trays
mu = 3;             % which rectifying tray is nominal


%% 2. Set-up figures
% create new figures or clear figures if they already exist
% by not deleting already created figures and just clearing, the positions
% do not change which is convenient for running the code multiple times
if exist('WCMAfig','var') && all(isvalid(WCMAfig))
    % already exists -> clear
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, WCMAfig));
    for i = 1:numel(WCMAfig)
        clf(WCMAfig(i))
    end
else
    % don't exist -> create
    for i = 1:3
        WCMAfig(i) = figure;
    end
end

% run setupFigs to set positions on screen and create axes
ax = setupFigs(WCMAfig,forPub,useNiceFigPos,3);


%% 3. Set-up axes
% lines
c_MA = brightRed;
m1_MA = {'-d','color',c_MA,'linewidth',3,'markersize',7,'MarkerFaceColor',c_MA};
m2_MA = {'--d','color',c_MA,'linewidth',3,'markersize',7,'MarkerFaceColor',c_MA};

% run set-up functions
setupPlotsDC(ax,m1_MA,m2_MA)
plotPlantOpt(ax,k_max,rpOpt,objpOpt,conpOpt)


%% 4.a Run WCMA (K=0.9)
% run at K=0.9;
[rk_WCMA,yk_WCMA,conk_WCMA,objk_WCMA] = runRobustMA(...
    'method','WCMA',...     % robust MA method
    'startingPoint',r0,...  % initial conditions for RTO
    'kmax',k_max,...        % maximum number of iterates
    'filter',0.9,...        % RTO input filter gain
    'conFun',conFun,...     % constraint function @(r,y)
    'objFun',objFun,...     % objective function @(r,y)
    'modelFun',model_ri,... % model function @(r,th)
    'plantFun',plant,...    % plant function @(r)
    'th',(1:5)',...         % model parameters
    'th_nom',mu,...         % nominal model parameters
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

plot(ax{1},0:k_max-1,(rk_WCMA(:,1)-70)/4,m1_MA{:})
plot(ax{1},0:k_max-1,(rk_WCMA(:,2)-82)/4,m2_MA{:})
plot(ax{2},0:k_max-1,-objk_WCMA,m1_MA{:})
plot(ax{3},0:k_max-1,conk_WCMA(:,1)*100,m1_MA{:})
plot(ax{3},0:k_max-1,conk_WCMA(:,2)*100,m2_MA{:})

setupLegendDC(ax)

saveas(fig{1},'plots/MAr1_DCshort.eps','epsc')
saveas(fig{2},'plots/MAobj1_DCshort.eps','epsc')
saveas(fig{3},'plots/MAcon1_DCshort.eps','epsc')

