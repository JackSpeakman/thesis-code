% makerrPlotsMA_DCshort makes the r-r plots for:    "Standard MA"
% on the case study:                                "Distillation Column" 
% used in:                                          "Chapter 6"
%
% 1. figure is set-up,
% 2. axe is set-up,
% 3. plant/model/obj/con functions are set-up,
% 4. MA schemes are run,
% 5. plot,
% 6. save.

%% 0. General set-up
% clear all vars (other than MAfig and above), allows for figures to be
% updated rather than made new for running multiple times (QOL thing)
clearvars -except rrFig

% variables for figure screen layout
forPub = 0;         % set to 1 for publication, and 0 to for screen use
useNiceFigPos = 1;  % set to 1 for nice layout (on my laptop)

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/diehlDistillationColumn/functions/')
addpath('../../../caseStudies/diehlDistillationColumn/analysis/')
addpath('../')


%% 1. Set-up figure
if exist('rrFig','var') && all(isvalid(rrFig))
    % already exists -> clear
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, rrFig));
    clf(rrFig(i))
else
    % don't exist -> create
    rrFig = figure;
end

%% 2. Set-up axes
% plot plant stuff (this is a script which makes a variable "l" used later)
plotPlantContour

% Remove contour
for i = 1:numel(pObjContour)
    delete(pObjContour{i})
end

if isvalid(l{4})
    delete(l{4})
end

%% 3. Set-up functions
% get plant optimum
fprintf('\nFinding plant optimum\n');
tic
[rpOpt,upOpt,ypOpt,objpOpt,conpOpt] = getPlantOpt;
t = toc;
fprintf('Plant optimum found [%4.2f s]\n\n',t);

% model/plant functions
[model_ri,plant,objFun,conFun] = setupDCshortFunctions;

% parameters
r0 = [70.001,85];   % starting point
k_max = 11;         % number of iterations
Ns = 14;            % number of stripping trays
Nr = 7:11;          % number of rectifying trays
mu = 3;             % which rectifying tray is nominal

% replace model_ri(r,i) with model(r)
model = @(r)model_ri(r,mu);


%% 4.a Run MA (K=0.9)
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

%% 4.b Run MA (K=0.6)
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

%% 5 Plot
% lines
c_MA = [0.4,0.4,0.4];
m1_MA = {'-o','color',c_MA,'linewidth',3,'markersize',7,'MarkerFaceColor',c_MA};
m2_MA = {'--o','color',c_MA,'linewidth',3,'markersize',7,'MarkerFaceColor',c_MA};

% Standard MA
pMA1 = plot(rk_MA1(:,1),rk_MA1(:,2),m1_MA{:},'color',[0.4,0.4,0.4],'MarkerFaceColor',[0.4,0.4,0.4]);
pMA2 = plot(rk_MA2(:,1),rk_MA2(:,2),m1_MA{:},'color',[0.7,0.7,0.7],'MarkerFaceColor',[0.7,0.7,0.7]);
plot(rpOpt(1),rpOpt(2),'x','LineWidth',3,'MarkerSize',15,'color',brightOrange)
plot(r0(1),r0(2),'x','LineWidth',3,'MarkerSize',15,'color','k')

%% 6. Save
set(l{7},'linestyle','-','marker','o',m1_MA{2:end},'color',[0.4,0.4,0.4],'MarkerFaceColor',[0.4,0.4,0.4])
set(l{8},'linestyle','-','marker','o',m1_MA{2:end},'color',[0.7,0.7,0.7],'MarkerFaceColor',[0.7,0.7,0.7])
leg = {'$G_{1,p}=0$','$G_{2,p}=0$','Infeasible','$\emph{\textbf{r}}_p^{\ast}$','$\emph{\textbf{r}}_0$','MA ($K=0.9$)','MA ($K=0.6$)'};
legend(axrr,leg{:},'Interpreter','latex','location','southeast')
saveas(figrr,'plots/MArr_DCshort.eps','epsc')