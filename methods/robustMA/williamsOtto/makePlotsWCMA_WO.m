% makePlotsWCMA_WO makes the plots for:     "Worst-case MA"
% on the case study:                        "Williams-Otto" 
% used in:                                  "Chapter 4"
%
% 1. plant/model/obj/con functions are set-up,
% 2. figures are set-up,
% 3. axes are set-up,
% 4. MA scheme (K=0.5) is run and plotted,
% 5. WCMA scheme (K=0.9) is run and plotted,
% 6. plots are saved.

%% 0. General set-up
% clear all vars (other than MAfig and above), allows for figures to be
% updated rather than made new for running multiple times (QOL thing)
clearvars -except WCMAfig

% variables for figure screen layout
forPub = 1;         % set to 1 for publication, and 0 to for screen use
useNiceFigPos = 1;  % set to 1 for nice layout (on my laptop)

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/williamsOttoCSTR/functions/')
addpath('../../modifierAdaptation/')
addpath('../')

%% 1. Set-up functions
yGuess = [0.141, 0.332, 0.023, 0.103, 0.284, 0.117];

% plant/model functions
plant = @(u)WOplantFun(u,yGuess);           % plant function
model_th = @(u,th)WOmodelFun(u,yGuess,th);  % model function

% parameters
th = [0,0; 70,160; -70,160; 70,-160; -70,-160]*3; % set of model parameters
th_nom = [0,0];                             % default nominal parameters

% nominal model
model = @(u)model_th(u,th_nom);             % nominal model function

% optimization functions
conFun = @(u,y)WOconFun2(u,y);              % constraint function
objFun = @(u,y)WOobjFun(u,y);               % objective function

% set up number of iterations
kmax = 11;

% find plant optimum (for plotting)
% fmincon set-up
uGuess = [3.9,9.3,91];
umin = [3,6,80];
umax = [4.5,11,105];
fminconopts = optimoptions('fmincon','Display','off');

% run optimization
uOptp = fmincon(@(u)objFun(u,plant(u)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,plant(u)),[]),fminconopts);

% run plant at optimum
yOptp = plant(uOptp);
objOptp = objFun(uOptp,yOptp);
conOptp = conFun(uOptp,yOptp);

%% 2. Set-up plots
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
    for i = 1:6
        WCMAfig(i) = figure;
    end
end

% run setupFigs to set positions on screen and create axes
ax = setupFigs(WCMAfig,forPub,useNiceFigPos,2);

%% 3. Set-up axes
% y-axis 
yAxisLabel = {'$J$','$X_G$','$X_B-X_E$','$F_A$','$F_B$','$T_R$'};
yAxisVal = [100,220;0.045,0.09;0.14,0.26;3,4.5;6,11;80,105];
yAxisTicks = {100:20:220,0.05:0.01:0.09,0.14:0.02:0.26,3:0.5:4.5,6:11,80:5:105};

% set up constraint zero (used for setting up plots)
con0 = conFun([],zeros(1,6));

% line styles
if forPub 
    lw = 2.5;
else
    lw = 1.5;
end

% MA line
cp1 = [0.7,0.7,0.7];
mp1 = {'^-','Color',cp1,'MarkerSize',5,'LineWidth',lw,'MarkerFaceColor',cp1};

% WCMA line
cp2 = brightRed;
mp2 = {'o-','Color',cp2,'MarkerSize',5,'LineWidth',lw,'MarkerFaceColor',cp2};

% plant lines
cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',lw};

% optimum plot
optVal = [-objOptp;conOptp(1)-con0(1);con0(2)-conOptp(2);uOptp(1);uOptp(2);uOptp(3)];

% plant infeasible area
area(ax(2),[0,kmax-1],-[con0(1),con0(1)],'BaseValue',1,ma{:})
area(ax(3),[0,kmax-1],[con0(2),con0(2)],'BaseValue',0,ma{:})

% run for each figure
for i = 1:numel(WCMAfig)
    % set ylabels
    ylabel(ax(i),yAxisLabel(i),'Interpreter','latex')
    ylim(ax(i),yAxisVal(i,:));
    yticks(ax(i),yAxisTicks{i});
    
    % x-axis limits
    xlim(ax(i),[0,kmax-1]);
    
    % plot these off the figure first in the order of the desired legend
    plot(ax(i),-1,-1,mp1{:});
    plot(ax(i),-1,-1,mp2{:});
    plot(ax(i),-1,-1,mp{:});
    area(ax(i),-1,-1,ma{:});
    
    % plot plant optimum
    plot(ax(i),[0,kmax-1],[optVal(i),optVal(i)],mp{:})
end 

drawnow

%% 4. Run MA
% set-up parameters
u0 = [];                                % initial condition
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

%% 5. Run WCMA
% run
[ukWCMA,ykWCMA,conkWCMA,objkWCMA] = runRobustMA(...
    'method','WCMA',...     % robust MA method
    'filter',1,...          % RTO input filter gain
    'kmax',kmax,...         % maximum number of iterates
    'conFun',conFun,...     % constraint function @(u,y)
    'objFun',objFun,...     % objective function @(u,y)
    'startingPoint',u0,...  % initial conditions for RTO
    'th',th,...             % model parameters
    'modelFun',model_th,... % model function @(u,th)
    'plantFun',plant);      % plant function @(u)

% plot
plot(ax(1),0:(kmax-1),-objkWCMA,mp2{:});
plot(ax(2),0:(kmax-1),conkWCMA(:,1)-con0(1),mp2{:});
plot(ax(3),0:(kmax-1),-conkWCMA(:,2)+con0(2),mp2{:});
plot(ax(4),0:(kmax-1),ukWCMA(:,1),mp2{:});
plot(ax(5),0:(kmax-1),ukWCMA(:,2),mp2{:});
plot(ax(6),0:(kmax-1),ukWCMA(:,3),mp2{:});

%%
% make legend
leg = {'Standard MA','WCMA','Plant Optimum'};
for i = 1:numel(WCMAfig)
    legend(ax(i),leg,'Interpreter','latex','Location','southeast')
end

% fix constraint legend
legend(ax(2),{'Infeasible','Standard MA','WCMA','Plant Optimum'},'Interpreter','latex','Location','southeast')
legend(ax(3),{'Infeasible','Standard MA','WCMA','Plant Optimum'},'Interpreter','latex','Location','northeast')

%% 6. Save figures
filenames = {'obj','con1','con2','u1','u2','u3'};
if forPub
    for i = 1:numel(WCMAfig)
        saveas(WCMAfig(i),['plots\WCMA' filenames{i} '_WO.eps'],'epsc')
        saveas(WCMAfig(i),['plots\WCMA' filenames{i} '_WO.fig'],'fig')
    end
end

