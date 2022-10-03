% makePlotsMA_WO makes the plots for:       "Standard MA"
% on the case study:                        "Williams-Otto" 
% used in:                                  "Chapter 2"
%
% 1. plant/model/obj/con functions are set-up,
% 2. figures are set-up,
% 3. axes are set-up,
% 4. RTO scheme is run,
% 5. output is plotted,
% 6. plots are saved.

%% 0. General set-up
% clear all vars (other than MAfig and above), allows for figures to be
% updated rather than made new for running multiple times (QOL thing)
clearvars -except MAfig

% variables for figure screen layout
forPub = 1;         % set to 1 for publication, and 0 to for screen use
useNiceFigPos = 1;  % set to 1 for nice layout (on my laptop)

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/williamsOttoCSTR/functions/')
addpath('../')

%% 1. Set-up functions
yGuess = [0.141, 0.332, 0.023, 0.103, 0.284, 0.117];

% plant/model functions
plant = @(u)WOplantFun(u,yGuess);       % plant function
model = @(u)WOmodelFun(u);              % model function

% optimization functions
conFun = @(u,y)WOconFun(u,y);           % constraint function
objFun = @(u,y)WOobjFun(u,y);           % objective function

% set up number of iterations
kmax = 21;

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

%% 2. Set-up figures
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

% run setupFigs to set positions on screen and create axes
ax = setupFigs(MAfig,forPub,useNiceFigPos,1);

%% 3. Set-up axes
% y-axis 
yAxisLabel = {'$J$','$X_G$','$F_A$','$F_B$','$T_R$'};
yAxisVal = [100,220;0.045,0.09;3,4.5;6,11;80,105];
yAxisTicks = {100:20:220,0.05:0.01:0.09,3:0.5:4.5,6:11,80:5:105};

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

% plant lines
cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',lw};

% optimum plot
optVal = [-objOptp;conOptp(1)-con0(1);uOptp(1);uOptp(2);uOptp(3)];

% plant infeasible area
area(ax(2),[0,kmax-1],-[con0(1),con0(1)],'BaseValue',1,ma{:})

% run for each figure
for i = 1:numel(MAfig)
    % set ylabels
    ylabel(ax(i),yAxisLabel(i),'Interpreter','latex')
    ylim(ax(i),yAxisVal(i,:));
    yticks(ax(i),yAxisTicks{i});
    
    % x-axis limits
    xlim(ax(i),[0,kmax-1]);
    
    % plot these off the figure first in the order of the desired legend
    plot(ax(i),-1,-1,mp1{:});
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

%% 5. Plot
% plot lines
plot(ax(1),0:(kmax-1),-objkMA,mp1{:});
plot(ax(2),0:(kmax-1),conkMA(:,1)-con0(1),mp1{:});
plot(ax(3),0:(kmax-1),ukMA(:,1),mp1{:});
plot(ax(4),0:(kmax-1),ukMA(:,2),mp1{:});
plot(ax(5),0:(kmax-1),ukMA(:,3),mp1{:});

% make legend
leg = {'Standard MA','Plant Optimum'};
for i = 1:numel(MAfig)
    legend(ax(i),leg,'Interpreter','latex','Location','southeast')
end

% fix constraint legend
legend(ax(2),{'Infeasible','Standard MA','Plant Optimum'},'Interpreter','latex')

%% 6. Save figure
filenames = {'obj','con','u1','u2','u3'};
if forPub
    for i = 1:numel(MAfig)
        saveas(MAfig(i),['plots\MA' filenames{i} '_WO.eps'],'epsc')
        saveas(MAfig(i),['plots\MA' filenames{i} '_WO.fig'],'fig')
    end
end

