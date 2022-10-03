% makePlotsAll_WO makes the plots for:      "All methods"
% on the case study:                        "Williams-Otto" 
% used in:                                  "Chapter 4"
%
% 1. plant/model/obj/con functions are set-up,
% 2. plots are set-up,
% 3. MA scheme (K=0.5) is run and plotted,
% 4. WCMA scheme (K=0.9) is run and plotted,
% 4. PMAi scheme (K=0.9) is run and plotted,
% 4. PMAj scheme (K=0.9) is run and plotted,
% 6. plots are saved.

%% 0. General set-up
% clear all vars (other than MAfig and above), allows for figures to be
% updated rather than made new for running multiple times (QOL thing)
clearvars -except Robustfig

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
if exist('Robustfig','var') && all(isvalid(Robustfig))
    % already exists -> clear
    all_figs = findobj(0, 'type', 'figure');
    delete(setdiff(all_figs, Robustfig));
    for i = 1:numel(Robustfig)
        clf(Robustfig(i))
    end
else
    % don't exist -> create
    for i = 1
        Robustfig(i) = figure;
    end
end

% set up axes for tiled layout
tL = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for i = 1:9
    ax(i) = nexttile;
    set(ax(i),'fontname','Times New Roman')
    set(ax(i),'fontsize',12)
    set(ax(i)','layer','top')
    set(ax(i)','linewidth',1)
    hold on
end

% set labels
ylabel(ax(1),'$J$','Interpreter','latex');
ylabel(ax(4),'$X_G$','Interpreter','latex');
ylabel(ax(7),'$X_B-X_E$','Interpreter','latex');

xlabel(ax(7),'Iteration, $k$','Interpreter','latex');
xlabel(ax(8),'Iteration, $k$','Interpreter','latex');
xlabel(ax(9),'Iteration, $k$','Interpreter','latex');

% set axis limits
% profit
ylim(ax(1),[100,220])
ylim(ax(2),[100,220])
ylim(ax(3),[100,220])
yticks(ax(1),100:20:220)
yticks(ax(2),100:20:220)
yticks(ax(3),100:20:220)

% con 1
ylim(ax(4),[0.05,0.09])
ylim(ax(5),[0.05,0.09])
ylim(ax(6),[0.05,0.09])
yticks(ax(4),0.05:0.01:0.09)
yticks(ax(5),0.05:0.01:0.09)
yticks(ax(6),0.05:0.01:0.09)

% con 2
ylim(ax(7),[0.14,0.26])
ylim(ax(8),[0.14,0.26])
ylim(ax(9),[0.14,0.26])
yticks(ax(7),0.14:0.02:0.26)
yticks(ax(8),0.14:0.02:0.26)
yticks(ax(9),0.14:0.02:0.26)

% set up constraint zero (used for setting up plots)
con0 = conFun([],zeros(1,6));

% linewidth
lw = 2;

% MA line
cp1 = [0.7,0.7,0.7];
mp1 = {'^-','Color',cp1,'MarkerSize',3,'LineWidth',lw,'MarkerFaceColor',cp1};

% WCMA line
cp2 = brightRed;
mp2 = {'o-','Color',cp2,'MarkerSize',4,'LineWidth',lw,'MarkerFaceColor',cp2};

% PMAi line
cp3 = brightBlue;
mp3 = {'s-','Color',cp3,'MarkerSize',4,'LineWidth',lw,'MarkerFaceColor',cp3};

% PMAj line
cp4 = brightGreen;
mp4 = {'d-','Color',cp4,'MarkerSize',4,'LineWidth',lw,'MarkerFaceColor',cp4};

% plant lines
cr = [1,0.8,0.8];
ma = {'LineStyle','none','FaceColor',cr,'ShowBaseLine',0};
mp = {'k--','LineWidth',lw};

% optimum plot
optVal = [-objOptp;conOptp(1)-con0(1);con0(2)-conOptp(2)];

% plot legend stuff
plot(ax(8),-1,-1,mp1{:});
plot(ax(8),-1,-1,mp2{:});
plot(ax(8),-1,-1,mp3{:});
plot(ax(8),-1,-1,mp4{:});
plot(ax(8),-1,-1,mp{:});
area(ax(8),-1,-1,ma{:});

for i = 1:3
    % plant infeasible area
    area(ax(3+i),[0,kmax-1],-[con0(1),con0(1)],'BaseValue',1,ma{:})
    area(ax(6+i),[0,kmax-1],[con0(2),con0(2)],'BaseValue',0,ma{:})
end

for i = 1:3
    for j = 1:3
        % plot plant optimum lines
        plot(ax(i*3+j-3),[0,kmax-1],[optVal(i),optVal(i)],mp{:})
    end
end

title(ax(1),'(a) WCMA')
title(ax(2),'(b) PMAi')
title(ax(3),'(c) PMAj')

drawnow

%% 3. Run MA
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
for i = 0:2
    plot(ax(1+i),0:(kmax-1),-objkMA,mp1{:});
    plot(ax(4+i),0:(kmax-1),conkMA(:,1)-con0(1),mp1{:});
    plot(ax(7+i),0:(kmax-1),-conkMA(:,2)+con0(2),mp1{:});
end

drawnow

%% 4. Run WCMA
% run
[ukWCMA,ykWCMA,conkWCMA,objkWCMA] = runRobustMA(...
    'method','WCMA',...     % robust MA method
    'filter',0.9,...          % RTO input filter gain
    'kmax',kmax,...         % maximum number of iterates
    'conFun',conFun,...     % constraint function @(u,y)
    'objFun',objFun,...     % objective function @(u,y)
    'startingPoint',u0,...  % initial conditions for RTO
    'th',th,...             % model parameters
    'modelFun',model_th,... % model function @(u,th)
    'plantFun',plant);      % plant function @(u)


plot(ax(1),0:(kmax-1),-objkWCMA,mp2{:});
plot(ax(4),0:(kmax-1),conkWCMA(:,1)-con0(1),mp2{:});
plot(ax(7),0:(kmax-1),-conkWCMA(:,2)+con0(2),mp2{:});

drawnow

%% 5. Run PMAi
% new set of parameters
sig = [70,160]*2;
corr = 0.85;
Sigma = sig'*sig.*([1,corr;corr,1]);

th = mvnrnd(th_nom,Sigma,50); % set of model parameters

% run RTO
[ukPMAi,ykPMAi,conkPMAi,objkPMAi] = runRobustMA(...
    'method','PMAi',...     % robust MA method
    'filter',0.9,...        % RTO input filter gain
    'kmax',kmax,...         % maximum number of iterates
    'conFun',conFun,...     % constraint function @(u,y)
    'objFun',objFun,...     % objective function @(u,y)
    'startingPoint',u0,...  % initial conditions for RTO
    'th',th,...             % model parameters
    'modelFun',model_th,...    % model function @(u,th)
    'plantFun',plant,...    % plant function @(u)
    'probChance',0.95);      % probability constraint chance


plot(ax(2),0:(kmax-1),-objkPMAi,mp3{:});
plot(ax(5),0:(kmax-1),conkPMAi(:,1)-con0(1),mp3{:});
plot(ax(8),0:(kmax-1),-conkPMAi(:,2)+con0(2),mp3{:});

drawnow

%% 6. Run PMAj
[ukPMAj,ykPMAj,conkPMAj,objkPMAj] = runRobustMA(...
    'method','PMAj',...     % robust MA method
    'filter',0.9,...        % RTO input filter gain
    'kmax',kmax,...         % maximum number of iterates
    'conFun',conFun,...     % constraint function @(r,y)
    'objFun',objFun,...     % objective function @(r,y)
    'startingPoint',u0,...  % initial conditions for RTO
    'th',th,...             % model parameters
    'modelFun',model_th,...    % model function @(u,th)
    'plantFun',plant,...    % plant function @(u)
    'probChance',0.95);      % probability constraint chance


plot(ax(3),0:(kmax-1),-objkPMAj,mp4{:});
plot(ax(6),0:(kmax-1),conkPMAj(:,1)-con0(1),mp4{:});
plot(ax(9),0:(kmax-1),-conkPMAj(:,2)+con0(2),mp4{:});

legend(ax(8),{'Standard MA$\quad$','WCMA$\quad$','PMAi$\quad$','PMAj$\quad$','Plant Optimum$\quad$','Infeasible'},'Interpreter','latex','Location','southoutside','Orientation','horizontal')

%% 7. Save
%saveas(Robustfig,'plots\robustAll_WO.eps','epsc')
%saveas(Robustfig,'plots\robustAll_WO.fig','fig')
