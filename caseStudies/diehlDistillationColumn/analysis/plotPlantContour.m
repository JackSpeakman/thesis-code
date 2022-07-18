% plotPlantContour
%
% Runs the plant across the whole setpoint space and plots the fesaible
% region/contour.

%% 1. Set up variables
% grid space
n_u = 21;           % number of grid points

[vv1,vv2] = ndgrid(linspace(0,0.75,n_u),linspace(0.25,1,n_u));
TT28 = vv1*4+70;
TT14 = vv2*4+82;

% Define feed
Fin = [14/3600,0.32,71];
Fdiff = [32/30,30/32,1];

% Define parameters
theta_p = [8.5, 0.17, 0.155, 0.62,   0.35,   0.166, 0.5, 93900, 250, 190, 47.2];

% Set-up cost/constraint/controller functions
Profit = [200*32,400*60]*3.6;
cost = @(x,u)((u(2)*0.12-x(end)*Profit(2)-x(end-1)*Profit(1))+3.7);
cons = @(x,u)([x(3)-0.005,0.986-x(44)])*100;

contp = @(x)(x(125+[28,14]));

tFun = @(t)(logspace(2,log10(t+100),10001)-100);
tspan = 60000000;

% Define plant
K = [-0.1,-0.0002,0.1,0.0002]; % Controller gain
M = diag([ones(1,84),zeros(1,83)]);
opts = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-8);
plant = @(r,x0)ode15s(@(t,x)plantPI(t,x,r,Fin.*Fdiff,K,theta_p),tFun(tspan),x0,opts);

% Run plant at some point
x0 = [linspace(0.001,0.01,5),linspace(0.012,0.988,32),linspace(0.99,0.999,5)];
n0 = linspace(0.0034,0.0044,40);
V0 = linspace(0.0028,0.0032,40);
T0 = linspace(95,65,41);
Fout0 = [4/3600,10/3600]/69;

r0 = [71.7,84];
p0 = plant(r0,[0,0,x0,n0,V0,T0,Fout0]');
yp0 = p0.y(:,end);

% Set-up variables for embedded functions used in the optimisation
r_last = 0;
y_last = 0;
u_last = 0;

%% 2. Run plant at each grid point
tic
cVal = zeros(n_u,n_u,3); % [obj,con1,con2]

fprintf('Running Plant\nProgress Report\n\n')
fprintf('   /%s\\\n',repmat('-',1,n_u))
for i = 1:n_u
    fprintf('   |')
    for ii = 1:n_u
        r = [TT28(i,ii),TT14(i,ii)];
        p = plant(r,yp0);
        y = p.y(:,end);
        u = PIcontrol(r,K,y(1:2),contp(y)).*[3600,1];
        
        cVal(i,ii,1) = cost(y,u);
        cVal(i,ii,2:3) = cons(y,u);
        fprintf('#')
    end
    fprintf('|\n')
end
fprintf('   \\%s/\n\n',repmat('-',1,n_u))
fprintf('Plant Complete  [%6.4fs]\n',toc)


%% 3. Plot contour
if exist('figrr','var')
    newFig = 0;
else
    newFig = 1;
end


if newFig || ~isvalid(figrr)
    figrr = figure;
else
    clf(figrr)
    axrr = [];
end
axrr = axes(figrr);
hold(axrr,'on')
fixAxis(figrr,axrr,'box','on')
set(axrr,'Layer','top')

xlabel('Set point 1, $T_{28,s}$','Interpreter','latex')
ylabel('Set point 2, $T_{14,s}$','Interpreter','latex')

xlim([70,73])
ylim([83,86])

cp = [1,0.9,0.9];
l{1} = plot(-1,-1,'-','color',[1,0.3,0.3],'LineWidth',2);
l{2} = plot(-1,-1,'--','color',[1,0.3,0.3],'LineWidth',2);
l{3} = patch([-1,-1,0,0],[-1,0,0,-1],cp,'linestyle','none');
[~,l{4}] = contour(vv1,vv2,vv1+vv2,'color',[0.7,0.7,0.7],'LineWidth',1);
l{5} = plot(-1,-1,'x','LineWidth',3,'MarkerSize',15,'color',brightOrange);
l{6} = plot(-1,-1,'x','LineWidth',3,'MarkerSize',15,'color','k');

l{7} = plot(-1,-1);
l{8} = plot(-1,-1);
l{9} = plot(-1,-1);

% infeasible region
[a1,b] = contour(TT28',TT14',cVal(:,:,2)',[0,0]);
delete(b)
[a2,b] = contour(TT28',TT14',cVal(:,:,3)',[0,0]);
delete(b)

patch([a2(1,2:end),73,73],[a2(2,2:end),83,86],cp,'linestyle','none')
patch([70,73,a1(1,2:end),70],[83,83,a1(2,2:end),83],cp,'linestyle','none')

plot(a1(1,2:end),a1(2,2:end),'color',[1,0.4,0.4],'LineWidth',2);
plot(a2(1,2:end),a2(2,2:end),'--','color',[1,0.4,0.4],'LineWidth',2);

% objective contour
if ~exist('objpOpt','var')
    fprintf('Finding Plant Optimum\n'), tic
    [rpOpt,upOpt,ypOpt,objpOpt,conpOpt] = getPlantOpt;
    fprintf('Plant Optimum Found  [%6.4fs]\n',toc)
end
[a,b] = contour(TT28',TT14',cVal(:,:,1)',objpOpt+[-0.018:0.006:0.012]);
delete(b)
j = 1;
for i = 1:4
    pObjContour{i} = plot(a(1,j+1:j+a(2,j)),a(2,j+1:j+a(2,j)),'color',[0.7,0.7,0.7].*cp,'LineWidth',1);
    j = j+a(2,j)+1;
end

for i = 5:2:7
    [x,y] = polyxpoly(a2(1,2:end),a2(2,2:end),a(1,j+1:j+a(2,j)),a(2,j+1:j+a(2,j)));
    l1x = a(1,j+1:j+a(2,j));
    l1y = a(2,j+1:j+a(2,j));
    pObjContour{i} = plot([l1x(l1x<x),x],[l1y(l1x<x),y],'color',[0.7,0.7,0.7],'LineWidth',1);
    pObjContour{i+1} = plot([x,l1x(l1x>x)],[y,l1y(l1x>x)],'color',[0.7,0.7,0.7].*cp,'LineWidth',1);
    
    j = j+a(2,j)+1;
end

% plant optimum/starting point
plot(rpOpt(1),rpOpt(2),'x','LineWidth',3,'MarkerSize',15,'color',brightOrange)
r0 = [70,85];
plot(r0(1),r0(2),'x','LineWidth',3,'MarkerSize',15,'color','k')

% legend
leg = {'$G_{1,p}=0$','$G_{2,p}=0$','Infeasible','Objective','$\emph{\textbf{r}}_p^{\ast}$','$\emph{\textbf{r}}_0$'};
legend(axrr,leg{:},'Interpreter','latex','location','southeast')

