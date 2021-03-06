% iC1DHMA makes the plots for DHMA using case study 1

%% Set up variables
% input space
rlim = [-2,2];
n_r = 201;
r = linspace(rlim(1),rlim(2),n_r);

% theta space
n_th = 20;
rng(104)

mu_th = [2,2.5];
s_th = [0.3,0.1];
th = rand(n_th,2).*s_th*2+mu_th-s_th;

% functions
Gfun = @(r,th)((r).^2-th(1)*sin(th(2)*(r+2))-1.1+th(1)*sin(2*th(2))+r*(-1.4+th(1)*th(2)*cos(2*th(2))));
Sfun = @(r,th)(2+th(1)*th(2)^2*sin(th(2)*(r+2)));

% SS min
epi = 0.1;

%% Plot distribution of S
if exist('figSS','var') && isvalid(figSS)
    set(0,'CurrentFigure',figSS);
    clf
else
    figSS = figure;
end
axSS = axes;
hold on
xlim([0,4])
ylim([-20,20])
xlabel('Input, $u$','Interpreter','latex')
ylabel('Spectral surface, $S^G$','Interpreter','latex')


plot(-1,-1,'-','color',cShift(brightBlue,0),'LineWidth',2)
plot(-1,-1,'-','color',cShift(brightBlue,0.5),'LineWidth',1)
plot(-1,-1,'--','color',cShift(brightRed,0),'LineWidth',2)

plot([0,4],[0,0],'k--','LineWidth',1)
for i = 1:n_th
    plot(r+2,Sfun(r,th(i,:)),'color',cShift(brightBlue,0.5))
end

plot(r+2,Sfun(r,mu_th),'color',cShift(brightBlue,0),'LineWidth',2)

%% plot analytical maximum
opts = optimoptions('fmincon','Display','off');
th1max = zeros(n_r,1);
th2max = zeros(n_r,1);
Smax = zeros(n_r,1);

for ii = 1:n_r
    th2max(ii) = fmincon(@(th2)2-Sfun(r(ii),[mu_th(1),th2]),mu_th(2),[],[],[],[],mu_th(2)-s_th(2),mu_th(2)+s_th(2),[],opts);
    th1max(ii) = mu_th(1)-s_th(1)*sign(2-Sfun(r(ii),[mu_th(1),th2max(ii)]));
    
    Smax(ii) = Sfun(r(ii),[th1max(ii),th2max(ii)]);
end
Smax(Smax<epi) = epi;
plot(r+2,Smax,'--','color',cShift(brightRed,0.1),'LineWidth',2)

fixAxis(figSS,axSS)
leg = {'Nominal model','Set of models','Analytical maximum'};
legend(axSS,leg,'Location','southeast','Interpreter','latex')
saveas(figSS,'plots/iC1DHMA_SS.eps','epsc')

%% plot G
if exist('figG','var') && isvalid(figG)
    set(0,'CurrentFigure',figG);
    clf
else
    figG = figure;
end

axG = axes;
hold on
ylim([-5,5])
xlim([0,4])
xlabel('Input, $u$','Interpreter','latex')
ylabel('Modified constraint, $G_k$','Interpreter','latex')
set(axG,'layer','top')
fixAxis(figG,axG)

plot(-1,-1,'-','color',cShift(brightBlue,0),'LineWidth',2)
plot(-1,-1,'-','color',cShift(brightBlue,0.5),'LineWidth',1)
plot(-1,-1,'--','color',cShift(brightRed,0),'LineWidth',2)
cal = 0.9;
cp = [1,cal,cal];
patch('XData',[-1,-1,-2,-2],'YData',[-1,-2,-2,-1],'facecolor',cp,'linestyle','none')


% find Gbar
delta = zeros(n_r,1);
delta(r>=0) = cumtrapz(r(r>=0),cumtrapz(r(r>=0),Smax(r>=0)));
Smaxi = Smax(end:-1:1);
ri = r(end:-1:1);
deltai = cumtrapz(ri(ri<=0),cumtrapz(ri(ri<=0),Smax(ri<=0)));
delta(r<=0) = deltai(end:-1:1);
Gbar = delta - 1.1 - 1.4*r';

% get feasible region
[x,y] = polyxpoly(r,Gbar,r,r*0);

patch('XData',[-2,x(1),x(1),-2]+2,'YData',[-5,-5,5,5],'facecolor',cp,'linestyle','none')
patch('XData',[2,x(2),x(2),2]+2,'YData',[-5,-5,5,5],'facecolor',cp,'linestyle','none')

% plot region 1
plot([-2,x(1)]+2,[0,0],'--','LineWidth',1,'color',(1-cal)*[1,0,0])
for i = 1:n_th
    G = Gfun(r,th(i,:));
    g1 = interp1(r,G,x(1));
    plot([r(r<x(1)),x(1)]+2,[G(r<x(1)),g1],'color',cShift(brightBlue,0.5)*cal+(1-cal)*[1,0,0])
end
G = Gfun(r,mu_th);
g1 = interp1(r,G,x(1));
plot([r(r<x(1)),x(1)]+2,[G(r<x(1)),g1],'color',cShift(brightBlue,0.1)*cal+(1-cal)*[1,0,0],'LineWidth',2)

% plot region 2
plot([x(1),x(2)]+2,[0,0],'--','LineWidth',1,'color','k')
for i = 1:n_th
    G = Gfun(r,th(i,:));
    g1 = interp1(r,G,x(1));
    g2 = interp1(r,G,x(2));
    plot([x(1),r(r>x(1) & r<x(2)),x(2)]+2,[g1,G(r>x(1) & r<x(2)),g2],'color',cShift(brightBlue,0.5))
end
G = Gfun(r,mu_th);
g1 = interp1(r,G,x(1));
g2 = interp1(r,G,x(2));
plot([x(1),r(r>x(1) & r<x(2)),x(2)]+2,[g1,G(r>x(1) & r<x(2)),g2],'color',cShift(brightBlue,0.1),'LineWidth',2)

% plot region 3
plot([x(2),2]+2,[0,0],'--','LineWidth',1,'color','k')
for i = 1:n_th
    G = Gfun(r,th(i,:));
    g2 = interp1(r,G,x(2));
    plot([x(2),r(r>x(2))]+2,[g2,G(r>x(2))],'color',cShift(brightBlue,0.5)*cal+(1-cal)*[1,0,0])
end
G = Gfun(r,mu_th);
g2 = interp1(r,G,x(2));
plot([x(2),r(r>x(2))]+2,[g2,G(r>x(2))],'color',cShift(brightBlue,0.1)*cal+(1-cal)*[1,0,0],'LineWidth',2)



plot(r+2,Gbar,'--','color',cShift(brightRed,0.1),'LineWidth',2);

leg = {'Nominal model','Set of models','DHMA','DHMA infeasible'};
legend(axG,leg,'Location','northeast','Interpreter','latex')

saveas(figG,'plots/iC1DHMA_G.eps','epsc')









