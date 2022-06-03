% script for producing the plots of illustrative case 1.
clearvars -except iC1f1 iC1f2 iC1f3 iC1f23

%% 1. Parameters
n_th = 20;
rng(104)

mu_th = [2,2.5];
s_th = [0.3,0.1];
th = rand(n_th,2).*s_th*2+mu_th-s_th;

%% 2. Plots
% figure 1
if exist('iC1f1','var') && isvalid(iC1f1)
    set(0,'CurrentFigure',iC1f1);
    clf
else
    iC1f1 = figure;
end

ax1 = axes(iC1f1);
fplot(@(u)(iC1con(u,mu_th)),[0,4],'Color',brightBlue,'LineWidth',2)
hold on
plot([0,4],[0,0],'k--','LineWidth',1)

ylim([-5,5])
xlim([0,4])
ylabel('Unmodified constraint, $G(u)$','Interpreter','latex')
xlabel('Input, $u$','Interpreter','latex')

fixAxis(iC1f1,ax1);

% figure 2
if exist('iC1f2','var') && isvalid(iC1f2)
    set(0,'CurrentFigure',iC1f2);
    clf
else
    iC1f2 = figure;
end

ax2 = axes;
plot(-1,-1,'Color',brightBlue,'LineWidth',2);
hold on
plot(-1,-1,'Color',mutedLightBlue,'LineWidth',2);
xlim([0,4])
for i = 1:n_th
    fplot(@(u)(iC1con(u,th(i,:))),[0,4],'Color',mutedLightBlue,'LineWidth',1)
end

copyobj(get(ax1,'Children'),ax2);
legend({'Nominal model','Set of models'},'Interpreter','latex')

ylim([-5,5])
xlim([0,4])
ylabel('Unmodified constraint, $G(u)$','Interpreter','latex')
xlabel('Input, $u$','Interpreter','latex')

fixAxis(iC1f2,ax2);
saveas(iC1f2,'plots\iC1f2.eps','epsc')

% figure 3
if exist('iC1f3','var') && isvalid(iC1f3)
    set(0,'CurrentFigure',iC1f3);
    clf
else
    iC1f3 = figure;
end

ax3 = axes(iC1f3);
plot(-1,-1,'Color',brightBlue,'LineWidth',2);
hold on
plot(-1,-1,'Color',mutedLightBlue,'LineWidth',2);

uk = 2;
gp = iC1con(uk,mu_th);
dgpdu = (iC1con(uk+0.0001,mu_th)-iC1con(uk-0.0001,mu_th))/0.0002;

xCon = @(u,th)(iC1con(u,th)+gp-iC1con(uk,th)+...
    (u-uk).*(dgpdu-((iC1con(uk+0.0001,th)-iC1con(uk-0.0001,th))/0.0002)));

plot([0,4],[0,0],'k--','LineWidth',1)
for i = 1:n_th
    fplot(@(u)(xCon(u,th(i,:))),[0,4],'Color',mutedLightBlue,'LineWidth',1)
end
fplot(@(u)xCon(u,mu_th),[0,4],'Color',brightBlue,'LineWidth',2)

ylim([-5,5])
xlim([0,4])
ylabel('Modified constraint, $G_k(u)$','Interpreter','latex')
xlabel('Input, $u$','Interpreter','latex')

%set(ax3,'Children',[ax3.Children(1:2);ax3.Children(end);ax3.Children(3:end-1)]);

legend({'Nominal model','Set of models'},'Interpreter','latex')
fixAxis(iC1f3,ax3);
saveas(iC1f3,'plots\iC1f3.eps','epsc')

% figure 2+3
if exist('iC1f23','var') && isvalid(iC1f23)
    set(0,'CurrentFigure',iC1f23);
    clf
else
    iC1f23 = figure;
end

tiledlayout(1,2);

ax23a = nexttile;
copyobj(get(ax2,'Children'),ax23a)
ylim([-5,5])
xlim([0,4])
ylabel('Unmodified constraint, $G(u)$','Interpreter','latex')
xlabel('Input, $u$','Interpreter','latex')
legend({'Nominal model','Set of models'},'Interpreter','latex')
fixAxis(iC1f23,ax23a);

ax23b = nexttile;
plot(-1,-1,'Color',brightBlue,'LineWidth',2);
hold on
plot(-1,-1,'Color',mutedLightBlue,'LineWidth',2);
copyobj(get(ax3,'Children'),ax23b)
ylim([-5,5])
xlim([0,4])
ylabel('Modified constraint, $G_k(u)$','Interpreter','latex')
xlabel('Input, $u$','Interpreter','latex')
legend({'Nominal model','Set of models'},'Interpreter','latex')
fixAxis(iC1f23,ax23b);

