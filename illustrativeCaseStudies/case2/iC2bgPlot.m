% script for producing the plots of illustrative case 2.
clearvars -except iC2f1 iC2f2 iC2f3 iC2f23

%% 1. Parameters
mu_th = [0,2];    % mean
s1 = 1;          % standard deviation 1
s2 = sqrt(0.25);          % standard deviation 2
p = 0;              % correlation

Sigma_th = [s1^2,p*s1*s2;p*s2*s1,s2^2];    % mv variance

n_th = 50;          % number of parameters

rng(100)
th = mvnrnd(mu_th,Sigma_th,n_th);           % parameters

%% 2. Plots
% figure 1
if exist('iC2f1','var') && isvalid(iC2f1)
    set(0,'CurrentFigure',iC2f1);
    clf
else
    iC2f1 = figure;
end

ax1 = axes(iC2f1);

fimplicit(@(u1,u2)(iC2con(u1,u2,mu_th(1),mu_th(2),1)),[0,1,0,1],'Color',brightBlue,'LineWidth',2)
hold on
fimplicit(@(u1,u2)(iC2con(u1,u2,mu_th(1),mu_th(2),2)),[0,1,0,1],'Color',brightOrange,'LineWidth',2)

ylim([0,1])
xlim([0,1])
% ylabel('Input 1, $u_1$','Interpreter','latex')
xlabel('Input 2, $u_2$','Interpreter','latex')

set(ax1,'Layer','top')
fixAxis(iC2f1,ax1,'box','on')

% figure 2
if exist('iC2f2','var') && isvalid(iC2f2)
    set(0,'CurrentFigure',iC2f2);
    clf
else
    iC2f2 = figure;
end
ax2 = copyobj(ax1,iC2f2);
plot(-1,-1,'Color',brightBlue,'LineWidth',2);
hold on
plot(-1,-1,'Color',brightOrange,'LineWidth',2);
plot(-1,-1,'Color',mutedLightBlue,'LineWidth',2);
plot(-1,-1,'Color',cShift(brightOrange,0.5),'LineWidth',2);
plot(-1,-1,'kx','MarkerSize',16,'LineWidth',3);

ylim([0,1])
xlim([0,1])
for i = 1:n_th
    fimplicit(@(u1,u2)(iC2con(u1,u2,th(i,1),th(i,2),1)),[0,1,0,1],'Color',mutedLightBlue,'LineWidth',1)
    fimplicit(@(u1,u2)(iC2con(u1,u2,th(i,1),th(i,2),2)),[0,1,0,1],'Color',cShift(brightOrange,0.5),'LineWidth',1)
end

%xlim([0,1])
ylabel('Input 1, $u_1$','Interpreter','latex')
xlabel('Input 2, $u_2$','Interpreter','latex')
set(ax2,'Children',[ax2.Children(end-1:end);ax2.Children(1:end-2)]);

legend({'$G_1(\emph{\textbf{u}}, \mu_\theta)=0$',...
'$G_2(\emph{\textbf{u}}, \mu_\theta)=0$',...
'$G_1(\emph{\textbf{u}}, \theta)=0$',...
'$G_2(\emph{\textbf{u}}, \theta)=0$'},...
'Location','southeast','LineWidth',1,'Interpreter','latex')

set(ax2,'Layer','top')
fixAxis(iC2f2,ax2,'box','on')

saveas(iC2f2,'plots\iC2f2.eps','epsc')

% figure 3
if exist('iC2f3','var') && isvalid(iC2f3)
    set(0,'CurrentFigure',iC2f3);
    clf
else
    iC2f3 = figure;
end

ax3 = axes(iC2f3);
plot(-1,-1,'Color',brightBlue,'LineWidth',2);
hold on
plot(-1,-1,'Color',brightOrange,'LineWidth',2);
plot(-1,-1,'Color',mutedLightBlue,'LineWidth',2);
plot(-1,-1,'Color',cShift(brightOrange,0.5),'LineWidth',2);
plot(-1,-1,'kx','MarkerSize',16,'LineWidth',3);
uk = [0.3,0.4];
gp = iC2con(uk(1),uk(2),mu_th(1),mu_th(2));

ylim([0,1])
xlim([0,1])
dgpdu = zeros(2,2);
du = diag([0.0001,0.0001]);

xdgpdu = @(u1,u2,th1,th2,i,j)((iC2con(u1+du(1,i),u2+du(2,i),th1,th2,j)-iC2con(u1-du(1,i),u2-du(2,i),th1,th2,j))/0.0002);

for i = 1:2
    for j = 1:2
        dgpdu(i,j) = xdgpdu(uk(1),uk(2),mu_th(1),mu_th(2),i,j);
    end
end

xCon = @(u1,u2,th1,th2,j)(iC2con(u1,u2,th1,th2,j)+gp(j)-iC2con(uk(1),uk(2),th1,th2,j)+...
    (u1-uk(1)).*(dgpdu(1,j)-xdgpdu(uk(1),uk(2),th1,th2,1,j))+...
    (u2-uk(2)).*(dgpdu(2,j)-xdgpdu(uk(1),uk(2),th1,th2,2,j)));

range = [0,1,0,1];

fimplicit(@(u1,u2)(xCon(u1,u2,mu_th(1),mu_th(2),1)),range,'Color',brightBlue,'LineWidth',2)
hold on
fimplicit(@(u1,u2)(xCon(u1,u2,mu_th(1),mu_th(2),2)),range,'Color',brightOrange,'LineWidth',2)

for i = 1:n_th
    fimplicit(@(u1,u2)(xCon(u1,u2,th(i,1),th(i,2),1)),range,'Color',mutedLightBlue,'LineWidth',1)
    fimplicit(@(u1,u2)(xCon(u1,u2,th(i,1),th(i,2),2)),range,'Color',cShift(brightOrange,0.5),'LineWidth',1)
end

ylabel('Input 1, $u_1$','Interpreter','latex')
xlabel('Input 2, $u_2$','Interpreter','latex')

plot(uk(1),uk(2),'kx','MarkerSize',16,'LineWidth',3)

set(ax3,'Children',[ax3.Children(end-6:end-5);ax3.Children(1:end-7);ax3.Children(end-4:end)]);

legend({'$G_{1,k}(\emph{\textbf{u}}, \mu_\theta)=0$',...
'$G_{2,k}(\emph{\textbf{u}}, \mu_\theta)=0$',...
'$G_{1,k}(\emph{\textbf{u}}, \theta)=0$',...
'$G_{2,k}(\emph{\textbf{u}}, \theta)=0$',...
'$\emph{\textbf{u}}_k$'},...
'Location','southeast','LineWidth',1,'Interpreter','latex')

set(ax3,'Layer','top')
fixAxis(iC2f3,ax3,'box','on')

saveas(iC2f3,'plots\iC2f3.eps','epsc')

% figure 2+3
if exist('iC2f23','var') && isvalid(iC2f23)
    set(0,'CurrentFigure',iC2f23);
    clf
else
    iC2f23 = figure;
end

t = tiledlayout(1,2);

ax23a = nexttile;
plot(-1,-1,'Color',brightBlue,'LineWidth',2);
hold on
plot(-1,-1,'Color',brightOrange,'LineWidth',2);
plot(-1,-1,'Color',mutedLightBlue,'LineWidth',2);
plot(-1,-1,'Color',cShift(brightOrange,0.5),'LineWidth',2);
plot(-1,-1,'kx','MarkerSize',16,'LineWidth',3);
copyobj(get(ax2,'Children'),ax23a)
ylim([0,1])
xlim([0,1])
ylabel('Input 1, $u_1$','Interpreter','latex')
xlabel('Input 2, $u_2$','Interpreter','latex')
title('(a) Unmodified constraints')
legend({'$G_{1}(\emph{\textbf{u}}, \mu_\theta)=0$',...
'$G_{2}(\emph{\textbf{u}}, \mu_\theta)=0$',...
'$G_{1}(\emph{\textbf{u}}, \theta)=0$',...
'$G_{2}(\emph{\textbf{u}}, \theta)=0$'},...
'Location','southeast','LineWidth',1,'Interpreter','latex')
set(ax23a,'Layer','top')
fixAxis(iC2f23,ax23a,'box','on');

ax23b = nexttile;
plot(-1,-1,'Color',brightBlue,'LineWidth',2);
hold on
plot(-1,-1,'Color',brightOrange,'LineWidth',2);
plot(-1,-1,'Color',mutedLightBlue,'LineWidth',2);
plot(-1,-1,'Color',cShift(brightOrange,0.5),'LineWidth',2);
plot(-1,-1,'kx','MarkerSize',16,'LineWidth',3);
copyobj(get(ax3,'Children'),ax23b)
ylim([0,1])
xlim([0,1])
ylabel('Input 1, $u_1$','Interpreter','latex')
xlabel('Input 2, $u_2$','Interpreter','latex')
title('(b) Modified constraints')
set(ax23b,'Layer','top')
fixAxis(iC2f23,ax23b,'box','on');

legend({'$G_{1,k}(\emph{\textbf{u}}, \mu_\theta)=0$',...
'$G_{2,k}(\emph{\textbf{u}}, \mu_\theta)=0$',...
'$G_{1,k}(\emph{\textbf{u}}, \theta)=0$',...
'$G_{2,k}(\emph{\textbf{u}}, \theta)=0$',...
'$\emph{\textbf{u}}$'},...
'Location','southeast','LineWidth',1,'Interpreter','latex')

saveas(iC2f23,'plots\iC2f23.eps','epsc')
