% script for producing the plots of illustrative case 3.
clearvars -except iC3f0 iC3f1 iC3f2 iC3f3 iC3f4 iC3f13

addpath('../../plotFunctions/')
addpath('../../plotFunctions/colours/')

%% 1. set up parameters
% plant
thp = [1,1.5];

% model 1
m1_muth = [1,1.8];
m1_Sigma = [0.2,0.2];

n_th = 30;          % number of parameters

for i = 103%1:200
    rng(i)
    m1th = mvnrnd(m1_muth,m1_Sigma,n_th);           % parameters
    mi(i) = max(m1th(:,2));
end
% model 2
m2_muth = [0.9,1.8];
m2_Sigma = [0.5,0.3];

rng(104)
m2th = mvnrnd(m2_muth,m2_Sigma,n_th);           % parameters

%% 1.5 plot zeta
if exist('iC3f0','var') && isvalid(iC3f0)
    set(0,'CurrentFigure',iC3f0);
    clf
else
    iC3f0 = figure;
end

ax0 = axes(iC3f0);

fplot(@(u1)(exp(u1)+exp(-u1)-1),[0,1.2],'Color',brightRed,'LineWidth',2)
hold on
fplot(@(u1)(u1.^2+1),[0,1.2],'Color',brightBlue,'LineWidth',2)
fplot(@(u1)(sqrt(3*u1.^2+1)),[0,1.2],'Color',brightOrange,'LineWidth',2)

ylim([1,3])
xlim([0,1.2])
ylabel('\zeta(u_1)')
xlabel('Input1, $u_1$','Interpreter','latex')

legend({'Plant','Model A','Model B'})
fixAxis(iC3f0,ax0)

%% 2. plot plant
if exist('iC3f1','var') && isvalid(iC3f1)
    set(0,'CurrentFigure',iC3f1);
    clf
else
    iC3f1 = figure;
end

ax1 = axes(iC3f1);
hold on
plot(-1,-1,'Color',brightRed,'LineWidth',2);
[xx,yy] = meshgrid([-1,-2],[-1,-2]);
contour(xx,yy,xx+yy,'Color',[0.7,0.7,0.7],'LineWidth',1);
plot(-1,-1,'x','MarkerSize',16,'LineWidth',3,'Color',brightRed);

ylim([1,2])
xlim([0,1.2])

% plant optimum
fminoptions = optimoptions('fmincon','Display','off');

upOpt = fmincon(@(u)(iC3obj(u(1),u(2),thp(1),thp(2),0)),[0.84,1.75],[],[],[],[],...
    [0,1.2],[1,2],@(u)deal(iC3con(u(1),u(2),thp(1),thp(2),0),[]),fminoptions);

% plot obj contour
objOpt = iC3obj(upOpt(1),upOpt(2),thp(1),thp(2),0);
objMin = iC3obj(1.2,2,thp(1),thp(2),0);
objMax = iC3obj(0,1,thp(1),thp(2),0);

objStep = (objMax-objMin)/4.5;

objCont = objOpt:-objStep:objMin;
objCont = [objCont(end:-1:2),objOpt:objStep:objMax];

for i = 1:numel(objCont)
    fimplicit(@(u1,u2)(iC3obj(u1,u2,thp(1),thp(2),0)-objCont(i)),[0,1.2,1,2],'Color',[0.7,0.7,0.7],'LineWidth',1)
    hold on
end

fimplicit(@(u1,u2)(iC3con(u1,u2,thp(1),thp(2),0)),[0,1.2,1,2],'Color',brightRed,'LineWidth',2)

plot(upOpt(1),upOpt(2),'x','MarkerSize',16,'LineWidth',3,'Color',brightRed)

xlabel('Input1, $u_1$','Interpreter','latex')
ylabel('Input2, $u_2$','Interpreter','latex')
fixAxis(iC3f1,ax1,'box','on')

l2 = legend({'$G_p(\emph{\textbf{u}})=0$',...
'$\Phi_p(\emph{\textbf{u}})$',...
'$\emph{\textbf{u}}_p^\ast$'},...
'Location','southeast','LineWidth',1,'Interpreter','latex');

saveas(iC3f1,'plots\iC3f1.eps','epsc')

%% 3. Add nominal model's
if exist('iC3f2','var') && isvalid(iC3f2)
    set(0,'CurrentFigure',iC3f2);
    clf
else
    iC3f2 = figure;
end

ax2 = copyobj(ax1,iC3f2);

% plot constraint
fimplicit(@(u1,u2)(iC3con(u1,u2,m1_muth(1),m1_muth(2),1)),[0,1.2,1,2],'Color',brightBlue,'LineWidth',2,'LineStyle','--')
fimplicit(@(u1,u2)(iC3con(u1,u2,m2_muth(1),m2_muth(2),2)),[0,1.2,1,2],'Color',brightOrange,'LineWidth',2,'LineStyle','--')

% find optima
upOpt1 = fmincon(@(u)(iC3obj(u(1),u(2),m1_muth(1),m1_muth(2),1)),[0.84,1.75],[],[],[],[],...
    [0,1.2],[1,2],@(u)deal(iC3con(u(1),u(2),m1_muth(1),m1_muth(2),1),[]),fminoptions);
plot(upOpt1(1),upOpt1(2),'x','Color',brightBlue,'MarkerSize',16,'LineWidth',3)

upOpt2 = fmincon(@(u)(iC3obj(u(1),u(2),m2_muth(1),m2_muth(2),2)),[0.84,1.75],[],[],[],[],...
    [0,1.2],[1,2],@(u)deal(iC3con(u(1),u(2),m2_muth(1),m2_muth(2),2),[]),fminoptions);
plot(upOpt2(1),upOpt2(2),'x','Color',brightOrange,'MarkerSize',16,'LineWidth',3)

xlabel('Input1, $u_1$','Interpreter','latex')
ylabel('Input2, $u_2$','Interpreter','latex')
fixAxis(iC3f2,ax2,'box','on')

%% 4. Plot modified functions
if exist('iC3f3','var') && isvalid(iC3f3)
    set(0,'CurrentFigure',iC3f3);
    clf
else
    iC3f3 = figure;
end

ax3 = axes(iC3f3);
plot(-1,-1,'Color',brightRed,'LineWidth',2);
plot(-1,-1,'--','Color',brightBlue,'LineWidth',2);
plot(-1,-1,'--','Color',brightOrange,'LineWidth',2);
plot(-1,-1,'x','MarkerSize',16,'LineWidth',3,'Color','k');

uk = [0.5,1.5];
conp = iC3con(uk(1),uk(2),thp(1),thp(2),0);
objp = iC3obj(uk(1),uk(2),thp(1),thp(2),0);

dconpdu = [0,0];
dobjpdu = [0,0];
du = diag([0.0001,0.0001]);

xdconpdu = @(u1,u2,th1,th2,i,j)((iC3con(u1+du(1,i),u2+du(2,i),th1,th2,j)-iC3con(u1-du(1,i),u2-du(2,i),th1,th2,j))/0.0002);
xdobjpdu = @(u1,u2,th1,th2,i,j)((iC3obj(u1+du(1,i),u2+du(2,i),th1,th2,j)-iC3obj(u1-du(1,i),u2-du(2,i),th1,th2,j))/0.0002);

for i = 1:2
    dconpdu(i) = xdconpdu(uk(1),uk(2),thp(1),thp(2),i,0);
    dobjpdu(i) = xdobjpdu(uk(1),uk(2),thp(1),thp(2),i,0);
end

xCon = @(u1,u2,th1,th2,j)(iC3con(u1,u2,th1,th2,j)+conp-iC3con(uk(1),uk(2),th1,th2,j)+...
    (u1-uk(1)).*(dconpdu(1)-xdconpdu(uk(1),uk(2),th1,th2,1,j))+...
    (u2-uk(2)).*(dconpdu(2)-xdconpdu(uk(1),uk(2),th1,th2,2,j)));
xObj = @(u1,u2,th1,th2,j)(iC3obj(u1,u2,th1,th2,j)+objp-iC3obj(uk(1),uk(2),th1,th2,j)+...
    (u1-uk(1)).*(dobjpdu(1)-xdobjpdu(uk(1),uk(2),th1,th2,1,j))+...
    (u2-uk(2)).*(dobjpdu(2)-xdobjpdu(uk(1),uk(2),th1,th2,2,j)));

% plot plant
fimplicit(@(u1,u2)(iC3con(u1,u2,thp(1),thp(2),0)),[0,1.2,1,2],'Color',brightRed,'LineWidth',2)
hold on

% plot constraint
fimplicit(@(u1,u2)(xCon(u1,u2,m1_muth(1),m1_muth(2),1)),[0,1.2,1,2],'Color',brightBlue,'LineWidth',2,'LineStyle','--')
fimplicit(@(u1,u2)(xCon(u1,u2,m2_muth(1),m2_muth(2),2)),[0,1.2,1,2],'Color',brightOrange,'LineWidth',2,'LineStyle','--')

plot(uk(1),uk(2),'kx','MarkerSize',16,'LineWidth',3)

xlabel('Input 1, $u_1$','Interpreter','latex')
ylabel('Input 2, $u_2$','Interpreter','latex')
fixAxis(iC3f3,ax3,'box','on')

l2 = legend({'$G_p(\emph{\textbf{u}})=0$',...
'$G_k(\emph{\textbf{u}},$ \boldmath$\theta$)=0 (Structure A)',...
'$G_k(\emph{\textbf{u}},$ \boldmath$\theta$)=0 (Structure B)',...
'$\emph{\textbf{u}}_k$'},...
'Location','southeast','LineWidth',1,'Interpreter','latex');

saveas(iC3f3,'plots\iC3f3.eps','epsc')

%% 4.5 Combined figure
if exist('iC3f13','var') && isvalid(iC3f13)
    set(0,'CurrentFigure',iC3f13);
    clf
else
    iC3f13 = figure;
end

t = tiledlayout(1,2);

ax13a = nexttile;
hold on
plot(-1,-1,'Color',brightRed,'LineWidth',2);
plot(-1,-1,'Color',[0.7,0.7,0.7],'LineWidth',2);
plot(-1,-1,'x','MarkerSize',16,'LineWidth',3,'Color',brightRed);

plot(-1,-1,'--','Color',brightBlue,'LineWidth',2);
plot(-1,-1,'--','Color',brightOrange,'LineWidth',2);
plot(-1,-1,'kx','MarkerSize',16,'LineWidth',3);

copyobj(get(ax1,'Children'),ax13a)
ylim([1,2])
xlim([0,1.2])
ylabel('Input2, $u_2$','Interpreter','latex')
xlabel('Input1, $u_1$','Interpreter','latex')
title('(a) Plant')
l1 = legend({'$G_p(\emph{\textbf{u}})=0$',...
'$\Phi_p(\emph{\textbf{u}})$',...
'$\emph{\textbf{u}}_p^\ast$',...
'$G(\emph{\textbf{u}},$ \boldmath$\theta$ )=0 (structure A)',...
'$G(\emph{\textbf{u}}, {\bf\theta})=0$ (structure B)',...
'$\emph{\textbf{u}}_k$'},...
'Location','southoutside','LineWidth',1,'Orientation','horizontal','Interpreter','latex');
l1.Visible = 0;
fixAxis(iC3f13,ax13a,'box','on');

ax13b = nexttile;
hold on
plot(-1,-1,'Color',brightRed,'LineWidth',2);
plot(-1,-1,'Color',[0.7,0.7,0.7],'LineWidth',2);
plot(-1,-1,'x','MarkerSize',16,'LineWidth',3,'Color',brightRed);

plot(-1,-1,'--','Color',brightBlue,'LineWidth',2);
plot(-1,-1,'--','Color',brightOrange,'LineWidth',2);
plot(-1,-1,'kx','MarkerSize',16,'LineWidth',3);

copyobj(get(ax3,'Children'),ax13b)
ylim([0,1])
xlim([0,1])
ylabel('Input 2, $u_2$','Interpreter','latex')
xlabel('Input 1, $u_1$','Interpreter','latex')
title('(a) Unmodified constraints')
legend({'G_p({\bfu})=0',...
'\Phi_p({\bfu})',...
'{\bfu}_p^\ast',...
'G_k({\bfu}, {\bf\theta})=0 (structure 1)',...
'G_k({\bfu}, {\bf\theta})=0 (structure 2)',...
'{\bfu}_k'},...
'Location','southeast','LineWidth',1)

ylim([1,2])
xlim([0,1.2])
ylabel('Input2, $u_2$','Interpreter','latex')
xlabel('Input1, $u_1$','Interpreter','latex')
title('(b) Modified constraints')
fixAxis(iC3f13,ax13b,'box','on');

l2 = legend({'G_p({\bfu})=0   ',...
'\Phi_p({\bfu})   ',...
'{\bfu}_p^\ast   ',...
'G_k({\bfu}, {\bf\theta})=0 (structure A)   ',...
'G_k({\bfu}, {\bf\theta})=0 (structure B)   ',...
'{\bfu}_k  '},...
'Location','southoutside','LineWidth',1,'Orientation','horizontal');

p = l2.Position;
l2.Position = [(1-p(3))/2,p(2:4)];


%% 5. Extra bit
if exist('iC3f4','var') && isvalid(iC3f4)
    set(0,'CurrentFigure',iC3f4);
    clf
else
    iC3f4 = figure;
end

ax4 = axes;
hold on
uSet = zeros(n_th,2,2);

m = {'s','color',cShift(brightBlue,0.1),'MarkerSize',10,'LineWidth',2;...
    'o','color',cShift(brightOrange,0.1),'MarkerSize',9,'LineWidth',2};

plot(-1,-1,m{1,:},'MarkerSize',11)
plot(-1,-1,m{2,:})
%plot(-1,-1,'x','LineWidth',3,'color',brightRed,'MarkerSize',15);
plot(-1,-1,'-','LineWidth',2,'color',brightRed);
% plot(-1,-1,'--','LineWidth',2,'color',brightBlue);
% plot(-1,-1,'--','LineWidth',2,'color',brightOrange);
contour(xx,yy,xx+yy,'LineWidth',1.5,'color',[0.7,0.7,0.7]);

copyobj(get(ax3,'Children'),ax4)

for i = 1:n_th
    % find optima
    upOpt1 = fmincon(@(u)(xObj(u(1),u(2),m1th(i,1),m1th(i,2),1)),[0.84,1.75],[],[],[],[],...
        [0,1],[1.2,2],@(u)deal(xCon(u(1),u(2),m1th(i,1),m1th(i,2),1),[]),fminoptions);
    plot(upOpt1(1),upOpt1(2),m{1,:})
    
    upOpt2 = fmincon(@(u)(xObj(u(1),u(2),m2th(i,1),m2th(i,2),2)),[0.84,1.75],[],[],[],[],...
        [0,1],[1.2,2],@(u)deal(xCon(u(1),u(2),m2th(i,1),m2th(i,2),2),[]),fminoptions);
    plot(upOpt2(1),upOpt2(2),m{2,:})
    
    % save optima
    uSet(i,:,1) = upOpt1;
    uSet(i,:,2) = upOpt2;
    
end

% fit dist 1
mu_dist(:,1) = mean(uSet(:,:,1));
mu_dist(:,2) = mean(uSet(:,:,2));
cov_dist(:,:,1) = cov(uSet(:,:,1));
cov_dist(:,:,2) = cov(uSet(:,:,2));

gm_dist = gmdistribution(mu_dist',cov_dist);

nd = 201;
[uu1d,uu2d] = meshgrid(linspace(0.7,1.1,nd),linspace(1.67,1.85,nd));
p = zeros(size(uu1d));
for i = 1:nd
    for j = 1:nd
        p(i,j) = pdf(gm_dist,[uu1d(i,j),uu2d(i,j)]);
    end
end

maxP = max(p,[],'all');

for i = 1:7
    % get range
    pVal = i*(maxP)/8;
    pA = p>(pVal);
    
    u1min = min(uu1d(pA))-uu1d(3,3)+uu1d(1,1);
    u1max = max(uu1d(pA))+uu1d(3,3)-uu1d(1,1);
    u2min = min(uu2d(pA))-uu2d(3,3)+uu2d(1,1);
    u2max = max(uu2d(pA))+uu2d(3,3)-uu2d(1,1);
    
    fimplicit(@(u1,u2)(reshape(pdf(gm_dist,[u1(:),u2(:)]),size(u1))-pVal),[u1min,u1max,u2min,u2max],'Color',[0.7,0.7,0.7],'LineWidth',0.5);
end

fixAxis(iC3f4,ax4,'box','on')
set(ax4,'layer','top')
xlabel('Input 1, $u_1$','Interpreter','latex')
ylabel('Input 2, $u_2$','Interpreter','latex')
xlim([0.7,1.05])
ylim([1.68,1.84])

% get all solutions
xsol = zeros(5,2);

modmag = zeros(n_th*2,5);

   
modmag(:,1) = objp-[iC3obj(uk(1),uk(2),m1th(:,1),m1th(:,2),1);iC3obj(uk(1),uk(2),m2th(:,1),m2th(:,2),2)];
modmag(:,2) = conp-[iC3con(uk(1),uk(2),m1th(:,1),m1th(:,2),1);iC3con(uk(1),uk(2),m2th(:,1),m2th(:,2),2)];
modmag(:,3) = dobjpdu(1)-[xdobjpdu(uk(1),uk(2),m1th(:,1),m1th(:,2),1,1);xdobjpdu(uk(1),uk(2),m2th(:,1),m2th(:,2),1,2)];
modmag(:,4) = dconpdu(1)-[xdconpdu(uk(1),uk(2),m1th(:,1),m1th(:,2),1,1);xdconpdu(uk(1),uk(2),m2th(:,1),m2th(:,2),1,2)];
modmag(:,5) = dconpdu(2)-[xdconpdu(uk(1),uk(2),m1th(:,1),m1th(:,2),2,1);xdconpdu(uk(1),uk(2),m2th(:,1),m2th(:,2),2,2)];

%     xCon = @(u1,u2,th1,th2,j)(iC3con(u1,u2,th1,th2,j)+conp-iC3con(uk(1),uk(2),th1,th2,j)+...
%     (u1-uk(1)).*(-)+...
%     (u2-uk(2)).*(-xdconpdu(uk(1),uk(2),th1,th2,2,j)));
% xObj = @(u1,u2,th1,th2,j)(iC3obj(u1,u2,th1,th2,j)++...
%     (u1-uk(1)).*)+...
%     (u2-uk(2)).*(dobjpdu(2)-xdobjpdu(uk(1),uk(2),th1,th2,2,j)));
% 
% [[abs((p0(1)-m0(:,1,1)));abs((p0(1)-m0(:,1,2)))],...
%     [abs((p0(2)-m0(:,2,1)));abs((p0(2)-m0(:,2,2)))],...
%     [abs(dobjp(1)-dobj(:,1,1));abs(dobjp(1)-dobj(:,1,2))],...
%     [abs(dconp-dcon(:,:,1));abs(dconp-dcon(:,:,2))]];
[~,modmagi] = sort(abs(modmag));
summodmag = sum(modmagi,2);
[~,maxmodi] = max(summodmag);
xsol(1,:) = uSet(mod(maxmodi-1,n_th)+1,:,ceil(maxmodi/n_th));

xsol(2,:) = mean(mean(uSet,3),1);

a = uSet(:,1,:);
b = uSet(:,2,:);
c = (uk(1)-a(:)).^2+(uk(2)-b(:)).^2;
xsol(3,:) = [a(c==min(c)),b(c==min(c))];

opts = optimoptions('fmincon','Display','off');
xsol(4,:) = fmincon(@(u)(-pdf(gm_dist,u)),mean(mu_dist,2)',[],[],[],[],[0.7,1.67],[1.1,1.85],[],opts);

xsol(5,:) = fmincon(@(y)(sum(sqrt((y(1)-a(:)).^2+(y(2)-b(:)).^2))),median([a(:),b(:)]),[],[],[],[],[0.7,1.67],[1.1,1.85],[],opts);

plot(ax4,xsol(:,1),xsol(:,2),'kx','markersize',12,'linewidth',2)

p = get(gca,'Position');

annotation('arrow',p(1)+p(3)*((xsol(1,1)-0.7)/0.35+[0.23,0.0046]),p(2)+p(4)*((xsol(1,2)-1.68)/0.16-[0.18,0.0036]))
annotation('arrow',p(1)+p(3)*((xsol(2,1)-0.7)/0.35+[0.23,0.004]),p(2)+p(4)*((xsol(2,2)-1.68)/0.16+[0.05,0.0006]))
annotation('arrow',p(1)+p(3)*((xsol(3,1)-0.7)/0.35+[0.05,0.001]),p(2)+p(4)*((xsol(3,2)-1.68)/0.16+[0.35,0.007]))
annotation('arrow',p(1)+p(3)*((xsol(4,1)-0.7)/0.35+[0.2,0.004]),p(2)+p(4)*((xsol(4,2)-1.68)/0.16-[0.04,0.0002]))
annotation('arrow',p(1)+p(3)*((xsol(5,1)-0.7)/0.35+[0.27,0.0054]),p(2)+p(4)*((xsol(5,2)-1.68)/0.16+[0.08,0.0016]))

annotation('textbox',[p(1)+p(3)*((xsol(1,1)-0.7)/0.35+0.23),p(2)+p(4)*((xsol(1,2)-1.68)/0.16-0.18),0,0.03],'String','Minimum modifier','Interpreter','latex','LineStyle','none','FontSize',18)
annotation('textbox',[p(1)+p(3)*((xsol(2,1)-0.7)/0.35+0.23),p(2)+p(4)*((xsol(2,2)-1.68)/0.16+0.05),0,0.03],'String','Mean','Interpreter','latex','LineStyle','none','FontSize',18)
annotation('textbox',[p(1)+p(3)*((xsol(3,1)-0.7)/0.35+0.06)-0.06,p(2)+p(4)*((xsol(3,2)-1.68)/0.16+0.42),0,0.0],'String','Closest','Interpreter','latex','LineStyle','none','FontSize',18)
annotation('textbox',[p(1)+p(3)*((xsol(4,1)-0.7)/0.35+0.2),p(2)+p(4)*((xsol(4,2)-1.68)/0.16-0.04),0,0.03],'String','Most likely','Interpreter','latex','LineStyle','none','FontSize',18)
annotation('textbox',[p(1)+p(3)*((xsol(5,1)-0.7)/0.35+0.27),p(2)+p(4)*((xsol(5,2)-1.68)/0.16+0.08),0,0.03],'String','Median','Interpreter','latex','LineStyle','none','FontSize',18)


leg = {'Model A','Model B','$G_p = 0$','Mixed Gaussian fit'};
legend(leg,'Location','northeast','Interpreter','latex')

saveas(iC3f4,'MSMAex2.eps','epsc')








