% Illustrative example
close all

% set up current operating point
uk = 2;
% set up model function and derivative
f = @(x,th1,th2)((x-2).^2-3-sin(x.*th1).*th2);
df = @(x,th1,th2)(2*(x-2)-th1.*th2.*cos(x.*th1));

nom = [2.5,2];
epi = 0.1;
xf = @(x,th1,th2)(f(x,th1,th2)+(x-uk)*(df(uk,nom(1),nom(2))-df(uk,th1,th2))+(f(uk,nom(1),nom(2))-f(uk,th1,th2)));
q = @(x)(epi*(x-uk).^2+(x-uk)*df(uk,nom(1),nom(2))+f(uk,nom(1),nom(2)));

% set up model parameters
rng(123)
% thr1 = normrnd(nom(1),0.05,1,100);
% thr2 = normrnd(nom(2),0.05,1,100);
[thr1,thr2] = meshgrid(linspace(nom(1)-0.1,nom(1)+0.1,2),linspace(nom(2)-0.3,nom(2)+0.3,2));
thr1 = thr1(:)';
thr2 = thr2(:)';

% allocate ymax
ymax = zeros(1,401)-5;

% set up grid
epi = 0.1;
nr = 101;
r = linspace(0,4,nr);
ymin = zeros(size(r));
ymax = zeros(size(r));

%% WCMA plots
% run modified constraint
figure
ax1 = axes;
hold on
opts = optimoptions('fmincon','Display','off');
for  i = 1:nr 
    [~,ymin(i)] = fmincon(@(th)xf(r(i),th(1),th(2)),nom,[],[],[],[],...
        [min(thr1),min(thr2)],[max(thr1),max(thr2)],[],opts);
    [~,ymax(i)] = fmincon(@(th)-xf(r(i),th(1),th(2)),nom,[],[],[],[],...
        [min(thr1),min(thr2)],[max(thr1),max(thr2)],[],opts);
end
ymax = -ymax;

% for i = 1:numel(thr1)
%     fp = fplot(@(x)(f(x,thr1(i),thr2(i))+(x-uk)*(-1.4-df(uk,thr1(i),thr2(i)))+(-1-f(uk,thr1(i),thr2(i)))),...
%         [0,4],'Color',[0.8,0.8,0.8]);
%     ymax = max(ymax,interp1(fp.XData,fp.YData,0:0.01:4));
%     hold on
% end

% legend stuff
plot(-1,-1,'b','LineWidth',1.75,'color',brightBlue);
plot(-1,-1,'b--','LineWidth',1.75,'color',brightBlue);
area(-1,-1,'FaceColor',[0.8,0.8,0.8],'LineStyle','--','EdgeColor',[1,1,1]*0.5,'LineWidth',1)
area(-1,-1,'FaceColor',[1,0.9,0.9],'LineStyle','none')

% plot infeasible region
xlim([0,4])
ylim([-5,5])
ymax2 = max(ymax,q(r));
xmax(1) = interp1(ymax2(ceil(end/2):end),r(ceil(end/2):end),0);
xmax(2) = interp1(ymax2(1:floor(end/2)),r(1:floor(end/2)),0);
xl = xlim;
yl = ylim;
area([0,xmax(2)],[yl(2),yl(2)],'BaseValue',-5,'FaceColor',[1,0.9,0.9],'LineStyle','none');
area([xmax(1),4],[yl(2),yl(2)],'BaseValue',-5,'FaceColor',[1,0.9,0.9],'LineStyle','none');

% plot shaded wc area
A = area(r,[ymin;ymax-ymin]','LineStyle','none');
A(1).Visible = 0;
A(2).FaceColor = [0.8,0.8,0.8];
plot(r,ymin,'--','Color',[1,1,1]*0.5);
plot(r,ymax,'--','Color',[1,1,1]*0.5);

% plot NLP
fp = fplot(@(x)q(x),[0,4],'b--','LineWidth',1,'color',brightBlue);
fplot(@(x)xf(x,nom(1),nom(2)),[0,4],'b','LineWidth',1,'color',brightBlue)
plot([0,4],[0,0],'k--','LineWidth',1)


box off
set(ax1 ,'Layer', 'Top')
set(gcf,'Position',[0,0,800,600]);
ylabel('Constraint, $$G_{k}(u)$$','Interpreter','latex')
xlabel('$$u$$','Interpreter','latex')
xticks([0,1,2,3,4])
set(gca,'FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontName','Times New Roman')
legend({'$$\,G_{k}(u,$$ \boldmath$\mu${}\unboldmath$${}_\theta$$)','$$Q_{k}(u)$$','$$G_{k}(u,$$ \boldmath$\theta$)','$$u\not\in\mathcal{F}_k$$'},'FontSize',18,'Interpreter','latex')

% save
saveas(gcf,'Figure/ex1modconNEW1.fig')
saveas(gcf,'Figure/ex1modconNEW1.eps','epsc')

% feasible region plot
figure
ax2 = axes;
set(gcf,'Position',[0,0,800,600]);
hold on

% legend stuff
g1 = brightGreen;
plot(-1,-1,'-','LineWidth',1.75,'Color',brightBlue);
plot(-1,-1,'--','LineWidth',1.75,'Color',brightBlue);
plot(-1,-1,'-','LineWidth',1.75,'Color',g1);
area(-1,-1,'FaceColor',[0.8,0.8,0.8],'LineStyle','--','EdgeColor',[1,1,1]*0.5,'LineWidth',1)
area(-1,-1,'FaceColor',[0.85,0.95,0.85],'LineStyle','none')



% new filter region
ddf = @(u)reshape(max(2 + thr1.*thr1.*thr2.*sin(u(:)*thr1),[],2),size(u));
r1 = linspace(uk,0,1001);
r2 = linspace(uk,4,1001);
x1NF = interp1(cumtrapz(r1,cumtrapz(r1,ddf(r1)))+(r1-uk)*(df(uk,nom(1),nom(2)))+f(uk,nom(1),nom(2)),r1,0);
x2NF = interp1(cumtrapz(r2,cumtrapz(r2,ddf(r2)))+(r2-uk)*(df(uk,nom(1),nom(2)))+f(uk,nom(1),nom(2)),r2,0);
area([x1NF,x2NF],[yl(2),yl(2)],'BaseValue',-5,'FaceColor',[0.85,0.95,0.85],'LineStyle','none');
% DHMA region
% ddfDH = @(u)reshape(max(max(2 + thr1.*thr1.*thr2.*sin(u(:)*thr1),[],2),u(:)*0+epi),size(u));
% g2 = [0.2,1,0.2];
% plot(r1,cumtrapz(r1,cumtrapz(r1,ddfDH(r1)))+(r1-uk)*(-1.4)+(-1),'Color',g2);
% plot(r2,cumtrapz(r2,cumtrapz(r2,ddfDH(r2)))+(r2-uk)*(-1.4)+(-1),'Color',g2);
% x1DH = interp1(cumtrapz(r1,cumtrapz(r1,ddfDH(r1)))+(r1-uk)*(-1.4)+(-1),r1,0);
% x2DH = interp1(cumtrapz(r2,cumtrapz(r2,ddfDH(r2)))+(r2-uk)*(-1.4)+(-1),r2,0);
%area([x1DH,x2DH],[yl(2),yl(2)],'BaseValue',0,'FaceColor',[0.8,1,0.8],'LineStyle','none');

% wc region
A = area(r,[ymin;ymax-ymin]','LineStyle','none');
A(1).Visible = 0;
A(2).FaceColor = [0.8,0.8,0.8];
plot(r,ymin,'--','Color',[1,1,1]*0.5);
plot(r,ymax,'--','Color',[1,1,1]*0.5);

% Gbar
plot(r1,cumtrapz(r1,cumtrapz(r1,ddf(r1)))+(r1-uk)*(df(uk,nom(1),nom(2)))+f(uk,nom(1),nom(2)),'Color',g1,'LineWidth',1);
plot(r2,cumtrapz(r2,cumtrapz(r2,ddf(r2)))+(r2-uk)*(df(uk,nom(1),nom(2)))+f(uk,nom(1),nom(2)),'Color',g1,'LineWidth',1);



% plot NLP
fp = fplot(@(x)q(x),[0,4],'b--','LineWidth',1,'color',brightBlue);
fplot(@(x)xf(x,nom(1),nom(2)),[0,4],'b','LineWidth',1,'color',brightBlue)
plot([0,4],[0,0],'k--','LineWidth',1)

% fix axes
box off
set(ax2 ,'Layer', 'Top')
%set(gcf,'Position',[-800,200,600,1000]);
ylabel('Constraint, $$G_{k}(u)$$','Interpreter','latex')
xlabel('$$u$$','Interpreter','latex')
xlim([0,4])
ylim([-5,5])
xticks([0,1,2,3,4])
set(gca,'FontSize',24)
set(gca,'LineWidth',2)
set(gca,'FontName','Times New Roman')
legend({'$$\,G_{k}(u,$$ \boldmath$\mu${}\unboldmath$${}_\theta$$)','$$Q_{k}(u)$$','$$\overline{G}_{k}(u)$$','$$G_{k}(u,$$ \boldmath$\theta$)','$$\overline{G}_{k}(u)\geq 0$$'},'FontSize',18,'Interpreter','latex')


% save
saveas(gcf,'plots/iC1WCMA.fig')
saveas(gcf,'plots/iC1WCMA.eps','epsc')

