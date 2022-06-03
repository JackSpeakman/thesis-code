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

%% PMA plots
rng(1234)
thr1 = normrnd(nom(1),0.05,1,100);
thr2 = normrnd(nom(2),0.15,1,100);

% Zoomed modified constraint
figure
for i = 1:100
    fp = fplot(@(x)(f(x,thr1(i),thr2(i))+(x-uk)*(-1.4-df(uk,thr1(i),thr2(i)))+(-1-f(uk,thr1(i),thr2(i)))),...
        [3,4],'Color',[0.8,0.8,0.8]);
    ymin = min(ymin,interp1(fp.XData,fp.YData,3:0.01:4));
    ymax = max(ymax,interp1(fp.XData,fp.YData,3:0.01:4));
    hold on
end
fplot(@(x)(f(x,nom(1),nom(2))+(x-uk)*(-1.4-df(uk,nom(1),nom(2)))+(-1-f(uk,nom(1),nom(2)))),[3,4],'b','LineWidth',1)
plot([3,4],[0,0],'k--')

plot([3.4,3.4],[ymin(41),ymax(41)],'k')
plot([3.5,3.5],[ymin(51),ymax(51)],'k')
plot([3.6,3.6],[ymin(61),ymax(61)],'k')

% Histogram plots
figure
for i = 1:3
    % calculate histogram
    subplot(3,1,i)
    x = 3.3+i/10;
    val = xf(x,thr1,thr2);
    
    % plot hist
    h(i) = histogram(val,-5:0.5:3,'FaceColor',[0.8,0.8,0.8],'EdgeColor',[0.7,0.7,0.7]);
    xlim([-5,3])
    ylim([0,25])
    hold on
    box off
    
    % plot fitted distribution
    pd(i) = fitdist(val','normal');
    xv = -5:0.01:3;
    pdval = pdf(pd(i),xv);
    plot(xv,pdval*100*0.5);
    
    % plot filled pdf/x
    pvV = cdf(pd(i),0);
    pvX = icdf(pd(i),0.9);
    plot([pvX,pvX],[0,interp1(xv,pdval,pvX)]*100*0.5,'k')
    p = patch([xv,-5],[pdval(xv<=0),pdval(xv>0)*0,0]*100*0.5,'r');
    p.EdgeAlpha = 0;
    p.FaceAlpha = 0.2;
end

% Combined plots

figure
% legend stuff
g1 = brightGreen;
plot(-1,-1,'b','LineWidth',4,'color',brightBlue)
hold on
plot(-1,-1,'Color',[0.8,0.8,0.8],'LineWidth',4)
plot(-1,-1,'r','LineWidth',4)
area(-1,-1,'FaceColor',[1,0.8,0.8],'LineStyle','none','ShowBaseLine',0)
plot(-1,-1,'kx-','MarkerSize',20,'LineWidth',4)
plot(-1,-1,'Color',g1,'LineWidth',4)

% plot modified constraint
for i = 1:100
    fp = fplot(@(x)(xf(x,thr1(i),thr2(i))),...
        [3.3,3.8],'Color',[0.8,0.8,0.8],'LineWidth',2);
    
end
fplot(@(x)(xf(x,nom(1),nom(2))),[3.3,3.8],'b','LineWidth',4,'color',brightBlue)
plot([3.3,3.8],[0,0],'k--','LineWidth',2)

% plot([3.4,3.4],[ymin(41),ymax(41)],'k','LineWidth',2)
% plot([3.5,3.5],[ymin(51),ymax(51)],'k','LineWidth',2)
% plot([3.6,3.6],[ymin(61),ymax(61)],'k','LineWidth',2)

ylim([-5,2])
xlim([3.35,3.8])

% plot histograms/distribution
for i = 1:3
    % histo
    maxh = 0.08;
    minx = 3.3+i/10;
    maxx = 3.3+i/10+maxh;
    for j = 1:numel(h(i).Values)
        p = patch(minx+[0,h(i).Values(j)*maxh/25,h(i).Values(j)*maxh/25,0],...
            [h(i).BinEdges(j),h(i).BinEdges(j),h(i).BinEdges(j+1),h(i).BinEdges(j+1)],[0.8,0.8,0.8]);
        p.EdgeColor = [0.7,0.7,0.7];
    end
    
    % pdf
    pdval = pdf(pd(i),xv);
    plot(pdval*100*0.5*maxh/25+minx,xv,'r','LineWidth',2)
    
    % filled pdf
    pvV = cdf(pd(i),0);
    fprintf('%f\n',pvV)
    pvX = icdf(pd(i),0.9);
    plot(minx,[pvX],'kx','LineWidth',4,'MarkerSize',20)
    p = patch([pdval(xv<=0),pdval(xv>0)*0,0]*100*0.5*maxh/25+minx,[xv,-5],'r');
    p.EdgeAlpha = 0;
    p.FaceAlpha = 0.2;
    
    
end

% plot Quantile
r = linspace(uk,3.8,nr);
Quan = zeros(size(r));
zp = sqrt(2)*erfinv(2*0.9-1);
for i = 1:nr
    fth = xf(r(i),thr1,thr2);
    Quan(i) = mean(fth) + zp*sqrt(var(fth));
end

plot(r,Quan,'k','LineWidth',4)

ddQuan = zeros(size(r));
% plot overbarG
ddf = @(u)(2 + thr1.*thr1.*thr2.*sin(u(:)*thr1));
for i = 1:nr
    ddfth = ddf(r(i));
    ddQuan(i) = mean(ddfth) + zp*sqrt(var(ddfth));
end

plot(r,cumtrapz(r,cumtrapz(r,ddQuan))+(r-uk)*df(uk,nom(1),nom(2))+f(uk,nom(1),nom(2)),'Color',g1,'LineWidth',4)
    

box off
set(gcf,'Position',[-1000,00,1600,1200]);
ylabel('Constraint, $$G_k(u)$$','Interpreter','latex')
xlabel('Input, $u$','Interpreter','latex')
set(gca,'FontSize',40)
set(gca ,'Layer', 'Top')
set(gca,'LineWidth',4)
set(gca,'FontName','Times New Roman')
legend({'$$\,G_{k}(u,$$ \boldmath$\mu${}\unboldmath$${}_\theta$$)','$$G_{k}(u,$$ \boldmath$\theta$)','$$G_k(u)\sim\mathcal{N}$$','$$P(G_k(u)\leq0)$$','$$\mathcal{Q}_{\mathcal{N}}(0.9,u)$$','$$\overline{G}_k(u)$$'},'Location','southeast','Interpreter','latex')
% save
saveas(gcf,'plots/iC1PMA.eps','epsc')