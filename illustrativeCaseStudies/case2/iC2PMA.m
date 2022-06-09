% iC2PMA creates the plots of the second case study using the PMA/WCMA
% approaches
%
% does the 4 plots for standard, WCMA, PMAi and PMAj, then the scatter
% plot

close all
clear

% get distribution of u (weighted towards uc)
uc = [0.28,0.83];
ulim1 = [0,1]; %[0.24,0.45]; %
ulim2 = [0,1]; %[0.65,0.9]; %
n_u = 201;

u1 = icdf('Beta',linspace(0.01,0.99,n_u),4,4);
u1 = (u1-u1(1))/(u1(end)-u1(1))*(ulim1(2)-ulim1(1))+ulim1(1);
u2 = icdf('Beta',linspace(0.01,0.98,n_u),5,4);
u2 = (u2-u2(1))/(u2(end)-u2(1))*(ulim2(2)-ulim2(1))+ulim2(1);

[uu1,uu2] = meshgrid(u1,u2);
du = diag([0.0001,0.0001]);

n_th = 50;
mu_th = [0,2];
p_th = 0;
s1 = 1;
s2 = 0.5;
rng(100);
sigma_th = [s1^2,p_th*s1*s2;p_th*s1*s2,s2^2];
th = mvnrnd(mu_th,sigma_th,n_th);
th1 = th(:,1);
th2 = th(:,2);

u = [0.3,0.4];

fig1 = figure;
ax1 = axes;
hold on
% legend stuff
plot(-1,-1,'rx','LineWidth',3,'MarkerSize',15)
plot(-1,-1,'LineWidth',2,'Color',[0.3,0.5,1])
plot(-1,-1,'LineWidth',1,'Color',[0.8,0.9,1])
plot(-1,-1,'LineWidth',2,'Color',[1,0.3,0.3])
plot(-1,-1,'LineWidth',1,'Color',[1,0.9,0.8])

ylabel('$$u_2$$','Interpreter','latex')
xlabel('$$u_1$$','Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'LineWidth',1)
set(gca,'FontName','Times New Roman')
set(gcf,'Position',[-1000,0,800,600]);
box on

plot(u(1),u(2),'rx','LineWidth',2,'MarkerSize',10)
ylim([min(u2),max(u2)]); xlim([min(u1),max(u1)]);
set(ax1,'Layer','top')



s = 'abcd';

for i = 1:4
    fig2(i) = figure;
    set(gcf,'Position',[-1800+i*100,300,800,600]);
    ax2(i) = axes;
    plot(u(1),u(2),'kx')
    hold on
    ylim([min(u2),max(u2)]); xlim([min(u1),max(u1)]);
    set(ax2(i),'Layer','top')
    
    ylabel('Input 2, $u_2$','Interpreter','latex')
    xlabel('Input 1, $u_1$','Interpreter','latex')
    xticks([0.3,0.4])
    yticks([0.7,0.8,0.9])
    set(gca,'FontSize',20)
    set(gca,'LineWidth',1)
    set(gca,'FontName','Times New Roman')
    
    xlim([0.24,0.45])
    ylim([0.65,0.9])
end

% f
f = @(u1,u2,th1,th2)((th1*((u1).^2+(u2).^2)-5*(u1)+(u2)+0.5));
df = @(u1,u2,th1,th2,i)((f(u1+du(1,i),u2+du(2,i),th1,th2)-...
        f(u1-du(1,i),u2-du(2,i),th1,th2))./(2*du(i,i)));
    
for i = 1:2
    dfdu(i) = df(u(1),u(2),mu_th(1),mu_th(2),i);
end

% g
g = @(u1,u2,th1,th2)(th2*((u1).^2+(u2).^2)+1.2*(u1)+1.2*(u2)-3);
dg = @(u1,u2,th1,th2,i)((g(u1+du(1,i),u2+du(2,i),th1,th2)-...
    g(u1-du(1,i),u2-du(2,i),th1,th2))./(2*du(i,i)));

for i = 1:2
    dgdu(i) = dg(u(1),u(2),mu_th(1),mu_th(2),i);
end

xf = @(th1,th2)(f(uu1,uu2,th1,th2)+...
    (uu1-u(1))*(dfdu(1)-df(u(1),u(2),th1,th2,1))+...
    (uu2-u(2))*(dfdu(2)-df(u(1),u(2),th1,th2,2))+...
    f(u(1),u(2),mu_th(1),mu_th(2))-f(u(1),u(2),th1,th2));

xg = @(th1,th2)(g(uu1,uu2,th1,th2)+...
    (uu1-u(1))*(dgdu(1)-dg(u(1),u(2),th1,th2,1))+...
    (uu2-u(2))*(dgdu(2)-dg(u(1),u(2),th1,th2,2))+...
    g(u(1),u(2),mu_th(1),mu_th(2))-g(u(1),u(2),th1,th2));

ff = zeros(n_u,n_u,n_th);
gg = zeros(n_u,n_u,n_th);

l1 = cell(1,n_th+1);
l2 = cell(1,n_th+1);

for j = 1:n_th
    % f1
    ff(:,:,j) = xf(th1(j),th2(j));    
    
    % f2
    gg(:,:,j) = xg(th1(j),th2(j));
    
    % plot
    
    for i = 0:4
        if i == 0
            ax = ax1;
        else
            ax = ax2(i);
        end
        [a,b] = contour(ax,uu1,uu2,ff(:,:,j),[0,0],'LineColor',[0.8,0.9,1]);
        delete(b)
        l1{j} = a(:,2:end);
        [a,b] = contour(ax,uu1,uu2,gg(:,:,j),[0,0],'LineColor',[1,0.8,0.8]);
        delete(b)
        l2{j} = a(:,2:end);
    end
    
end

% plot nominal
xff = xf(mu_th(1),mu_th(2));
xgg = xg(mu_th(1),mu_th(2));
for i = 0:4
    if i == 0
        ax = ax1;
    else
        ax = ax2(i);
    end
    [a,b] = contour(ax,uu1,uu2,xff,[0,0],'LineColor',[0.3,0.5,1],'LineWidth',1);
    %delete(b)
    l1{n_th+1} = a(:,2:end);
    [a,b] = contour(ax,uu1,uu2,xgg,[0,0],'LineColor',[1,0.3,0.3],'LineWidth',1);
    %delete(b)
    l2{n_th+1} = a(:,2:end);
end

% get feasible region for PMAi/PMAj
p = 0.9;
zp = sqrt(2)*(erfinv(2*p-1));
opts = optimoptions('fsolve','Display','off');
ci1 = zeros(n_u,n_u);
ci2 = zeros(n_u,n_u);
cm = zeros(n_u,n_u);
cm0 = zeros(n_u,n_u);

for i = 1:n_u
    for j = 1:n_u
        fd = [permute(ff(i,j,:),[3,2,1]),permute(gg(i,j,:),[3,2,1])];
        
        % Pi
        pd = fitdist(fd(:,1),'Normal');
        ci1(i,j) = pd.mu+zp*pd.sigma;
        
        pd = fitdist(fd(:,2),'Normal');
        ci2(i,j) = pd.mu+zp*pd.sigma;
        
        % Pj
        try
            %pd = fitgmdist(fd,1);
            pd = [];
            pd.mu = mean(fd);
            pd.Sigma = cov(fd);
            a = max(pd.mu);
            estcm = a+zp*sqrt(pd.Sigma(diag(a==max(a))));
            if min(diag(pd.Sigma))<eps
                cm(i,j) = estcm;
            else
                cm(i,j) = fsolve(@(x)(mvncdf_new([-inf,-inf],[0,0],pd.mu-x,pd.Sigma)-p),estcm,opts);
            end
        catch
            cm(i,j) = max(mean(fd));
            
        end

        
    end
end

cp = [0.4,0.7,0.4];
ca = 0.3;

for i = 1:4
    if i == 1 % nom
        cmap = max(xff,xgg);
        [yy1,a] = contour(ax2(i),uu1,uu2,xff,[0,0],'LineColor',[0.4,0.7,0.4],'LineWidth',1);
        delete(a)
        [yy2,a] = contour(ax2(i),uu1,uu2,xgg,[0,0],'LineColor',[0.4,0.7,0.4],'LineWidth',1);
        delete(a)
        [xi,yi] = polyxpoly(yy1(1,2:end),yy1(2,2:end),yy2(1,2:end),yy2(2,2:end));
    elseif i == 2 % WC
        cmap = max(max(ff,[],3),max(gg,[],3));
        [yy1,a] = contour(ax2(i),uu1,uu2,max(ff,[],3),[0,0],'LineColor',[0.4,0.7,0.4],'LineWidth',1);
        delete(a)
        [yy2,a] = contour(ax2(i),uu1,uu2,max(gg,[],3),[0,0],'LineColor',[0.4,0.7,0.4],'LineWidth',1);
        delete(a)
        [xi,yi] = polyxpoly(yy1(1,2:end),yy1(2,2:end),yy2(1,2:end),yy2(2,2:end));
    elseif i == 3 % Pi
        cmap = max(ci1,ci2);
        [yy1,a] = contour(ax2(i),uu1,uu2,ci1,[0,0],'LineColor',[0.4,0.7,0.4],'LineWidth',1);
        delete(a)
        [yy2,a] = contour(ax2(i),uu1,uu2,ci2,[0,0],'LineColor',[0.4,0.7,0.4],'LineWidth',1);
        delete(a)
        [xi,yi] = polyxpoly(yy1(1,2:end),yy1(2,2:end),yy2(1,2:end),yy2(2,2:end));
    else % Pj
        cmap = cm;
        xi = [];
        yi = [];
    end
    
    [yy,a] = contour(ax2(i),uu1,uu2,cmap,[0,0],'LineColor',[0.4,0.7,0.4],'LineWidth',1);
    delete(a)
    ym(i) = max(yy(2,2:end));
    xm(i) = yy(1,yy(2,2:end)==ym(i));
    
    yy = sortrows([yy(:,2:end),[xi;yi]]',1)';
    if i == 2 || i == 3
        a = find(yy(1,:)==xi);
        yy(:,[a-1,a+1]) = [];
    end
    
    % deal with comeback
    area(ax2(i),yy(1,:),yy(2,:),'FaceColor',1-(1-cp)*ca,'LineStyle','none');
    
    yy_all{i} = yy;
end

% plot the lines
c1 = [0.8000    0.8894    0.9067];
c2 = [0.9608    0.8706    0.8149];

for i = 0:4
    if i == 0
        for j = 1:n_th
            plot(ax1,l1{j}(1,:),l1{j}(2,:),'Color',c1,'LineWidth',1)
            plot(ax1,l2{j}(1,:),l2{j}(2,:),'Color',c2,'LineWidth',1)
        end
        % plot nominal
        plot(ax1,l1{n_th+1}(1,:),l1{n_th+1}(2,:),'Color',[0    0.4471    0.5333],'LineWidth',2);
        plot(ax1,l2{n_th+1}(1,:),l2{n_th+1}(2,:),'Color',[0.8039    0.3529    0.0745],'LineWidth',2);
    else
        for j = 1:n_th
            % find crossing point 1
            [xi,yi] = polyxpoly(l1{j}(1,:),l1{j}(2,:),yy_all{i}(1,:),yy_all{i}(2,:));
            if numel(xi)~=1
                xi = 0;
                yi = 0;
            end
            
            % plot x -> xi
            plot(ax2(i),[l1{j}(1,l1{j}(1,:)<xi),xi],[l1{j}(2,l1{j}(1,:)<xi),yi],'Color',c1-(c1-cp)*ca,'LineWidth',1)
            % plot xi -> x
            plot(ax2(i),[xi,l1{j}(1,l1{j}(1,:)>xi)],[yi,l1{j}(2,l1{j}(1,:)>xi)],'Color',c1,'LineWidth',1)
            
            % find crossing point 1
            [xi,yi] = polyxpoly(l2{j}(1,:),l2{j}(2,:),yy_all{i}(1,:),yy_all{i}(2,:));
            if numel(xi)~=1
                xi = 1;
                yi = 1;
            end
            
            % plot x -> xi
            plot(ax2(i),[xi,l2{j}(1,l2{j}(1,:)>xi)],[yi,l2{j}(2,l2{j}(1,:)>xi)],'Color',c2-(c2-cp)*ca,'LineWidth',1)
            % plot xi -> x
            plot(ax2(i),[l2{j}(1,l2{j}(1,:)<xi),xi],[l2{j}(2,l2{j}(1,:)<xi),yi],'Color',c2,'LineWidth',1)
            
        end
        % plot nominal
        plot(ax2(i),l1{n_th+1}(1,:),l1{n_th+1}(2,:),'Color',[0    0.4471    0.5333],'LineWidth',2);
        plot(ax2(i),l2{n_th+1}(1,:),l2{n_th+1}(2,:),'Color',[0.8039    0.3529    0.0745],'LineWidth',2);
        
        % plot boundary
        plot(ax2(i),yy_all{i}(1,:),yy_all{i}(2,:),'Color',[0.4,0.7,0.4],'LineWidth',2);
    end
end
        
legend(ax1,{'$$\emph{\textbf{u}}_k$$','$$G_{1,k}(\emph{\textbf{u}},$$ \boldmath$\mu${}\unboldmath$${}_\theta$$)','$$G_{1,k}(\emph{\textbf{u}},$$ \boldmath$\theta$)','$$G_{2,k}(\emph{\textbf{u}},$$ \boldmath$\mu${}\unboldmath$${}_\theta$$)','$$G_{2,k}(\emph{\textbf{u}},$$ \boldmath$\theta$)'},'FontSize',18,'Interpreter','latex')


for i = 1:4
    saveas(fig2(i),sprintf('plots/iC2robustMA%i.eps',i),'epsc')
end

x = [0.312,0.814];
j = find((diag(uu1)-x(1))>0,1);
i = find((diag(uu2)-x(2))>0,1);
figure
% legend stuff
ax3 = axes;
hold on
plot(-1,-1,'bo','LineWidth',2,'MarkerSize',10,'color',brightBlue)
contour(-uu1-1,-uu2-1,uu1,'LineColor',[0.9,0.9,0.9],'LineWidth',1)
plot(-1,-1,'k-','LineWidth',1)

fd = [permute(ff(i,j,:),[3,2,1]),permute(gg(i,j,:),[3,2,1])];
scatter(fd(:,1),fd(:,2),'bo','LineWidth',1,'markeredgecolor',brightBlue)
hold on
plot([0,0],[-1,1],'k','LineWidth',1)
plot([-1,1],[0,0],'k','LineWidth',1)

pd = [];
pd.mu = mean(fd);
pd.Sigma = cov(fd);

n_g = 101;
g1 = linspace(pd.mu(1)-sqrt(pd.Sigma(1,1))*3,pd.mu(1)+sqrt(pd.Sigma(1,1))*3,n_g);
g2 = linspace(pd.mu(2)-sqrt(pd.Sigma(2,2))*3,pd.mu(2)+sqrt(pd.Sigma(2,2))*3,n_g);

[gg1,gg2] = meshgrid(g1,g2);

pp = mvnpdf([gg1(:),gg2(:)],pd.mu,pd.Sigma);
pp =reshape(pp,n_g,n_g);

pmax = mvnpdf(pd.mu,pd.mu,pd.Sigma);

contour(gg1,gg2,pp,pmax*(0.1:0.1:0.9),'LineColor',[0.8,0.8,0.8],'LineWidth',1)

pdi1 = fitdist(fd(:,1),'Normal');
pi1 = (erf((-pdi1.mu/pdi1.sigma)/sqrt(2))+1)/2;
pdi2 = fitdist(fd(:,2),'Normal');
pi2 = (erf((-pdi2.mu/pdi2.sigma)/sqrt(2))+1)/2;

pdj.mu = mean(fd);
pdj.Sigma = cov(fd);
pj = mvncdf_new([-inf,-inf],[0,0],pdj.mu,pdj.Sigma);



ylabel('$$G_{2,k}(\emph{\textbf{u}}_0,$$ \boldmath$\theta${}\unboldmath$$)$$','Interpreter','latex')
xlabel('$$G_{1,k}(\emph{\textbf{u}}_0,$$ \boldmath$\theta${}\unboldmath$$)$$','Interpreter','latex')
xlim([-0.75,0.3])
ylim([-0.4,0.2])
set(gca,'FontSize',24)
set(gca,'LineWidth',1)
set(gca,'FontName','Times New Roman')
set(gcf,'Position',[-1000,0,800,600]);
box on

legend({'Model Instance','Normal Distribution Fit','Constraint Feasibility Limit'},'FontSize',18,'Interpreter','latex')

saveas(gcf,'plots/iC2PMAscatter.eps','epsc')