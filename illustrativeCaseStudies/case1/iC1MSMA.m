% runs the code relating MSMA for case 1
close all
clear

%% 1. Set-up variables

for i = 1:3
    fig(i) = figure;
    ax(i) = axes;
    hold on
    
    xlim([0,4]);
    xticks(0:4);
    xlabel(ax(i),'Input, $u$','Interpreter','latex')
    ylabel(ax(i),'Constaint, $G$','Interpreter','latex')
    
    ylim([-5,5]);
    yticks([-5,0,5]);
    
    fixAxis(fig(i),ax(i),'linewidth',2)
    set(ax(i),'Layer','Top')
end

% constraint function
G = @(u,th)((u-2).^2-3-th(:,1).*sin(th(:,2)*u));

% parameters
th = [2,2.5;
    1.7,2.4;
    1.7,2.6;
    2.3,2.4;
    2.3,2.6];

% inputs
u_range = linspace(-0.1,4.1,101);

% quad fun
quad = @(Q,u)(Q(1)*u.^2+Q(2)*u+Q(3));

% lines
m0 = {'color',brightBlue,'linewidth',2.5};
m1 = {'color',[0.7,0.7,0.7],'linewidth',1.5};

%% 2. Plot 1
plot(ax(1),u_range,G(u_range,th(1,:)),m0{:});

for i = 2:5
    plot(ax(1),u_range,G(u_range,th(i,:)),m1{:});
end

%% 3. Plot 2
% legend lines
plot(-1,-1,m0{:})
plot(-1,-1,m1{:})

% plot
plot(ax(2),u_range,G(u_range,th(1,:)),m0{:});

for i = 1:5
    Q = fmincon(@(q)sum((quad(q,u_range)-G(u_range,th(i,:))).^2),[0.25,0,-3],[],[]);
    plot(ax(2),u_range,quad(Q,u_range),m1{:})
end

% legend
legend(ax(1),{'Nominal model\quad','Set of models'},'Location','northeast','Interpreter','latex')

%% 4. Plot 3
plot(ax(3),u_range,G(u_range,th(1,:)),m0{:});
u_reg = [-0.55:0.5:2.95;1.05:0.5:4.55]';

for i = 1:size(u_reg,1)
    u_r = linspace(u_reg(i,1),u_reg(i,2),21);
    Q = fmincon(@(q)sum((quad(q,u_r)-G(u_r,th(1,:))).^2),[0.25,0,-3],[],[],[],[],[],[],@(q)deal(0.2-q(1),[]));
    plot(ax(3),u_r,quad(Q,u_r),m1{:})
end


%% 5. Save
for i = 1:3
    s = sprintf('iC1MSMA%i',i);
    saveas(fig(i),['plots\',s,'.eps'],'epsc')
end





