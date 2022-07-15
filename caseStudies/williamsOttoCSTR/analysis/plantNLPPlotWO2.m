% Plots the characteristic functions for the plant of the 2-var WO CSTR
addpath('../WO_functions/')

%% 1. Set-up variables
% grid space
n_u = 51;
umin = [3,6];
umax = [4.5,11];

u1 = linspace(umin(1),umax(1),n_u);
u2 = linspace(umin(2),umax(2),n_u);
[uu1,uu2] = meshgrid(u1,u2);

% plant
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
uGuess = [3.6,9];
T = 90;
model = @(u)WOplantFun([u,T],yGuess,[6600,8400,11400]);
plant = @(u)WOplantFun([u,T],yGuess);

objFun = @(u,y)WOobjFun([u,T],y);
conFun = @(u,y)WOconFun2([u,T],y);

% NLP functions
obj = zeros(size(uu1));
con = zeros([size(uu1),2]);

%% 2. Run plant
for i = 1:n_u
    for ii = 1:n_u
        % run plat at each u
        u = [uu1(i,ii),uu2(i,ii)];
        y = plant(u);
        
        % calc NLP
        obj(i,ii) = objFun(u,y);
        con(i,ii,:) = conFun(u,y);
    end
end

%% 3. Get optima
fminopts = optimoptions('fmincon','Display','off','Algorithm','interior-point','MaxFunctionEvaluations',20000,'MaxIterations',10000);

uk = fmincon(@(u)objFun(u,model(u)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,model(u)),[]),fminopts);
uk = [3.45,10.5];

up = fmincon(@(u)objFun(u,plant(u)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,plant(u)),[]),fminopts);

%% 4. Plot
close all

% create figure
fig = figure;
ax = axes;
hold on
xlim([umin(1),umax(1)])
ylim([umin(2),umax(2)])
box on
ylabel('$u_2$','Interpreter','latex')
xlabel('$u_1$','Interpreter','latex')
set(ax,'linewidth',1)
set(ax,'fontsize',16)

% plot contours for data
[a1,b] = contour(ax,uu1,uu2,con(:,:,1),[0,0]);
delete(b)
[a2,b] = contour(ax,uu1,uu2,con(:,:,2),[0,0]);
delete(b)

% plot patches
patch('XData',[umax(1),a1(1,2:end),umin(1)],'YData',[umin(2),a1(2,2:end),umin(2)],...
    'facecolor','r','facealpha',0.2,'linestyle','none')
% patch('XData',[umax(1),a2(1,2:end),umin(1)],'YData',[umin(2),a2(2,2:end),umin(2)],...
%     'facecolor','r','facealpha',0.2,'linestyle','none')

% plot objective
contour(ax,uu1,uu2,obj);

% plot optima
plot(uk(1),uk(2),'kx','markersize',12,'LineWidth',3)
plot(up(1),up(2),'bo','markersize',12,'LineWidth',3)

% legend
legend({'infeasible region','objective','$u_k$','$u_p$'},'Interpreter','latex')