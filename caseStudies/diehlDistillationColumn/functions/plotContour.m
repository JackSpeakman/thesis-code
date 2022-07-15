% plotContour
% plot the contour of the plant and the model to help determine the value
% of the objective function and constraint functions

%% 1. Set-up variables
tic

% Define efficiency
eff0 = [0.62,0.35];                 
eff = [1,1];%eff0 + -[0.1,-0.1];

% Define feed
Fin = [14/3600,0.32,71];

% Define parameters
theta_p = [8.5, 0.17, 0.155, 0.62,   0.35,   0.166, 0.5, 93900, 250, 190, 47.2];
theta_m = [8.5, 0.17, 0.155, eff(1), eff(2), 0.166, 0.51, 93900, 250/0.62, 190/0.35, 47.2];

% Set-up cost/constraint/controller functions
tFun = @(t)(logspace(2,log10(t+100),10001)-100);
tspan = 60000000;

x0 = [linspace(0.001,0.01,5),linspace(0.012,0.988,32),linspace(0.99,0.999,5)];
n0 = linspace(0.0034,0.0044,40);
V0 = linspace(0.0028,0.0032,40);
T0 = linspace(95,65,41);
Fout0 = [4/3600,10/3600]/69;

r0 = [71.7,84];
xp0 = [0,0,x0,n0,V0,T0,Fout0]';

% Define plant
K = [-0.1,-0.0002,0.1,0.0002]; % Controller gain
M = diag([ones(1,84),zeros(1,83)]);
opts = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-8);
plant = @(r,xp0)ode15s(@(t,x)plantPI(t,x,r,Fin.*[32/30,30/32,1],K,theta_p),tFun(tspan),xp0,opts);

p0 = plant(r0,xp0);
yp0 = p0.y(:,end);

% Define model
M = diag([ones(1,44),zeros(1,85)]);
opts = odeset('Mass',M,'RelTol',1e-10,'AbsTol',1e-10);
model = @(r,xm0)ode15s(@(t,x)modelPI(t,x,r,Fin,K,theta_m),tFun(tspan),xm0,opts);

% Define input space
rmin = [71,83];
rmax = [72,85];

n_v = 11;
v1 = linspace(0,1,n_v);
v2 = linspace(0,1,n_v);

[vv1,vv2] = meshgrid(v1,v2);

v2r = @(v)(v.*(rmax-rmin)+rmin);
r2v = @(r)((r-rmin)./(rmax-rmin));

fprintf('Set-up Complete       [%6.4fs]\n',toc)

%% 2. Run plant/model
tic
p = cell(n_v,n_v);
m = cell(n_v,n_v);

fprintf('Running Plant and Model\nProgress Report\n\n')
fprintf('    |||||||||||\n')

for i = 1:n_v
    fprintf('   -')
    for j = 1:n_v
        r = v2r([vv1(i,j),vv2(i,j)]);
        p{i,j} = plant(r,yp0);
        
        x0 = p{i,j}.y(3:44,end);
        V0 = [p{i,j}.y(end-82,end);p{i,j}.y(end-82:end-43,end)];
        T0 = [p{i,j}.y(end-42:end-2,end);p{i,j}.y(end-2,end)];
        Fout0 = p{i,j}.y(end-1:end,end);
        xm0 = [0;0;x0;V0;T0;Fout0];
        
        m{i,j} = model(r,xm0);
        
        fprintf('#')
    end
    fprintf('-\n')
end
fprintf('    |||||||||||\n')
fprintf('Plant/Model Complete  [%6.4fs]\n',toc)

%% 3. Analyse results
% functions
Profit = [200*32,0*60]*3.6;
cost = @(x,u)(u(2)*0.1-x(end)*Profit(2)-x(end-1)*Profit(1));
cons = @(x,u)([x(44),x(3)])*100;

contm = @(x)(x(44+[17,7]));
contp = @(x)(x(125+[28,14]));

% Pre-allocate
val.setpoint1 = zeros(n_v);
val.setpoint2 = zeros(n_v);
val.contVar1m = zeros(n_v);
val.contVar2m = zeros(n_v);
val.contVar1p = zeros(n_v);
val.contVar2p = zeros(n_v);

val.u1m = zeros(n_v);
val.u2m = zeros(n_v);
val.u1p = zeros(n_v);
val.u2p = zeros(n_v);

val.x3m = zeros(n_v);
val.x44m = zeros(n_v);
val.x3p = zeros(n_v);
val.x44p = zeros(n_v);

val.Distm = zeros(n_v);
val.Bottm = zeros(n_v);
val.Distp = zeros(n_v);
val.Bottp = zeros(n_v);

for i = 1:n_v
    for j = 1:n_v
        r = v2r([vv1(i,j),vv2(i,j)]);
        val.setpoint1(i,j) = r(1);
        val.setpoint2(i,j) = r(2);
        
        val.contVar1m(i,j) = m{i,j}.y(84+28,end);
        val.contVar2m(i,j) = m{i,j}.y(84+14,end);
        val.contVar1p(i,j) = p{i,j}.y(125+28,end);
        val.contVar2p(i,j) = p{i,j}.y(125+14,end);
        
        um = PIcontrol(r,K,m{i,j}.y(1:2,end),contm(m{i,j}.y(:,end)));
        up = PIcontrol(r,K,p{i,j}.y(1:2,end),contp(p{i,j}.y(:,end)));
        val.u1m(i,j) = um(1);
        val.u2m(i,j) = um(2);
        val.u1p(i,j) = up(1);
        val.u2p(i,j) = up(2);
        
        val.x3m(i,j) = m{i,j}.y(3,end);
        val.x44m(i,j) = m{i,j}.y(44,end);
        val.x3p(i,j) = p{i,j}.y(3,end);
        val.x44p(i,j) = p{i,j}.y(44,end);
        
        val.Distm(i,j) = m{i,j}.y(end-1,end);
        val.Bottm(i,j) = m{i,j}.y(end,end);
        val.Distp(i,j) = p{i,j}.y(end-1,end);
        val.Bottp(i,j) = p{i,j}.y(end,end);
    end
end


%% 4. Plot
Q = @(u)(u(2)-(u(2)<1)*(u(2)-1)-(u(2)>3)*(u(2)-3));
profitp = @(r,p)(p.y(end,end)*60*400*3.6+p.y(end-1,end)*32*200*3.6-Q(PIcontrol(r,K,p.y(1:2,end),contp(p.y(:,end))))*0.12);
profitm = @(r,p)(p.y(end,end)*60*400*3.6+p.y(end-1,end)*32*200*3.6-Q(PIcontrol(r,K,p.y(1:2,end),contm(p.y(:,end))))*0.12);
constraint = @(r,p)(p.y(3,end)-0.005);
opts = optimoptions('fmincon','Display','final-detailed');

vp_Opt = fmincon(@(v)(-profitp(v2r(v),plant(v2r(v),yp0))),[0.5,0.5],[],[],[],[],[0,0],[1,1],...
    @(v)deal(constraint(v2r(v),plant(v2r(v),yp0)),[]),opts);

QQ = val.u2p;
QQ(QQ<1) = 1;
QQ(QQ>3) = 3;

% plant
contour(vv1+71,vv2*2+83,val.x3p,[0,0]+0.005)
hold on
prof1 = interp2(vv1,vv2,val.Bottp*60*400*3.6+val.Distp*32*200*3.6-QQ.*0.12,vp_Opt(1),vp_Opt(2));
contour(vv1+71,vv2*2+83,val.Bottp*60*400*3.6+val.Distp*32*200*3.6-QQ.*0.12,[0,0]+prof1)
plot(vp_Opt(1)+71,vp_Opt(2)*2+83,'rx')

% model
dv = diag([0.00001,0.00001]);
dcostmdv = [0,0];
dconsmdv = [0,0];
dcostpdv = [0,0];
dconspdv = [0,0];

for i = 1:2
    v = vp_Opt + dv(i,:);
    dcostmdv(i) = (profitm(v2r(v),model(v2r(v),xm0))-profitm(v2r(vp_Opt),model(v2r(vp_Opt),xm0)))/dv(i,i);
    dconsmdv(i) = (constraint(v2r(v),model(v2r(v),xm0))-constraint(v2r(vp_Opt),model(v2r(vp_Opt),xm0)))/dv(i,i);
    
    dcostpdv(i) = (profitp(v2r(v),plant(v2r(v),yp0))-profitp(v2r(vp_Opt),plant(v2r(vp_Opt),yp0)))/dv(i,i);
    dconspdv(i) = (constraint(v2r(v),plant(v2r(v),yp0))-constraint(v2r(vp_Opt),plant(v2r(vp_Opt),yp0)))/dv(i,i);
    
end

costm_Opt = profitm(v2r(vp_Opt),model(v2r(vp_Opt),xm0));
consm_Opt = constraint(v2r(vp_Opt),model(v2r(vp_Opt),xm0));
costp_Opt = profitp(v2r(vp_Opt),plant(v2r(vp_Opt),yp0));
consp_Opt = constraint(v2r(vp_Opt),plant(v2r(vp_Opt),yp0));

QQ = val.u2m;
QQ(QQ<1) = 1;
QQ(QQ>3) = 3;

contour(vv1+71,vv2*2+83,val.x3m+(vv1-vp_Opt(1)).*(dconspdv(1)-dconsmdv(1))+(vv2-vp_Opt(2)).*(dconspdv(2)-dconsmdv(2))+consp_Opt-consm_Opt,[0,0]+0.005)
contour(vv1+71,vv2*2+83,val.Bottm*60*400*3.6+val.Distm*32*200*3.6-QQ.*0.12+(vv1-vp_Opt(1)).*(dcostpdv(1)-dcostmdv(1))+(vv2-vp_Opt(2)).*(dcostpdv(2)-dcostmdv(2))+costp_Opt-costm_Opt,[0,0]+prof1)





