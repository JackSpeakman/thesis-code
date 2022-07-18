function [rpOpt,upOpt,ypOpt,objpOpt,conpOpt] = getPlantOpt
% Define feed
Fin = [14/3600,0.32,71];
Fdiff = [32/30,30/32,1];

% Define parameters
theta_p = [8.5, 0.17, 0.155, 0.62,   0.35,   0.166, 0.5, 93900, 250, 190, 47.2];

% Set-up cost/constraint/controller functions
Profit = [200*32,400*60]*3.6;
cost = @(x,u)((u(2)*0.12-x(end)*Profit(2)-x(end-1)*Profit(1))+3.7);
cons = @(x,u)([x(3)-0.005,0.986-x(44)])*100;

contp = @(x)(x(125+[28,14]));

tFun = @(t)(logspace(2,log10(t+100),10001)-100);
tspan = 60000000;

% Define plant
K = [-0.1,-0.0002,0.1,0.0002]; % Controller gain
M = diag([ones(1,84),zeros(1,83)]);
opts = odeset('Mass',M,'RelTol',1e-8,'AbsTol',1e-8);
plant = @(r,x0)ode15s(@(t,x)plantPI(t,x,r,Fin.*Fdiff,K,theta_p),tFun(tspan),x0,opts);

% Set-up variables for embedded functions used in the optimisation
r_last = 0;
y_last = 0;
u_last = 0;

%% 2. Find Plant Optimum
% Run plant at some point
x0 = [linspace(0.001,0.01,5),linspace(0.012,0.988,32),linspace(0.99,0.999,5)];
n0 = linspace(0.0034,0.0044,40);
V0 = linspace(0.0028,0.0032,40);
T0 = linspace(95,65,41);
Fout0 = [4/3600,10/3600]/69;

r0 = [71.7,84];
p0 = plant(r0,[0,0,x0,n0,V0,T0,Fout0]');
yp0 = p0.y(:,end);

% Find plant optimum
resetHoldVar
objF = @(r)(cOptFun(r,@(r)(plant(r,yp0)),@(r,y)(PIcontrol(r,K,y(1:2)',contp(y')).*[3600,1]),@(u,y)(cost(y,u))));
conF = @(r)(cOptFun(r,@(r)(plant(r,yp0)),@(r,y)(PIcontrol(r,K,y(1:2)',contp(y')).*[3600,1]),@(u,y)(cons(y,u))));

rpOpt = fmincon(@(r)objF(r),r0,[],[],[],[],[68,82],[74,86],@(r)(deal(conF(r),[])),opts);

p = plant(rpOpt,yp0);
ypOpt = p.y(:,end);
upOpt = PIcontrol(rpOpt,K,ypOpt(1:2)',contp(ypOpt')).*[3600,1];
objpOpt = cost(ypOpt,upOpt);
conpOpt = cons(ypOpt,upOpt);

    function [c] = cOptFun(r,odeF,uF,cF)
        % Calculates the cost/constraints but does not re-run the
        % model/plant if the last run was for r (used in fmincon)
        
        %if any(r ~= r_last) 
            % run column
            r_last = r;
            ode = odeF(r_last);
            y_last = ode.y(:,end);
            u_last = uF(r_last, y_last);
        %end
        
        c = cF(u_last, y_last);
    end

    function resetHoldVar
        % results the embedded variables
        r_last = 0;
        y_last = 0;
        u_last = 0;
    end

end