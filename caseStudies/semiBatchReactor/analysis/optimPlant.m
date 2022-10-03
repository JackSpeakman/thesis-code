% optimPlant
% runs an optimizer to solve the plant

%% solve
% set-up plant
plant = @(u)(batchPlant(u));

% set up NLP functions
con_y = @(u,y)([y(:,3)./(y(:,1))-0.011;y(:,2)-800;y(end,2)-150]);
obj_y = @(u,y)(-y(end,3));

con = @(u)con_y(u,plant(u));
obj = @(u)obj_y(u,plant(u));

% solver initial guess
uG = [400,400,400,154,120,400,10,20,30,40,40,40];
yG = plant(uG);

% solver options
opts = optimoptions('fmincon','Algorithm','interior-point','Display','iter','MaxFunctionEvaluations',1e5,'MaxIterations',1e4,'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-7);

% run solver
runSolver = 1;
if runSolver
    tic
    uOpt = fmincon(@(u)(obj(u)),uG,[],[],[],[],...
        [120,120,120,120,120,120,0,0,0,0,0,0],...
        [400,400,400,400,400,400,40,40,40,40,40,40],@(u)deal(con(u),[]),opts);
    toc
else
    uOpt = uG;
end

%% results
% print output
fprintf('Plant optimum (I): [%5.0f, %5.0f, %5.0f, %5.0f, %5.0f, %5.0f]\n',uOpt(1:6))
fprintf('            (F_N): [%5.1f, %5.1f, %5.1f, %5.1f, %5.1f, %5.1f]\n\n',uOpt(7:12))
fprintf('        Objective: %6.4f\n',-obj(uOpt))

% plot optimal batch
[y,t] = plant(uOpt);
close all
figure
tiledlayout(3,1,'Padding','compact','TileSpacing','compact')
hold on
title('optimal batch')

% plot 1
nexttile
hold on
plot(t,y(:,1))
plot(t,y(:,3)/0.011,'--','color',[0.7,0.7,0.7])
ylabel('$C_X$','Interpreter','latex')
xlim([0,240])

% plot 2
nexttile
hold on
ylabel('$C_N$','Interpreter','latex')
xlim([0,240])
plot(t,y(:,2))
plot(t,t*0+800)
plot(t(end),150,'rx')

% plot 3
nexttile
plot(t,y(:,3))
hold on
plot(t,y(:,1)*0.011,'--','color',[0.7,0.7,0.7])
ylabel('$C_P$','Interpreter','latex')
xlim([0,240])


