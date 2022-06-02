% calculates the lagrangian of the model at the plant optimum
clear

% add path
addpath('../WO_functions/')

%% 1. Find plant optimum
uGuess = [3.9,9.3,91];
umin = [3,6,80];
umax = [4.5,11,105];

plant = @(u)(WOplantFun(u));
objFun = @(u,y)WOobjFun(u,y);
conFun = @(u,y)(y(:,6)-0.08);

fminconopts = optimoptions('fmincon','Display','off','Algorithm','sqp','FunctionTolerance',1e-10,'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10,'StepTolerance',1e-10);

uOptp = fmincon(@(u)objFun(u,plant(u)),uGuess,[],[],[],[],umin,umax,...
    @(u)deal(conFun(u,plant(u)),[]),fminconopts);

yOptp = plant(uOptp);
objOptp = objFun(uOptp,yOptp);
conOptp = conFun(uOptp,yOptp);

%% 2. Solve modified model NLP at plant optimum
% set up model function
model = @(u)WOmodelFun(u);
du = 0.00001;

% modify the objective function
epio = -objFun(uOptp,model(uOptp))+objFun(uOptp,plant(uOptp));
lamo(1) = ((objFun(uOptp+[du,0,0],plant(uOptp+[du,0,0]))-objFun(uOptp-[du,0,0],plant(uOptp-[du,0,0]))-objFun(uOptp+[du,0,0],model(uOptp+[du,0,0]))+objFun(uOptp-[du,0,0],model(uOptp-[du,0,0])))/(2*du));
lamo(2) = ((objFun(uOptp+[0,du,0],plant(uOptp+[0,du,0]))-objFun(uOptp-[0,du,0],plant(uOptp-[0,du,0]))-objFun(uOptp+[0,du,0],model(uOptp+[0,du,0]))+objFun(uOptp-[0,du,0],model(uOptp-[0,du,0])))/(2*du));
lamo(3) = ((objFun(uOptp+[0,0,du],plant(uOptp+[0,0,du]))-objFun(uOptp-[0,0,du],plant(uOptp-[0,0,du]))-objFun(uOptp+[0,0,du],model(uOptp+[0,0,du]))+objFun(uOptp-[0,0,du],model(uOptp-[0,0,du])))/(2*du));

xmodelo = @(u)(objFun(u,model(u))+epio+...
    (u(1)-uOptp(1))'*lamo(1)+...
    (u(2)-uOptp(2))'*lamo(2)+...
    (u(3)-uOptp(3))'*lamo(3));

% modify the constraint function
epic = -conFun(uOptp,model(uOptp))+conFun(uOptp,plant(uOptp));
lamc(1) = ((conFun(uOptp+[du,0,0],plant(uOptp+[du,0,0]))-conFun(uOptp-[du,0,0],plant(uOptp-[du,0,0]))-conFun(uOptp+[du,0,0],model(uOptp+[du,0,0]))+conFun(uOptp-[du,0,0],model(uOptp-[du,0,0])))/(2*du));
lamc(2) = ((conFun(uOptp+[0,du,0],plant(uOptp+[0,du,0]))-conFun(uOptp-[0,du,0],plant(uOptp-[0,du,0]))-conFun(uOptp+[0,du,0],model(uOptp+[0,du,0]))+conFun(uOptp-[0,du,0],model(uOptp-[0,du,0])))/(2*du));
lamc(3) = ((conFun(uOptp+[0,0,du],plant(uOptp+[0,0,du]))-conFun(uOptp-[0,0,du],plant(uOptp-[0,0,du]))-conFun(uOptp+[0,0,du],model(uOptp+[0,0,du]))+conFun(uOptp-[0,0,du],model(uOptp-[0,0,du])))/(2*du));

xmodelc = @(u)(conFun(u,model(u))+epic+...
    (u(1)-uOptp(1))'*lamc(1)+...
    (u(2)-uOptp(2))'*lamc(2)+...
    (u(3)-uOptp(3))'*lamc(3));

% run optimization from plant optimum and far from optimum
uOpt1 = fmincon(@(u)xmodelo(u),uOptp,[],[],[],[],umin,umax,...
    @(u)deal(xmodelc(u),[]),fminconopts);

uOpt2 = fmincon(@(u)xmodelo(u),[3,8.8,94.2],[],[],[],[],umin,umax,...
    @(u)deal(xmodelc(u),[]),fminconopts);

% print KKT points of modified model
fprintf('KKT point 1: u=[%4.1f,%4.1f,%4.1f], con=%5.2f, obj=%5.1f (@plant optimum)\n',uOpt1,xmodelc(uOpt1),xmodelo(uOpt1))
fprintf('KKT point 2: u=[%4.1f,%4.1f,%4.1f], con=%5.2f, obj=%5.1f (@global optimum)\n',uOpt2,xmodelc(uOpt2),xmodelo(uOpt2))

%% 3. Find model hessian

val = zeros(3,3,3,2);

for i = 1:3
    for j = 1:3
        for k = 1:3
            u = uOptp + [i-2,j-2,k-2]*du;
            y = model(u);
            val(i,j,k,1) = objFun(u,y);
            val(i,j,k,2) = conFun(u,y);
        end
    end
end

H = zeros(3,3,2);

for i = 1:2
    H(1,1,i) = (val(3,2,2,i)-2*val(2,2,2,i)+val(1,2,2,i))/(du)/(du);
    H(2,2,i) = (val(2,3,2,i)-2*val(2,2,2,i)+val(2,1,2,i))/(du)/(du);
    H(3,3,i) = (val(2,2,3,i)-2*val(2,2,2,i)+val(2,2,1,i))/(du)/(du);
    H(1,2,i) = (val(3,3,2,i)-val(3,1,2,i)-val(1,3,2,i)+val(1,1,2,i))/(2*du)/(2*du);
    H(1,3,i) = (val(3,2,3,i)-val(3,2,1,i)-val(1,2,3,i)+val(1,2,1,i))/(2*du)/(2*du);
    H(2,3,i) = (val(2,3,3,i)-val(2,3,1,i)-val(2,1,3,i)+val(2,1,1,i))/(2*du)/(2*du);
    H(2,1,i) = H(1,2,i);
    H(3,1,i) = H(1,3,i);
    H(3,2,i) = H(2,3,i);
end

%% 4. Find plant hessian
val = zeros(3,3,3,2);

for i = 1:3
    for j = 1:3
        for k = 1:3
            u = uOptp + [i-2,j-2,k-2]*du;
            y = plant(u);
            val(i,j,k,1) = objFun(u,y);
            val(i,j,k,2) = conFun(u,y);
        end
    end
end

Hp = zeros(3,3,2);

for i = 1:2
    Hp(1,1,i) = (val(3,2,2,i)-2*val(2,2,2,i)+val(1,2,2,i))/(du)/(du);
    Hp(2,2,i) = (val(2,3,2,i)-2*val(2,2,2,i)+val(2,1,2,i))/(du)/(du);
    Hp(3,3,i) = (val(2,2,3,i)-2*val(2,2,2,i)+val(2,2,1,i))/(du)/(du);
    Hp(1,2,i) = (val(3,3,2,i)-val(3,1,2,i)-val(1,3,2,i)+val(1,1,2,i))/(2*du)/(2*du);
    Hp(1,3,i) = (val(3,2,3,i)-val(3,2,1,i)-val(1,2,3,i)+val(1,2,1,i))/(2*du)/(2*du);
    Hp(2,3,i) = (val(2,3,3,i)-val(2,3,1,i)-val(2,1,3,i)+val(2,1,1,i))/(2*du)/(2*du);
    Hp(2,1,i) = Hp(1,2,i);
    Hp(3,1,i) = Hp(1,3,i);
    Hp(3,2,i) = Hp(2,3,i);
end

%% 5. Lagrangian calculation
% plant lagrandian modifier
dplant = zeros(3,2);
u1 = uOptp + [du,0,0];
u2 = uOptp - [du,0,0];
dplant(1,1) = (objFun(u1,plant(u1))-objFun(u2,plant(u2)))/(2*du);
dplant(1,2) = (conFun(u1,plant(u1))-conFun(u2,plant(u2)))/(2*du);

u1 = uOptp + [0,du,0];
u2 = uOptp - [0,du,0];
dplant(2,1) = (objFun(u1,plant(u1))-objFun(u2,plant(u2)))/(2*du);
dplant(2,2) = (conFun(u1,plant(u1))-conFun(u2,plant(u2)))/(2*du);

u1 = uOptp + [0,0,du];
u2 = uOptp - [0,0,du];
dplant(3,1) = (objFun(u1,plant(u1))-objFun(u2,plant(u2)))/(2*du);
dplant(3,2) = (conFun(u1,plant(u1))-conFun(u2,plant(u2)))/(2*du);

lv = -mean(dplant(:,1)./dplant(:,2));

% get legrangian of plant and model
L = H(:,:,1)+lv*H(:,:,2);
Lp = Hp(:,:,1)+lv*Hp(:,:,2);

fprintf('Plant Lagrangian:\n')
fprintf('%6.2f %6.2f %6.2f\n',Lp)
fprintf('Model Lagrangian:\n')
fprintf('%6.2f %6.2f %6.2f\n',L)

% check SOSC
th = linspace(0,2*pi,101);

v = dplant(:,1)/norm(dplant(:,1));
dv1 = (v([3,1,2])-v([2,3,1]));
dv1 = dv1/norm(dv1);
dv2 = cross(dv1,v);
dv2 = dv2/norm(dv2);

dv = dv1*sin(th).^2+dv2*cos(th).^2;

SOSC = zeros(1,numel(th));
SOSCp = zeros(1,numel(th));

for i = 1:numel(th)
    SOSC(i) = dv(:,i)'*L*dv(:,i);
    SOSCp(i) = dv(:,i)'*Lp*dv(:,i);
end

