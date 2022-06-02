% Plots the characteristic functions for the plant of the 2-var WO CSTR

%% 1. Set-up variables
% grid space
n_u = 21;
umin = [3,6];
umax = [4.5,11];

u1 = linspace(umin(1),umax(1),n_u);
u2 = linspace(umin(2),umax(2),n_u);
[uu1,uu2] = meshgrid(u1,u2);

% plant
yGuess = [0.08746, 0.38962, 0, 0.29061, 0.10945, 0.10754];
T = 90;
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
%% 3. Plot

contour(uu1,uu2,con(:,:,1),[0,0])
hold on
contour(uu1,uu2,con(:,:,2),[0,0])