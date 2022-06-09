% script for approximating the parameters.
clear
close all

%% 1. Set-up u values
n_u = 11;
u_range = [0,1.2,1,2];

[uu1,uu2] = meshgrid(linspace(u_range(1),u_range(2),n_u),linspace(u_range(3),u_range(4),n_u));

%% 2. Get plant values
thp = [1,1.5];

conp = iC3con(uu1,uu2,thp(1),thp(2),0);
objp = iC3obj(uu1,uu2,thp(1),thp(2),0);

%% 3. Get model 1 values
f = fittype(@(a,uu1,uu2)(iC3obj(uu1,uu2,a,[],1)),'independent',{'uu1','uu2'});
fitm1th1 = fit([uu1(:),uu2(:)],objp(:),f,'StartPoint',1)

f = fittype(@(b,uu1,uu2)(iC3con(uu1,uu2,[],b,1)),'independent',{'uu1','uu2'});
fitm1th2 = fit([uu1(:),uu2(:)],conp(:),f,'StartPoint',1.5)

%% 3. Get model 2 values
f = fittype(@(a,uu1,uu2)(iC3obj(uu1,uu2,a,[],2)),'independent',{'uu1','uu2'});
fitm2th1 = fit([uu1(:),uu2(:)],objp(:),f,'StartPoint',1)

f = fittype(@(b,uu1,uu2)(iC3con(uu1,uu2,[],b,2)),'independent',{'uu1','uu2'});
fitm2th2 = fit([uu1(:),uu2(:)],conp(:),f,'StartPoint',1.5)