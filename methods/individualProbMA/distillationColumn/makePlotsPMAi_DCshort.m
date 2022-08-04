% script to run worst case MA on WO and plot results
clearvars -except PMAifig

% add path
addpath('../../../plotFunctions/');
addpath('../../../plotFunctions/colours/');
addpath('../../../caseStudies/diehlDistillationColumn/functions/')
addpath('../../../caseStudies/diehlDistillationColumn/analysis/')
addpath('../../modifierAdaptation/distillationColumn/')

%% 0.a Set-up
r0 = [70.001,85];   % starting point
k_max = 11;          % number of iterations
Ns = 14;            % number of stripping trays
Nr = 7:11;          % number of rectifying trays
mu = 3;             % which rectifying tray is nominal

% lines
c_PMAi = brightBlue;
m1_PMAi = {'-s','color',c_PMAi,'linewidth',3,'markersize',7,'MarkerFaceColor',c_PMAi};
m2_PMAi = {'--s','color',c_PMAi,'linewidth',3,'markersize',7,'MarkerFaceColor',c_PMAi};

%% 0.b Get plant optimum
[rpOpt,upOpt,ypOpt,objpOpt,conpOpt] = getPlantOpt;


%% 1 Run WCMA
% run at K=0.9;
[rk_MA1,uk_MA1,yk_MA1,objk_MA1,conk_MA1] = runPMAi_DCshort('r_start',r0,'k_max',k_max,'Ns',Ns,'Nr',Nr,'filter',0.8,'th_nomi',mu);

%% 2 Plot 
if exist('fig','var')
    newFig = 0;
else
    newFig = 1;
end

for i = 1:3
    if newFig || ~isvalid(fig{i})
        fig{i} = figure;
    else
        clf(fig{i})
        ax{i} = [];
    end
    ax{i} = axes(fig{i});
    hold(ax{i},'on')
    fixAxis(fig{i},ax{i})
    p = fig{i}.Position;
    set(fig{i},'Position',[p(1:2),p(3:4)-100])
end

setupPlotsDC(ax,m1_PMAi,m2_PMAi)
plotPlantOpt(ax,k_max,rpOpt,objpOpt,conpOpt)

plot(ax{1},0:k_max-1,(rk_MA1(:,1)-70)/4,m1_PMAi{:})
plot(ax{1},0:k_max-1,(rk_MA1(:,2)-82)/4,m2_PMAi{:})
plot(ax{2},0:k_max-1,-objk_MA1,m1_PMAi{:})
plot(ax{3},0:k_max-1,conk_MA1(:,1),m1_PMAi{:})
plot(ax{3},0:k_max-1,conk_MA1(:,2),m2_PMAi{:})

setupLegendDC(ax)

saveas(fig{1},'plots/PMAir_DCshort.eps','epsc')
saveas(fig{2},'plots/PMAiobj_DCshort.eps','epsc')
saveas(fig{3},'plots/PMAicon_DCshort.eps','epsc')
