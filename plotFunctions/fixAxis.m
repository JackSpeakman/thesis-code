function fixAxis(fig,ax,varargin)
% fix1dAxis modifies the axis
if nargin<2
    ax = gca;
    if nargin<1
        fig = gcf;
    end
end

% set default values
lw = 2;
fs = 20;
bx = 'off';

% change values
for i = 1:numel(varargin)/2
    if strcmp(varargin{i},'linewidth')
        lw = varargin{i+1};
    elseif strcmp(varargin{i},'fontsize')
        fs = varargin{i+1};
    elseif strcmp(varargin{i},'box')
        bx = varargin{i+1};
    end
end

% set stuff
box(ax,bx)
set(ax,'fontname','Times New Roman')

% move
ch = get(fig,'Children');
try
    gs = ch.GridSize(end:-1:1);
catch
    gs = [1,1];
end

figSize = [200,200] + [600,400].*gs;

pos = get(fig,'position');
if all(pos(1:2) == [489,343])
    pos(1:2) = [300,100];
end
set(fig,'position',[pos(1:2),figSize])

% set other stuff
set(ax,'linewidth',lw)
set(ax,'fontsize',fs)



end

