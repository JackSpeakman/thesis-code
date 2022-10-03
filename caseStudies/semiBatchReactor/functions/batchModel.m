function [outputs,t_total] = batchModel(inputs,varargin)
% runs the batch reactor for the model
%
% inputs    1-by-12     decision vars (6xlight intensity,6xN flow)
% varargin  cell        cell of new parameters
%
% outputs   n-by-3      concentrations of (X, N, P)
% t_total   1-by-n      time stamps

% set time steps
if nargout == 2
    tStep = 1; % for plotting
else
    tStep = 10; % for measurements
end

% parameters
u_m = 0.0572;
u_d = 0;
K_N = 393.1;
Y_NX = 504.5;
k_m = 0.00016;
k_d = 0.281;
k_s = 178.9;
k_sq = 23.51;
K_Np = 16.89;

% replace certain parameter values
n_in = floor(numel(varargin));
for i = 1:2:n_in
    if strcmp(varargin{i},'u_m')
        u_m = varargin{i+1};
    elseif  strcmp(varargin{i},'u_d')
        u_d = varargin{i+1};
    elseif  strcmp(varargin{i},'K_N')
        K_N = varargin{i+1};
    elseif  strcmp(varargin{i},'Y_NX')
        Y_NX = varargin{i+1};
    elseif  strcmp(varargin{i},'k_m')
        k_m = varargin{i+1};
    elseif  strcmp(varargin{i},'k_d')
        k_d = varargin{i+1};
    elseif  strcmp(varargin{i},'k_s')
        k_s = varargin{i+1};
    elseif  strcmp(varargin{i},'k_sq')
        k_sq = varargin{i+1};
    elseif  strcmp(varargin{i},'K_Np')
        K_Np = varargin{i+1};
    end
end

% pre-allocate outputs
t_total = 0:tStep:240;
outputs = zeros(numel(t_total),3);
outputs(1,:) = [1,150,0];

for i = 0:5 % run for each 40h step
    % set the inputs (I,F_N)
    in = inputs(i+[1,7]);
    
    % pre-calc value used in dae
    proc = [u_m*(in(1))/(in(1)+k_s),...
        k_m*(in(1))/(in(1)+k_sq)];
    
    %     % run ode solver
    %     [~,y] = ode45(@(t,y)dae(t,y,in),0:40,outputs(1+i*40/tStep,:),options);
    %
    %     % save outputs
    %     outputs((1:40/tStep+1)+40/tStep*i,:) = y((0:tStep:40)+1,:);
    tfinal = (i+1)*40;
    if i == 0
        tinitial = i*40;
        sol = ode45(@(t,y)dae(t,y,in),[tinitial,tfinal],outputs(1+i*40/tStep,:));
    else
        sol = odextend(sol,@(t,y)dae(t,y,in),tfinal);
    end
end
outputs = deval(sol,t_total)';

    function dydt = dae(~,y,in)
        % dae to be solved
        dydt = [proc(1)*y(1)*(y(2))/(y(2)+K_N) - u_d*y(1);...
            -Y_NX*proc(1)*y(1)*(y(2))/(y(2)+K_N) + in(2);...
            proc(2)*y(1) - (k_d*y(3))/(y(2)+K_Np)];
    end

end

