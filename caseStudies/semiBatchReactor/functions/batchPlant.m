function [outputs,t_total] = batchPlant(inputs)
% runs the batch reactor for the plant
%
% inputs    1-by-12     decision vars (6xlight intensity,6xN flow)
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
k_i = 447.1;
k_sq = 23.51;
k_iq = 800;
K_Np = 16.89;

% pre-allocate outputs
t_total = 0:tStep:240;
outputs = zeros(numel(t_total),3);
outputs(1,:) = [1,150,0];

for i = 0:5 % run for each 40h step
    % set the inputs (I,F_N)
    in = inputs(i+[1,7]);
    
    % pre-calc value used in dae
    proc = [(in(1))/(in(1)+k_s+(in(1)^2/k_i)),...
        k_m*(in(1))/(in(1)+k_sq+(in(1)^2/k_iq))];
    
    % run ode solver
    [~,y] = ode45(@(t,y)dae(t,y,in),(0:tStep:40)+40*i,outputs(1+i*40/tStep,:));
    
    % save outputs
    outputs((1:40/tStep+1)+40/tStep*i,:) = y;
    
end

    function dydt = dae(~,y,in)
        % dae to be solved
        dydt = [u_m*proc(1)*y(1)*(y(2))/(y(2)+K_N) - u_d*y(1);...
            -Y_NX*u_m*proc(1)*y(1)*(y(2))/(y(2)+K_N) + in(2);...
            proc(2)*y(1) - (k_d*y(3))/(y(2)+K_Np)];
    end

end

