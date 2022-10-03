function [uk,yk,conk,objk] = runDHMA(varargin)
% The master code for running the DH filter to standard MA.
% The DHMA methods are run as defined in the Journal of Process Control article 
% titled "A Robust Modifier Adaptation Method via Hessian Augmentation using 
% Model Uncertainties" by Jack Speakman, Aris Papasavvas and Gregory Francois.
% 
% ------ INPUT VARIABLES ------
% varargin                                  cell of inputs
%   'filter'            1-by-1              Input filter
%   'kmax'              1-by-1              Number of RTO iterations
%   'stratingPoint'     1-by-n_u            RTO starting point
%   'uGuess'            1-by-n_u            Input guess (for initial runs)
%   'th_nom'            i-by-n_th           Nominal model parameters
%   'conFun'            @(u,y)              Constraint function
%   'objFun'            @(u,y)              Objective function
%   'plantFun'          @(u,th)             Plant function
%   'umin'              1-by-n_u            Optimization minimum limit
%   'umax'              1-by-n_u            Optimization maximum limit
%   'datafile'          string              Data filename
% 
% ------ OUTPUT VARIABLES ------
% uk        kmax-by-n_u         Inputs for iterations 1 to kmax
% yk        kmax-by-n_y         Outputs for iterations 1 to kmax
% conk      kmax-by-n_c         Constraints for iterations 1 to kmax
% objk      kmax-by-1           Objective for iterations 1 to kmax
% 
% ------ EXAMPLES ------
% >> addpath('../../caseStudies/williamsOttoCSTR/functions/')
% >> runDHMA
%       This runs the default argument which is DHMA applied to WO
% 
% 
% See 'distillationColumn/makePlots___.m' for distillation colunm examples
% See 'semiBatchReactor/makePlots___.m' for semi-batch reactor examples

% ~~~~~~~ TO DO ~~~~~~~~
end