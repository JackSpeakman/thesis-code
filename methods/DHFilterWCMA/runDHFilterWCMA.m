function [uk,yk,conk,objk] = runDHFilterWCMA(varargin)
% The master code for running the DH filter to WCMA.
% The DH filter approaches are described in the I&EC article titled "Robust 
% Modifier Adaptation via Worst-Case and Probabilistic Approaches" by Jack 
% Speakman and Gregory Francois, working off the DHMA constraints as defined
% in the Journal of Process Control article titled "A Robust Modifier 
% Adaptation Method via Hessian Augmentation using Model Uncertainties" by
% Jack Speakman, Aris Papasavvas and Gregory Francois.
% 
% ------ INPUT VARIABLES ------
% varargin                                  cell of inputs
%   'filter'            1-by-1              Input filter
%   'kmax'              1-by-1              Number of RTO iterations
%   'stratingPoint'     1-by-n_u            RTO starting point
%   'uGuess'            1-by-n_u            Input guess (for initial runs)
%   'th'                i-by-n_th           Set of model parameters
%   'th_nom'            i-by-n_th           Nominal model parameters
%   'conFun'            @(u,y)              Constraint function
%   'objFun'            @(u,y)              Objective function
%   'modelFun'          @(u,th)             Model function
%   'umin'              1-by-n_u            Optimization minimum limit
%   'umax'              1-by-n_u            Optimization maximum limit
% 
% ------ OUTPUT VARIABLES ------
% uk        kmax-by-n_u         Inputs for iterations 1 to kmax
% yk        kmax-by-n_y         Outputs for iterations 1 to kmax
% conk      kmax-by-n_c         Constraints for iterations 1 to kmax
% objk      kmax-by-1           Objective for iterations 1 to kmax
% 
% ------ EXAMPLES ------
% >> addpath('../../caseStudies/williamsOttoCSTR/functions/')
% >> runDHFilterWCMA
%       This runs the default argument which is DH filter applied to WO
% 
% 
% See 'distillationColumn/makePlots___.m' for distillation colunm examples
% See 'semiBatchReactor/makePlots___.m' for semi-batch reactor examples

% ~~~~~~~ TO DO ~~~~~~~~
end