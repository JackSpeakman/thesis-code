function [uk,yk,conk,objk] = runDHFilterMA(varargin)
% The master code for running the MSMA approaches.
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
%   'plantFun'          @(u,th)             Plant function
%   'method'            string              ['minmod','mean','closest']
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
% >> runMSMA
%       This runs the default argument which is MSMA using modifiers applied to WO

% >> addpath('../../caseStudies/williamsOttoCSTR/functions/')
% >> runMSMA('method','closest','filter',0.7,'startingPoint',[3.2,6.5,90])
%       This runs the closest solution method to WO with an input filter gain of 0.7, 
%       and the initial point of [3.2,6.5,90]
% 
% 
% See 'distillationColumn/makePlots___.m' for distillation colunm examples
% See 'semiBatchReactor/makePlots___.m' for semi-batch reactor examples

% ~~~~~~~ TO DO ~~~~~~~~
end