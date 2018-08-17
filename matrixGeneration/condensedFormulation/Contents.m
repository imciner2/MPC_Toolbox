% Condensed Formulation for LTI MPC
%
% This folder contains functions related to the creation of the matrices for
% the condensed linear time-invariant MPC formulation.
% 
% This formulation can be solved via the following QP (optimize over the 
% u variables):
%   min  u'Hu + x0'Ju
%   s.t. Gu <= Fx0 + g
%
% Optimization problem generation:
%   condensed_cost_hessian_gen - Generate the Hessian matrix H
%   condensed_cost_linear_gen  - Generate the linear term matrix J
%
% Intermediate matrix generation:
%   condensed_initial_gen    - Generate the initial state mapping matrix
%   condensed_prediction_gen - Generate the full condensed prediction matrix

