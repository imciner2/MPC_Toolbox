% Condensed Formulation for LTI MPC
%
% This folder contains functions related to the creation of the matrices for
% the condensed linear time-invariant MPC formulation.
% 
% This formulation can be solved via the following QP (optimize over the 
% u variables) in the primal form:
%   min  u'Hu + x0'Ju
%   s.t. Gu <= Fx0 + g
%
% Or in the dual form (optimize over y)
%   min  y'Hdy + ( Jdxo + g)'y
%   s.t. 0 <= y
%
% To recover the primal solution from the dual solution:
%   u = -inv(Hp)*(G'y + Jx0)

% Condensed problem formulation:
%   condensed_primal_cost_gen       - Generate the matrices for the condensed primal cost
%   condensed_primal_constraint_gen - Generate the matrices for the condensed constraints

% Dual problem formulation:
%   condensed_dual_cost_gen         - Generate the cost matrices for the dual problem
%   condensed_dual_primal_recovery  - Recover the primal solution from the dual optimal

% Stacked system matrix generation:
%   condensed_initial_gen    - Generate the initial state mapping matrix
%   condensed_prediction_gen - Generate the full condensed prediction matrix