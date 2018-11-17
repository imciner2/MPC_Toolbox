function [ u ] = condensed_dual_primal_recovery( y, x0, Hp, J, G )
%CONDENSED_DUAL_PRIMAL_RECOVERY Recover the primal solution from the dual optimal
%
% Recover the primal optimal solution from the dual optimal solution.
%
% Primal problem:
%   min  u'Hpu + x0'Ju
%   s.t. Gu <= Fx0 + g
%
% Dual problem:
%   min  y'Hdy + ( Jdxo + g)'y
%   s.t. 0 <= y
%
%
% Usage:
%   [ u ] = CONDENSED_DUAL_COST_GEN( y, x0, Hp, J, G )
%
% Inputs:
%   y  - The dual optimal solution
%   x0 - The initial state
%   Hp - The primal Hessian matrix
%   J  - The primal linear matrix
%   G  - The primal constraint coefficient matrix
%
% Outputs:
%   u - The primal optimal solution
%
%
% See also CONDENSED_PRIMAL_COST_GEN,  CONDENSED_PRIMAL_CONSTRAINT_GEN,
% CONDENSED_DUAL_COST_GEN
%
% Created by: Ian McInerney
% Created on: November 17, 2018
% Version: 1.0
% Last Modified: November 17, 2018
%
% Revision History
%   1.0 - Initial release


%% Make sure the primal Hessian is positive definite
[~,p] = chol(Hp);

if (p ~= 0)
    error('Hp must be positive definite');
end


%% Compute the solution
T = G'*y + J'*x0;
u = Hp\(T);

end
