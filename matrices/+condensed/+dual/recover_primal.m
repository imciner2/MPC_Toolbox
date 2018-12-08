function [ u ] = recover_primal( y, x0, Hp, J, G )
%RECOVER_PRIMAL Recover the primal solution from the dual optimal
%
% Recover the primal optimal solution from the dual optimal solution.
%
% Primal problem:
%   min  0.5*u'*Hp*u + x0'*J*u
%   s.t. Gu <= Fx0 + g
%
% Dual problem:
%   min  0.5*y'*Hd*y + ( Jd*xo + g)'*y
%   s.t. 0 <= y
%
%
% Usage:
%   [ u ] = RECOVER_PRIMAL( y, x0, Hp, J, G )
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
% See also CONDENSED.PRIMAL.GEN_COST,  CONDENSED.PRIMAL.GEN_CON_INEQ,
% CONDENSED.DUAL.GEN_COST
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
T = G'*y + J*x0;
u = Hp\(T);


%% Negate the solution
u = -u;

end
