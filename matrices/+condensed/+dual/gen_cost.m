function [ Hd, Jd ] = gen_cost( Hp, J, G, F )
%GEN_COST Generate the cost matrices for the dual problem
%
% Create the condensed matrices for the cost function of the dual form
% of the condensed linear time-invariant MPC problem.
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
%   [ Hd, Jd ] = GEN_COST( Hp, J, G, F )
%
% Inputs:
%   Hp - The primal Hessian matrix
%   J  - The primal linear matrix
%   G  - The primal constraint coefficient matrix
%   F  - The primal constraint initial state matrix
%
% Outputs:
%   Hd - The dual Hessian matrix
%   Jd - The dual initial state matrix
%
%
% See also CONDENSED.PRIMAL.GEN_COST,  CONDENSED.PRIMAL.GEN_CON_INEQ
%
% Created by: Ian McInerney
% Created on: September 19, 2018
% Version: 1.1
% Last Modified: November 17, 2018
%
% Revision History
%   1.0 - Initial release
%   1.1 - Added linear term


%% Make sure the primal Hessian is positive definite
[~,p] = chol(Hp);

if (p ~= 0)
    error('Hp must be positive definite');
end


%% Compute the dual Hessian
Hd = G*(Hp\(G'));


%% Compute the dual initial state matrix
Jd = G*(Hp\(J)) + F;

end
