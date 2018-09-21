function [ Hd ] = condensed_dual_cost_hessian_gen( Hp, G )
%CONDENSED_DUAL_COST_HESSIAN_GEN Generate the Hessian matrix H of the dual problem
%
% Create the condensed Hessian matrix for the cost function of the dual form
% of the condensed linear time-invariant MPC problem.
%
% Usage:
%   [ Hd ] = CONDENSED_DUAL_COST_HESSIAN_GEN( Hp, G )
%
% Inputs:
%   Hp - The condensed primal Hessian matrix
%   G  - The condensed primal constraint matrix
%
% Outputs:
%   Hd - The dual Hessian matrix
%
% See also CONDENSED_PRIMAL_COST_HESSIAN_GEN,  CONDENSED_PRIMAL_CONSTRAINT_COEFFICIENT_GEN
%
% Created by: Ian McInerney
% Created on: September 19, 2018
% Version: 1.0
% Last Modified: September 19, 2018
%
% Revision History
%   1.0 - Initial release  


%% Make sure the primal Hessian is positive definite
[~,p] = chol(Hp);

if (p ~= 0)
    error('Hp must be positive definite');
end


%% Compute the dual Hessian
Hd = G*(Hp\(G'));

end
