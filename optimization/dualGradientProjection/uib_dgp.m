function [ uib ] = uib_dgp( D, Lp, Ld, eps_g, eps_V, eps_z, eps_xi )
%UIB_DGP Compute the upper iteration bound for Dual Gradient Projection
%
% This function will compute the number of iterations required for a
% fixed-point implementation of the dual gradient projection algorithm to
% reach a desired error tolerance.
%
% The bound computation is based on:
%   P. Patrinos, A. Guiggiani, and A. Bemporad, “Fixed-point dual gradient
%   projection for embedded model predictive control,” in 2013 European
%   Control Conference (ECC), 2013.
%
%
% Usage:
%   [ uib ] = UIB_DGP( D, N, Hd, Jd, g, Lp, Ld, eps_g, eps_V, eps_z, eps_xi );
%
% Inputs:
%   D      - The Upper Dual Bound for the QP problem
%   Lp     - Largest eigenvalue of the primal Hessian
%   Ld     - Largest eigenvalue of the dual Hessian
%   eps_g  - Upper bound on the constraint satisfaction error
%   eps_V  - Suboptimality level for the primal problem
%   eps_z  - Suboptimality level for the dual problem
%   eps_xi - Largest error present in the gradient computation for the dual problem
%
% Output:
%   uib - The upper iteration bound
%
%
% Created by: Ian McInerney
% Created on: November 17, 2018
% Version: 1.0
% Last Modified: November 17, 2018
%
% Revision History
%   1.0 - Initial release


%% Check two conditions
c1 = 2*D*eps_xi;
c2 = eps_g*(Lp*eps_z^2 + 2*D*eps_xi) / (eps_g - 2*D*eps_xi);

if ( (eps_g <= c1) || (eps_V <= c2) )
    warning('Conditions on the error have not been met');
    uib = NaN;
    return
end


%% Compute the values of alpha
a = 2*(eps_g + Lp*eps_z^2) / (eps_g - 2*D*eps_xi);
b = (eps_V - Lp*eps_z^2) / ( 2*D*eps_xi );

alpha = min([a, b]);


%% Compute the upper iteration bound
num = Ld*D^2*alpha^2;
den = 2*(eps_g - 2*D*eps_xi)*alpha - 2*(eps_g + Lp*eps_z^2);
uib = (num / den) - 1;


end
