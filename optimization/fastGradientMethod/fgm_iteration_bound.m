function [ k, delta ] = fgm_iteration_bound( cond, L, N, eps, G, g, varargin )
%FGM_ITERATION_BOUND Compute the iteration bound for the FGM method applied to MPC
%
% Compute the lower iteration bound for the Fast Gradient Method when it is
% applied to a linear Model Predictive Control problem with input constraints
% of the form Au <= b.
%
% There are two methods of calculating the Delta parameter
%   * 'Conservative' - Utilize the upper-bound for delta
%   * 'Optimal-Cold' - Solve an optimization problem to find delta
%
% Using 'Optimal' requires the YALMIP toolbox to be installed.
%
%
% This lower bound is derived in 
%   S. Richter, C. N. Jones, and M. Morari, “Computational Complexity
%   Certification for Real-Time MPC With Input Constraints Based on the
%   Fast Gradient Method,” IEEE Transactions on Automatic Control,
%   vol. 57, no. 6, pp. 1391–1403, 2012.
%
%
% Usage:
%   [ k, delta ] = FGM_ITERATION_BOUND( cond, L, eps, A, b N);
%   [ k, delta ] = FGM_ITERATION_BOUND( cond, L, eps, A, b, N, 'Conservative' );
%   [ k, delta ] = FGM_ITERATION_BOUND( cond, L, eps, A, b, N, 'Optimal-Cold', H, M );
%
% Inputs:
%   cond - The condition number of the Hessian matrix
%   L    - The largest eigenvalue of the Hessian matrix
%   eps  - The suboptimality level to solve to
%   G    - The input constraint matrix
%   g    - The input constraint vector
%   N    - The horizon length
%   H    - The Hessian matrix (only needed when using 'Optimal')
%   M    - The linear term from the cost (only needed when using 'Optimal')
%
% Output:
%   k     - The lower iteration bound
%   delta - The delta parameter used to compute the bound
%
%
% Created by: Ian McInerney
% Created on: May 18, 2018
% Version: 1.1
% Last Modified: August 13, 2018
%
% Revision History
%   1.0 - Initial release
%   1.1 - Added optimal cold-starting

%% Parse the input arguments
p = inputParser;
addOptional(p, 'type', 'Conservative', @(x) isstring(x) || ischar(x));
addOptional(p, 'H', NaN);
addOptional(p, 'M', NaN);
parse(p,varargin{:});

type = p.Results.type;
H = p.Results.H;
M = p.Results.M;


%% Compute the Delta parameter
switch (type)
    case 'Conservative'
        % Call the conservative bounder
        delta = fgm_coldStart_consv_delta(L, N, G, g);
    case 'Optimal-Cold'
        % Call the optimal solver
        if ( max(max(isnan(H))) || max(max(isnan(M))) )
            error('Must supply H and M if optimal delta computation is selected.');
        end
        
        delta = fgm_coldStart_optim_delta(L, N, G, g, H, M);
    otherwise
        error('Unknown choice for delta computation.');
end


%% Compute two possible bounds
a1 = ceil( ( log(eps) - log(delta) ) / ( log(1 - sqrt(1/cond) ) ) );
a2 = ceil( 2*sqrt(delta/eps) - 2);


%% Choose the right bound
k = min( [a1, a2] );


%% Floor the iteration count to 0
k = max( [0, k] );

end

