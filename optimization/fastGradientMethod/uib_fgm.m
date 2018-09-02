function [ k, delta ] = uib_fgm( cond, L, N, eps, G, g, varargin )
%UIB_FGM Compute the iteration bound for the FGM method applied to primal condensed MPC
%
% Compute the iteration bound for the Fast Gradient Method when it is
% applied to a linear Model Predictive Control problem in condensed-primal
% form with input constraints of the form Gu <= g. This is the iteration
% number by which the Fast Gradient Method will be guaranteed to converge
% to a solution with suboptimality level eps.
%
% There are multiple methods of calculating the Delta parameter:
%   Cold-starting MPC:
%   * 'Cold-Conservative' - Utilize the upper-bound for delta
%   * 'Cold-Optimal'      - Solve an optimization problem to find delta
%   Warm-starting MPC:
%   * 'Warm-Lyapunov'     - Utilize the bound for when P is from the
%                           discrete-Lyapunov function
%
% Using 'Cold-Optimal' requires the YALMIP toolbox to be installed.
%
%
% This lower bound is presented in 
%   S. Richter, C. N. Jones, and M. Morari, “Computational Complexity
%   Certification for Real-Time MPC With Input Constraints Based on the
%   Fast Gradient Method,” IEEE Transactions on Automatic Control,
%   vol. 57, no. 6, pp. 1391–1403, 2012.
%
%   S. Richter, “Computational complexity certification of gradient
%   methods for real-time model predictive control,” PhD Thesis,
%   ETH Zurich, 2012.
%
%
% Usage:
%   [ k, delta ] = UIB_FGM( cond, L, eps, G, g, N);
%   [ k, delta ] = UIB_FGM( cond, L, eps, G, g, N, 'Cold-Conservative' );
%   [ k, delta ] = UIB_FGM( cond, L, eps, G, g, N, 'Cold-Optimal', H, M );
%   [ k, delta ] = UIB_FGM( cond, L, eps, G, g, N, 'Warm-Lyapunov', sys, Q, R );
%
% Inputs:
%   cond - The condition number of the Hessian matrix
%   L    - The largest eigenvalue of the Hessian matrix
%   eps  - The suboptimality level to solve to
%   G    - The input constraint matrix
%   g    - The input constraint vector
%   N    - The horizon length
%
%   For the Cold-Optimal case
%   H    - The Hessian matrix (only needed when using 'Optimal')
%   M    - The linear term from the cost (only needed when using 'Optimal')
%
%   For the Warm-Lyapunov case
%   sys  - The discrete-time system
%   Q    - The state weighting matrix
%   R    - The input weighting matrix
%
% Output:
%   k     - The lower iteration bound
%   delta - The delta parameter used to compute the bound
%
%
% Created by: Ian McInerney
% Created on: May 18, 2018
% Version: 1.2
% Last Modified: August 30, 2018
%
% Revision History
%   1.0 - Initial release
%   1.1 - Added optimal cold-starting
%   1.2 - Added lyapunov warm-starting


%% Parse the input arguments
p = inputParser;
addOptional(p, 'type', 'Cold-Conservative', @(x) isstring(x) || ischar(x));
addOptional(p, 'o1', NaN);
addOptional(p, 'o2', NaN);
addOptional(p, 'o3', NaN);
addOptional(p, 'o4', NaN);
parse(p,varargin{:});

type = p.Results.type;
o1 = p.Results.o1;
o2 = p.Results.o2;
o3 = p.Results.o3;
o4 = p.Results.o4;


%% Compute the Delta parameter
switch (type)
    case 'Cold-Conservative'
        % Use a conservative bound when cold-starting
        delta = fgm_coldStart_consv_delta(L, N, G, g);
        
    case 'Cold-Optimal'
        % Solve the optimization problem for cold-starting
        % Take the inputs and map the proper matrices
        H = o1;
        M = o2;
        if ( max(max(isnan(H))) || max(max(isnan(M))) )
            error('Must supply H and M if cold-start optimal delta computation is selected.');
        end
        delta = fgm_coldStart_optim_delta(L, N, G, g, H, M);
        
    case 'Warm-Lyapunov'
        % Figure out the asymptotic condition number and use it as the
        % upper bound for delta when warm-starting MPC
        sys = o1;
        Q = o2;
        R = o3;
        if ( max(max(isnan(Q))) || max(max(isnan(R))) )
            error('Must supply sys, Q and R if warm-start lyapunov delta computation is selected.');
        end
        delta = condensed_primal_hessian_cond_lyap(sys, Q, R);
        delta = delta*eps;
        
    otherwise
        % Well, this is bad
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

