function [ i, bits, bounds ] = fgm_fp_intlength( H, J, G, g, Gx0, gx0, beta, Linv, varargin )
%FGM_FP_INTLENGTH Compute the bounds on the integer portion for fixed-point representation
%   
% This function will compute the minimum number of integer bits needed for
% the fast gradient method to solve a linear MPC problem.
%
% The matrices needed to compute the bounds are those of the QP, which is
% as follows:
%   min  0.5*u'*H*u + x0'*J*u
%   s.t. G*u <= g
% Note that there is no dependence of the constraints on the initial state
% in this formulation. This implies that only the inputs can be constrained
% (due to the need for the simple projection operation).
%
% Additionally, a constraint set for the initial state x0 must be supplied
% as:
%   Gx0*x0 <= gx0
%
% There are two methods for computing the variable bounds:
%   * 'interval' - Use internval arithmetic to bound the variables. Some of
%                  the bounds may be tight, others may not be.
%   * 'optim' - Solve optimization probles to find tighter bounds on the
%               variables
%
% All of the numbers passed in should be in exact arithmetic.
%
% Note: This function requires YALMIP to run.
%
%
% Usage:
%   [ i, bits, bounds ] = FGM_FP_INTLENGTH( H, J, G, g, beta, Linv );
%   [ i, bits, bounds ] = FGM_FP_INTLENGTH( H, J, G, g, beta, Linv, method );
%
% Inputs:
%   H      - The Hessian matrix
%   J      - The linear term of cost function
%   G      - The constraint coefficient matrix
%   g      - The constraint bound matrix
%   Gx0    - The initial condition constraint coefficient matrix
%   gxo    - The initial condition bound matrix
%   beta   - The parameter in the algorithm
%   Linv   - The inverse of the largest eigenvalue of H
%   method - Method to use to compute the bounds
%
% Outputs:
%   i      - The overall number of bits needed to represent the numbers
%   bits   - The number of bits needed to represent each quantity
%   bounds - The bound on each quantity
%
%
% See also: FGM_FP_FRACLENGTH
%
% Created by: Ian McInerney
% Created on: March 10, 2019
% Version: 1.0
% Last Modified: March 10, 2019
%
% Revision History
%   1.0 - Initial release


%% Parse the input
p = inputParser;
addOptional(p, 'method', 'interval', @(x) isstring(x) || ischar(x));
parse(p,varargin{:});
method = p.Results.method;


%% Make sure that YALMIP is installed
if ( exist('yalmiptest', 'file') ~= 2 )
    error('YALMIP not found. Please install YALMIP and make sure it is on the path to run this function.');
end

[n, m] = size(J);


%% Compute the simple bounds
bounds.H    = max( max( abs(H) ) );
bounds.beta = max( abs(beta) );
bounds.Linv = max( abs(Linv) );


%% Compute bounds on the product h = J*x0
x0 = sdpvar(m,1);
con = [ Gx0*x0 <= gx0 ];
obj = norm(J*x0, inf);
ops = sdpsettings('verbose',0);
opt = optimize(con, -obj, ops);

if (opt.problem ~= 0)
    yalErr = yalmiperror(opt.problem);
    error(['YALMIP error: ', yalErr]);
end

bounds.h = abs( value(obj) );


%% Compute bounds on the variable z
z = sdpvar(n,1);
con = [ G*z <= g ];
obj = norm(z, inf);
ops = sdpsettings('verbose',0);
opt = optimize(con, -obj, ops);

if (opt.problem ~= 0)
    yalErr = yalmiperror(opt.problem);
    error(['YALMIP error: ', yalErr]);
end

bounds.z = abs( value(obj) );


%% Compute the intermediate bounds
switch method
    case 'interval'
        % The largest possible change in z
        bounds.zdiff = 2*bounds.z;
        
        % The largest possible term beta*(z_{i+1} - z_i)
        bounds.betamult = bounds.beta*bounds.zdiff;
        
        % The largest possible value of y = z_i + beta*(z_{i+1} - z_i)
        bounds.y = bounds.betamult + bounds.z;
        
        % The largest possible value of the term (H*y)
        bounds.Hy = norm(H, inf)*bounds.y;
        
        % The largest possible value of the term (H*y + h)
        bounds.gradY = bounds.Hy + bounds.h;
        
        % The largest possible value of the term (1/L)*(H*y + h)
        bounds.gradYscale = bounds.gradY*bounds.Linv;
        
        % The largest possible value of t = (1/L)*(H*y + h) + y
        bounds.t = bounds.gradYscale + bounds.y;
    
    case 'optim'
        warning('Not implemented yet');
    otherwise
        error('Unknown method.');
end


%% Compute the bits needed for each bound
bits = structfun( @(x) ceil(log2(x))+1, bounds, 'UniformOutput', 0);


%% Compute the most bits needed
i = max( structfun(@(x) max(x(:)), bits) );


end

