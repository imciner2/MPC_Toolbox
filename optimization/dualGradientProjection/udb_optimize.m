function [ Delta, varargout ] = udb_optimize( Hd, Jd, g, varargin )
%UDB_OPTIMIZE Solve an LPCC to find the upper dual bound
%
% This function will solve an LPCC given in [1] to find an upper bound on
% the size of the 1-norm for the dual variable vector.
%
% This method approximates the complimentarity condition (y'*z) in the
% optimization problem by adding the quantity rho*y'*z into the cost
% function instead of the constraint. The value of rho is then a penalty
% term to force the complimentarity condition to be 0.
%
% This function requires the YALMIP toolbox to be installed, since it is
% used to interface to the optimization solvers
%
% [1] P. Patrinos and A. Bemporad, “An Accelerated Dual Gradient-Projection
% Algorithm for Embedded Linear Model Predictive Control,” IEEE Transactions
% on Automatic Control, vol. 59, no. 1, pp. 18–33, 2014.
%
%
% Usage:
%   [ Delta ] = UDB_OPTIMIZE( Hd, Jd, g );
%   [ Delta, cce ] = UDB_OPTIMIZE( Hd, Jd, g, rho );
%
% Inputs:
%   Hd  - The dual Hessian matrix
%   Jd  - The dual linear matrix multiplying the initial state
%   g   - The vector from the RHS of the primal constraints
%   rho - Penalty parameter for the complimentarity condition
%
% Output:
%   Delta - The upper dual bound
%   cce   - The complimentarity error (y'*z)
%
%
% Created by: Ian McInerney
% Created on: November 17, 2018
% Version: 1.0
% Last Modified: November 17, 2018
%
% Revision History
%   1.0 - Initial release


%% Make sure thayt YALMIP is installed
if ( exist('yalmiptest', 'file') ~= 2 )
    error('YALMIP not found. Please install YALMIP and make sure it is on the path to run this function.');
end


%% Parse the input arguments
p = inputParser;
addOptional(p,  'rho', 1000);
parse(p,varargin{:});

% Extract the parameter
rho  = p.Results.rho;


%% Extract variable size information
[nd,  ~] = size(Hd);    % The number of dual variables
[ ~, np] = size(Jd);    % The number of primal variables

%% Create the variables
p   = sdpvar(np, 1);    % The initial state variables
y   = sdpvar(nd, 1);    % The dual variables
z   = sdpvar(nd, 1);
w   = sdpvar(nd, 1);    % A variable introduced to break-up the complementarity constraint
lam = sdpvar(nd, 1);


%% The objective function
c = ones(nd, 1);
Obj = c'*y - rho*y'*z;


%% Create the constraint functions
% The linear constraints
Con = [Hd*y + Jd*p + g - z == 0];

% The non-negativity constraints
Con = [Con;
       0 <= z;
       0 <= y];

% The complimentarity constraints
temp = (y + Hd*lam);
for i=1:1:nd
    Con = [Con;
           y(i)*temp(i) == 0];
end


%% Solve the optimization problem
ops = sdpsettings('verbose',0);
opt = optimize(Con, -Obj, ops);

if (opt.problem ~= 0)
    yalErr = yalmiperror(opt.problem);
    error(['YALMIP error: ', yalErr]);
end


%% Return the bound
Delta = sum( value(y) );
cce = value(y)'*value(z);

if nargout == 2
    varargout{1} = cce;
end

end
