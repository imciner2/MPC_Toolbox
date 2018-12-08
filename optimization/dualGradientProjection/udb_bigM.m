function [ Delta, lam ] = udb_bigM( Hp, Jp, G, F, g, D0, Cd0, N, pCon, varargin )
%UDB_BIGM Approximate the Upper Dual Bound using the Big-M MILP formulation
%
% This function will use the Big-M reformulation given in [Section 4.1, 1]
% to find an upper bound on the size of the 1-norm for the dual variable
% vector. The upper/lower bounds on y and w are given by bnd, which
% defaults to 1000. If the problem fails to solve properly, try increasing
% the bound.
%
% This function requires the YALMIP toolbox to be installed, since it is
% used to interface to the optimization solvers.
%
%   [1] A. Bemporad and P. Patrinos, “Simple and Certifiable Quadratic
%   Programming Algorithms for Embedded Linear Model Predictive Control,”
%   in 4th IFAC Nonlinear Model Predictive Control Conference, 2012,
%   pp. 14–20.
%
%
% Usage:
%   [ Delta ] = UDB_BIGM( Hd, Jd, g, L, plb, pub );
%   [ Delta ] = UDB_BIGM( Hd, Jd, g, L, plb, pub, bnd );
%
% Inputs:
%   Hd  - The dual Hessian matrix
%   Jd  - The dual linear matrix multiplying the initial state
%   g   - The vector from the RHS of the primal constraints
%   L   - The largest eigenvalue of the Hessian matrix
%   D0  - LHS for the polyhedral constraint on the initial state
%   Cd0 - RHS for the polyhedral constraint on the initial state
%   bnd - The bound to use for the variables (defaults to 1000)
%
% Output:
%   Delta - The upper dual bound
%
%
% see also: UDB_EPS
%
% Created by: Ian McInerney
% Created on: November 19, 2018
% Version: 1.0
% Last Modified: November 19, 2018
%
% Revision History
%   1.0 - Initial release


%% Make sure thayt YALMIP is installed
if ( exist('yalmiptest', 'file') ~= 2 )
    error('YALMIP not found. Please install YALMIP and make sure it is on the path to run this function.');
end


%% Extract the bound
p = inputParser;
addOptional(p, 'bnd', 1000);
parse(p,varargin{:});
bnd = p.Results.bnd;


%% Extract variable size information
[nd,  ~] = size(G);     % The number of dual variables
[ ~, np] = size(Hp);    % The number of primal variables
[ ~,  n] = size(D0);    % The number of states and initial state constraints

jl = nd/N;
m  = np/N;


%% Create the variables
lam = sdpvar(nd, 1);    % The dual variables
p   = sdpvar( n, 1);    % The initial state variables
z   = sdpvar(np, 1);    % The primal variables
s   = binvar(nd, 1);    % The indicator variables


%% Create some suitable bounds
mlam = 100;
mz = 100;

Mlam = diag( ones(1, nd).*mlam );
Mz   = diag( ones(1, nd).*mz   );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the LICQ enforcement constraint matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

numPar = length(pCon);

% Iterate over each horizon point to create the constraint
licq = [];
ub   = [];
for ( i=0:1:N-1 )
    startInd = i*jl+1;
    endInd   = (i+1)*jl;
    
    % licqRow1 = Enforce LICQ for this horizon point
    % licqRow2 = Enforce LICQ for this and all prior horizon points
    licqRow1 = zeros(1, nd);
    licqRow2 = zeros(1, nd);
    
    licqRow1( startInd:endInd ) = 1;
    licqRow2( 1:endInd ) = 1;
    
    licq = [ licq;
             licqRow1;
             licqRow2];
    
    ub = [ub;
          min([jl, (i+1)*m]);
          (i+1)*m];
    
    % Create the constraints to only allow one of two parallel constraints
    % to be active
    for (j=1:1:numPar)
        symRow  = zeros(1, nd);
        symRow( (i*jl)+pCon{j}(1) ) = 1;
        symRow( (i*jl)+pCon{j}(2) ) = 1;
        
        licq = [licq;
                symRow];
        ub = [ub;
               1];
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Con = [
       Hp*z + Jp*p + G'*lam == 0;         % Lagrangian derivative
                  G*z - F*p <= g;         % Primal inequality constraint
                        lam >= 0;         % Dual variable
                        lam <= Mlam*s;    % Big-M for setting lambda == 0
               G*z - F*p- g >= -Mz*(1-s); % Big-M for setting inequality to equality
                       D0*p <= Cd0;       % Constrain the initial state
                     licq*s <= ub;        % Enforce LICQ
      ];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = ones(nd, 1);
Obj = c'*lam;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the optimization problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ops = sdpsettings('verbose',0);
opt = optimize(Con, -Obj, ops);

if (opt.problem ~= 0)
    yalErr = yalmiperror(opt.problem);
    error(['YALMIP error: ', yalErr]);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve the solution and dual vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delta = value(Obj);
lam = value( lam );

end
