function [ Delta, lam ] = udb_cplex_milp( Hp, Jp, G, F, g, D0, Cd0, N, pCon )
%UDB_CPLEX_MILP Solve a MILP to find the upper dual bound
%
% This function will solve a MILP given in [1] to find an upper bound on
% the size of the 1-norm for the dual variable vector. Note that this
% function operates on a infinity-norm bounded set for p.
%
% The MILP solved is:
%   max     sum(y)
%   s.t.  D0*p <= Cd0
%         del(i) = 0 ==> y(i) = 0
%         del(i) = 1 ==> Hd(i,:)*y + Jd(i,:)*p = -g(i)
%
% This function requires that CPLEX is installed, since it is used to
% solve the optimization problem.
%
%   [1] A. Bemporad and P. Patrinos, “Simple and Certifiable Quadratic
%   Programming Algorithms for Embedded Linear Model Predictive Control,”
%   in 4th IFAC Nonlinear Model Predictive Control Conference, 2012,
%   pp. 14–20.
%
%
% Usage:
%   [ Delta ] = UDB_CPLEX_MILP( Hd, Jd, g, L, plb, pub, N, parCon );
%
% Inputs:
%   Hd   - The dual Hessian matrix
%   Jd   - The dual linear matrix multiplying the initial state
%   g    - The vector from the RHS of the primal constraints
%   L    - The largest eigenvalue of the Hessian matrix
%   D0  - LHS for the polyhedral constraint on the initial state
%   Cd0 - RHS for the polyhedral constraint on the initial state
%   N    - The horizon length
%   pcon - 
%
% Output:
%   Delta - The upper dual bound
%
%
% Created by: Ian McInerney
% Created on: November 19, 2018
% Version: 1.0
% Last Modified: November 19, 2018
%
% Revision History
%   1.0 - Initial release


%% Make sure that CPLEX is installed
if ( exist('Cplex', 'class') ~= 8 )
    error('CPLEX not found. Please install CPLEX and make sure it is on the path to run this function.');
end


%% Extract variable size information
[nd,  ~] = size(G);     % The number of dual variables
[ ~, np] = size(Hp);    % The number of primal variables
[n0,  n] = size(D0);    % The number of states and initial state constraints

jl = nd/N;
m  = np/N;


%% Create some useful quantities
% Matrices
Zdd = zeros(nd, nd);
Zpd = zeros(np, nd);
Idd = eye(nd);
Zdn = zeros(nd, n);
Zdn0 = zeros(nd, n0);
Zpn0 = zeros(np, n0);

% Vectors of the dual size
zd = zeros(nd, 1);
od = ones(nd, 1);
id_p = inf(nd, 1);
id_n = -id_p;

% Vectors of the primal size
zp = zeros(np, 1);
op = ones(np, 1);
ip_p = inf(np, 1);
ip_n = -ip_p;

% Vectors of the initial state size
zn = zeros(n, 1);
on = ones(n, 1);
in_p = inf(n, 1);
in_n = -in_p;

% Vectors of the inital constraint set size
in0_p = inf(n0, 1);
in0_n = -in0_p;


%% Create the CPLEX object
prob = Cplex('prob');
prob.Model.sense = 'maximize';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the variables to the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The variable order is:
%   [ ind, lambda, z, p]

% Add the indicator variables to the problem
bin = char( od'*'B' );
for (i=1:1:nd)
    ind{i} = ['ind', num2str(i)];
end
prob.addCols( zd, [], zd, od, bin, char(ind) );

% Add the lambda variables (the dual multipliers)
con = char( od'*'C' );
for (i=1:1:nd)
    lamName{i} = ['lam', num2str(i)];
end
prob.addCols( od, [], zd, 100*od, con, char(lamName) );
%id_p

% Add the z variables (the primal variables)
con = char( op'*'C' );
for (i=1:1:np)
    zName{i} = ['z', num2str(i)];
end
prob.addCols( zp, [], ip_n, ip_p, con, char(zName) );


% Add the p variables (the initial condition)
con = char( on'*'C' );
for (i=1:1:n)
    pName{i} = ['p', num2str(i)];
end
prob.addCols( zn, [], in_n, in_p, con, char(pName) );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the constraints to the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The variable order is:
%   [ ind, lambda, z, p]
%

% The Lagrangian constraint is:
%   Hp*z + Jp'*p + G'*lam == 0
prob.addRows( zp, [Zpd, G', Hp, Jp], zp );

% The inequality constraint is:
%   G*z - F*p <= g
prob.addRows( id_n, [Zdd, Zdd, G, -F], g );

% The constraints on the initial state
%   D0*p <= Cd0
prob.addRows( in0_n, [Zdn0', Zdn0', Zpn0', D0], Cd0 );


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the indicator constraints to the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The variable order is:
%   [ ind, lambda, z, p]
%
% The indicator constraints:
%   Constraint is inactive:
%    ind == 0    ===>   lam = 0
%    ind == 0    ===>   G*z - F*p <= g (implied through other constraint)
%   Constraint is active:
%    ind == 1    ===>   lam >= 0 (implied through the lower bound)
%    ind == 1    ===>   G*z - F*p = g


% Create the matrices that will be used for the constraints
Aind1 = [ Zdd, Zdd,    G,  -F ];
Aind3 = [ Zdd, Idd, Zpd', Zdn ];

% Actually create the indicator
for ( i=1:1:nd )
    % The 2nd term is if the indicator is complimented (e.g. 1 means it
    % takes the value if the indicator is 0).
    prob.addIndicators( i, 1, Aind3(i,:)', 'E', 0 );
    
    prob.addIndicators( i, 0, Aind1(i,:)', 'E', g(i) );
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the LICQ enforcement constraints to the problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

numPar = length(pCon);

% Iterate over each horizon point to create the constraint
for ( i=0:1:N-1 )
    startInd = i*jl+1;
    endInd   = (i+1)*jl;
    
    % Create the constraint enforcing the LICQ for this horizon point
    sumRow  = zeros(1, 2*nd+np+n);
    sumRow( startInd:endInd ) = 1;
    prob.addRows( 0, sumRow, min([jl, (i+1)*m]) );
    
    % Create the constraint enforcing LICQ for this and all prior points
    sumRow  = zeros(1, 2*nd+np+n);
    sumRow( 1:endInd ) = 1;
    prob.addRows( 0, sumRow, (i+1)*m );
    
    % Create the constraints to only allow one of two parallel constraints
    % to be active
    for (j=1:1:numPar)
        sumRow  = zeros(1, 2*nd+np+n);
        sumRow( (i*jl)+pCon{j}(1) ) = 1;
        sumRow( (i*jl)+pCon{j}(2) ) = 1;
        prob.addRows( 0, sumRow, 1 );
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the optimization problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure so that no output happens
prob.Param.mip.display.Cur = 0;
prob.Param.output.clonelog.Cur = 0;
prob.DisplayFunc = [];

% The problem should have a solution that is greater than m
%prob.Param.mip.tolerances.lowercutoff.Cur = nd;

% Solve the problem
prob.solve();

% Make sure the solution is valid
% 101 means the full solution was found
% 102 means that the problem has terminated with a solution very close to
% the optimal
if ( (prob.Solution.status ~= 101) && (prob.Solution.status ~= 102) )
   warning(['CPLEX error: ', prob.Solution.statusstring]);
   Delta = NaN;
   return;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v = prob.Solution.x;

% Variable ordering:
%   [ ind, lambda, z, p]
indw = v(         1:nd      );
lam  = v(      nd+1:2*nd    );
z    = v(    2*nd+1:2*nd+np );
p    = v( 2*nd+np+1:end  );

% Retrieve the upper dual bound from the objective value
Delta = prob.Solution.objval;

end
