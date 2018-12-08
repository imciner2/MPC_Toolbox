function [ Delta, y ] = udb_cplex_comp( Hd, Jd, g, D0, Cd0, N, m, pCon )
%UDB_CPLEX_COMP Solve an LPCC to find the upper dual bound
%
% This function will solve an LPCC given in [1] to find an upper bound on
% the size of the 1-norm for the dual variable vector.
%
% This method makes the complimentarity condition (y'*z) into an
% exclusivity constraint (e.g. only one of y or z can be non-zero).
%
% The full optimization problem it solves is:
%   max  sum(y)
%   s.t. Hd*y + Jd*p - z == -g
%         y + Hd*lam - w == 0
%                   D0*p <= Cd0
%                 ind.*y == 0
%           (1 - ind).*w == 0
%           (1 - ind).*z == 0
%
% This function requires that CPLEX is installed, since it is used to
% solve the optimization problem.
%
% [1] P. Patrinos and A. Bemporad, “An Accelerated Dual Gradient-Projection
% Algorithm for Embedded Linear Model Predictive Control,” IEEE Transactions
% on Automatic Control, vol. 59, no. 1, pp. 18–33, 2014.
%
%
% Usage:
%   [ Delta ] = UDB_CPLEX_COMP( Hd, Jd, g, plb, pub );
%
% Inputs:
%   Hd  - The dual Hessian matrix
%   Jd  - The dual linear matrix multiplying the initial state
%   g   - The vector from the RHS of the primal constraints
%   D0  - LHS for the polyhedral constraint on the initial state
%   Cd0 - RHS for the polyhedral constraint on the initial state
%
% Output:
%   Delta - The upper dual bound
%
%
% Created by: Ian McInerney
% Created on: November 17, 2018
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
[nd,  ~] = size(Hd);    % The number of dual variables
[n0,  n] = size(D0);   % The number of states and initial state constraints

jl = nd/N;


%% Create some useful quantities
% Matrices
Zdd = zeros(nd, nd);
Zpd = zeros(n, nd);
Idd = eye(nd);
Zdn0 = zeros(nd, n0);

% Vectors of the dual size
zd = zeros(nd, 1);
od = ones(nd, 1);
id_p = inf(nd, 1);
id_n = -id_p;

% Vectors of the state size
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

% Add the indicator variables to the problem
bin = char( od'*'B' );
for (i=1:1:nd)
    Ind{i} = ['ind', num2str(i)];
end
prob.addCols( zd, [], zd, od, bin, char(Ind) );


% Add the y variables (the dual multipliers)
con = char( od'*'C' );
for (i=1:1:nd)
    yName{i} = ['y', num2str(i)];
end
prob.addCols( od, [], zd, id_p, con, char(yName) );


% Add the lambda variables
for (i=1:1:nd)
    lamName{i} = ['lam', num2str(i)];
end
prob.addCols( zd, [], id_n, id_p, con, char(lamName) );


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

% The constraints on the initial state
%   D0*p <= Cd0
prob.addRows( in0_n, [Zdn0', Zdn0', Zdn0', D0], Cd0 );

% The complimentarity constraints y(i)*w(i) == 0 and y(i)*z(i) == 0 can be
% modeled as an indicator constraint
%         ind == 1 --> y == 0
%         ind == 1 --> Hd*y + Jd*p >= -g
%         ind == 0 --> Hd*y + Jd*p == -g
%         ind == 0 --> y + Hd*lam == 0

% Create the matrices that will be used for the constraints
Acompy = [ Zdd, Idd, Zdd, Zpd' ];
Acompz = [ Zdd,  Hd, Zdd,  Jd ];
Acompw = [ Zdd, Idd,  Hd, Zpd' ];

% Actually create the constraints
for ( i = 1:1:nd )
    % If the indicator is 1, then y == 0 and Hd*y + Jd*p >= -g
    prob.addIndicators(   i, 1, Acompy(i,:)', 'E',     0 );   % Force y to 0 if the indicator is 1
    prob.addIndicators(   i, 1, Acompz(i,:)', 'G', -g(i) );   % Force z to > -g if the indicator is 1
    
    % If the indicator is 0, then Hd*y + Jd*p == -g and y + Hd*lam == 0
    prob.addIndicators(   i, 0, Acompz(i,:)', 'E', -g(i) );   % Force z to -g if the indicator is 0
    prob.addIndicators(   i, 0, Acompw(i,:)', 'E',     0 );   % Force w to 0 if the indicator is 0
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
    sumRow  = zeros(1, 3*nd+n);
    sumRow( startInd:endInd ) = 1;
    prob.addRows( 0, sumRow, min([jl, (i+1)*m]) );
    
    % Create the constraint enforcing LICQ for this and all prior points
    sumRow  = zeros(1, 3*nd+n);
    sumRow( 1:endInd ) = 1;
    prob.addRows( 0, sumRow, (i+1)*m );
    
    % Create the constraints to only allow one of two parallel constraints
    % to be active
    for (j=1:1:numPar)
        sumRow  = zeros(1, 3*nd+n);
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


if (prob.Solution.status ~= 101)
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
%  v = [ind; y; z; w; lam; p];
ind = v(      1:nd   );
y   = v(   nd+1:2*nd );
lam = v( 2*nd+1:3*nd );
p   = v( 3*nd+1:end  );

% Retrieve the upper dual bound from the objective value
Delta = prob.Solution.objval;

end
