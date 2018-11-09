function [ P, E, k ] = RichterOptimalPrecond( H, m )
%RICHTEROPTIMALPRECOND Compute the optimal preconditioner for an MPC problem
%
% This function computes the optimal preconditioner for the condensed MPC
% problem as described by Ricther in Lemma 9 in
%   S. Richter, C. N. Jones, and M. Morari, “Computational Complexity
%   Certification for Real-Time MPC With Input Constraints Based on the
%   Fast Gradient Method,” IEEE Transactions on Automatic Control,
%   vol. 57, no. 6, pp. 1391–1403, 2012.
%
% The preconditioner matrix P is a left/right reconditioner (e.g. P'*H*P)
% with a minimum eigenvalue of 1 (this is different then specified by
% Richter). 
%
% This function requires the YALMIP toolbox (along with a SDP solver) to
% solve the optimization problem. If the optimization problem errors, then
% NaN is returned.
%
%
% Usage:
%   [ P ] = RICHTEROPTIMALPRECOND( H, m )
%   [ P, E ] = RICHTEROPTIMALPRECOND( H, m )
%   [ P, E, k ] = RICHTEROPTIMALPRECOND( H, m )
%
% Inputs:
%   H - The Hessian matrix
%   m - The number of inputs to the dynamical system that is described by H
%
% Outputs:
%   P - The preconditioner matrix
%   E - The raw result of the optimization problem
%   k - The condition number of the preconditioned matrix
%
%
% Created by: Ian McInerney
% Created on: June 5, 2018
% Version: 1.2
% Last Modified: August 17, 2018
%
% Revision History
%   1.0 - Initial release  
%   1.1 - Made E block diagonal in the solver
%   1.2 - Renamed file


%% Make sure thayt YALMIP is installed
if ( exist('yalmiptest', 'file') ~= 2 )
    error('YALMIP not found. Please install YALMIP and make sure it is on the path to run this function.');
end


%% Set the smallest eigenvalue of the preconditioner
mu = 1;


%% Get the size of H
[n, ~] = size(H);
numBlocks = (n/m);


%% Create the variables
t = sdpvar(1,1);
E = sdpvar(m, m);
E = kron(eye(numBlocks), E);
I = eye(n);


%% Find the lower Cholesky decomposition of H
C = chol(H, 'lower');


%% Create the constraints
F = [ H - mu*E >= 0,
      [ E, C;
       C', t*I] >= 0,
      E >= 0];


%% Call the optimize routine
ops = sdpsettings('verbose',0);
diagStruct = optimize(F, t, ops);

if (diagStruct.problem ~= 0)
    yalErr = yalmiperror(diagStruct.problem);
    warning(['YALMIP error: ', yalErr]);
    
    P = NaN;
    E = NaN;
    k = NaN;
    
    return
end


%% Compute the preconditioning matrix
E = value(E);
P = [];

for ( i=1:1:numBlocks )
    % Extract the next block on the diagonal for analysis
    startInd = ((i-1)*m)+1;
    stopInd = (i*m);
    Esub = E(startInd:stopInd, startInd:stopInd);
    [V, D] = eig(Esub);
    
    % Manipulate the eigenvalues
    lam = sqrt( 1./diag(D) );
    
    % Put the matrix back together
    Part = V*diag(lam)*V';
    P = blkdiag(P, Part);
end

end
