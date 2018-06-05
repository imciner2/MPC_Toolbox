function [ P, E ] = optimalPrecond( H, m )
%OPTIMALPRECOND Compute the optimal preconditioner for an MPC problem
%
% This function computes the optimal preconditioner for the condensed MPC
% problem as described by Ricther in Lemma 9 in
%   S. Richter, C. N. Jones, and M. Morari, “Computational Complexity
%   Certification for Real-Time MPC With Input Constraints Based on the
%   Fast Gradient Method,” IEEE Transactions on Automatic Control,
%   vol. 57, no. 6, pp. 1391–1403, 2012.
%
% This function requires the YALMIP toolbox (along with a SDP solver) to
% solve the optimization problem.
%
%
% Usage:
%   [ P ] = OPTIMALPRECOND( H, m )
%   [ P, E ] = OPTIMALPRECOND( H, m )
%
% Inputs:
%   A - The Hessian matrix
%   m - The number of inputs to the dynamical system that is described by H
%
% Outputs:
%   P - The preconditioner matrix
%   E - The raw result of the optimization problem
%
%
% Created by: Ian McInerney
% Created on: June 5, 2018
% Version: 1.0
% Last Modified: June 5, 2018
%
% Revision History
%   1.0 - Initial release  


%% Make sure thayt YALMIP is installed
if ( exist('yalmiptest', 'file') ~= 2 )
    error('YALMIP not found. Please install YALMIP and make sure it is on the path to run this function.');
end


%% Find the smallest eigenvalue of H
mu = min( eig(H) );


%% Get the size of H
[n, ~] = size(H);


%% Create the variables
t = sdpvar(1,1);
E = sdpvar(n, n);
I = eye(n);


%% Find the lower Cholesky decomposition of H
C = chol(H, 'lower');


%% Create the constraints
F = [ H - mu*E >= 0,
      [ E, C;
       C', t*I] >= 0,
      E >= 0,];


%% Call the optimize routine
ops = sdpsettings('verbose',0);
optimize(F, t, ops);


%% Compute the preconditioning matrix
E = double(E);
P = [];

for ( i=1:1:(n/m) )
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