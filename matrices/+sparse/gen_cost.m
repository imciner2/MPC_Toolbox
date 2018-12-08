function [ H ] = gen_cost( N, Q, R, varargin )
%GEN_COST Summary of this function goes here
%   Detailed explanation goes here
%
% The variable vector has the form of:
%   [ x_0, u_0, x_1, u_1, ..., u_(N-1), x_N ]
%
%


[n, ~] = size(Q);
[m, ~] = size(R);

%% Parse the input arguments
p = inputParser;
addOptional(p, 'P', Q);
addOptional(p, 'S', zeros(n,m));
parse(p,varargin{:});

% Extract the matrices
P = p.Results.P;
S = p.Results.S;

% See if P was provided by the user
if (isempty(P))
    P = Q;
end


%% Form the Hessian
T = [ Q, S;
     S', R];

H = kron( speye(N), T);
H = blkdiag( H, P );


end
