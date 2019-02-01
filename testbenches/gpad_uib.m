%% Test using the double integrator example

%% Create the physical system
A = [1, 1;
     0, 1];
B = [0;
     1];

Q = [1, 0;
     0, 0];
R = 0.01;

P = dare(A, B, Q, R);

[n, m] = size(B);

%% Create the constraints
ulb = [ 1;
       -1];

E  = kron(eye(m), ulb );
Ce = [ 1
       1];

% Normalize the constraints so that Ce and Cd are +-1
[~, ni] = size(E);
E  = E./(kron( ones(1, ni), abs(Ce) ));
Ce = Ce./abs(Ce);


%% Create the MPC matrices
N = 2;

[Hp, Jp]  = condensed_primal_cost_gen(N, A, B, Q, R, Q);
[G, F, g] = condensed_primal_constraint_gen(N, A, B, E, Ce);
[Hd, Jd]  = condensed_dual_cost_gen( Hp, Jp, G, F );


L_e = eigs( Hd, 1, 'LM')
L_f = norm( Hd, 'fro')
L_1 = norm( Hd, 1)
L_2 = norm( Hd, 2)

L_a = 2*norm(G, 2)^2 / eigs( Hp, 1, 'SM' )

L = L_e

%L = 411.5


%% Configure the error tolerances
eps_g = 10^(-3);
eps_V = 10^(-2);

ub = 15*ones(n, 1);
lb = -15*ones(n, 1);

jl = ni;
pCon = { [1,2] };%, [3,4], [5,6] };

D_cplex_comp = udb_cplex_comp( Hd, Jd, g, lb, ub, N, m, pCon )
D_cplex_milp = udb_cplex_milp_orig( Hd, Jd, g, L, lb, ub, N, m, jl )
D_cplex_full = udb_cplex_milp( Hp, Jp, G, F, g, lb, ub, N, pCon )
D_bm         = udb_bigM( Hd, Jd, g, L, lb, ub )
% D_eps        = udb_eps( Hd, Jd, g, lb, ub )
% D_penalty    = udb_penalty( Hd, Jd, g, lb, ub )
%D_approx     = udb_quadApprox( Hp, Jp, G, F, g, lb, ub, 1e-4 )

D = D_bm;
%D = 34.90   % This is the value for D given in the paper


%% Compute the GPAD UIB value for comparison
N_V = ceil( sqrt( (2*L) / (eps_V) ) * D ) - 2
N_g = ceil( sqrt( (8*L*D) / (eps_g) ) ) - 2

% This is supposed to be 10718
N_nu = max( N_V, N_g )