%% Test using the ball on plate example

%% Create the physical system
A = [1, 0.01;
     0,    1];
B = [-0.0004;
     -0.0701];

Q = [100,  0;
       0, 10];
R = 1;

P = dare(A, B, Q, R);

[n, m] = size(B);

%% Create the constraints
ulb = [1;
      -1];

E  = kron(eye(m), ulb );
Ce = [0.0524
      0.0524];
  
D  = kron(eye(n), ulb );
Cd = [0.01;
      0.2
      0.1
      0.1];

% Normalize the constraints so that Ce and Cd are 1
[ni, ~] = size(E);
E  = E./(kron( ones(1, m), abs(Ce) ));
Ce = Ce./abs(Ce);

[ns, ~] = size(D);
D  = D./(kron( ones(1, n), abs(Cd) ));
Cd = Cd./abs(Cd);


%% Bound the initial state
ub = [ Cd(1);
       Cd(3)];
lb = [ -Cd(2);
       -Cd(4)];
    

%% Create the MPC matrices
N = 5;

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

jl = ni;
pCon = { [1,2] };%, [3,4], [5,6] };

D_cplex_comp = udb_cplex_comp( Hd, Jd, g, lb, ub, N, m, pCon )
D_cplex_milp = udb_cplex_milp_orig( Hd, Jd, g, L, lb, ub, N, m, jl )
D_cplex_full = udb_cplex_milp( Hp, Jp, G, F, g, lb, ub, N, pCon )
D_bm         = udb_bigM( Hd, Jd, g, L, lb, ub )
%D_eps        = udb_eps( Hd, Jd, g, lb, ub )
%D_penalty    = udb_penalty( Hd, Jd, g, lb, ub )
%D_approx     = udb_quadApprox( Hp, Jp, G, F, g, lb, ub, 1e-4 )

D = D_bm;
%D = 34.90   % This is the value for D given in the paper

UDB = D_cplex_full;


%% Configure the error tolerances
eps_g  = 1e-2;    % Constraint satisfaction error
eps_V  = 1e-2;    % Primal suboptimality level
eps_z  = 0;       % Dual suboptimality level
eps_xi = 0;       % Dual gradient error computation


%% Compute the DGP UIB value for comparison
UIB = uib_dgp( Hd, Jd, g, Lp, Ld, plb, pub, eps_g, eps_V, eps_z, eps_xi )
