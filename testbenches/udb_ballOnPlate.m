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

[Hp, Jp]  = condensed.primal.gen_cost(N, A, B, Q, R, Q);
[G, F, g] = condensed.primal.gen_con_ineq(N, A, B, E, Ce, D, Cd);
[Hd, Jd]  = condensed.dual.gen_cost( Hp, Jp, G, F );


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

jl = ni + ns;
pCon = { [1,2], [3,4], [5,6] };

D_cplex_comp = udb_cplex_comp( Hd, Jd, g, D, Cd, N, m, pCon )
D_cplex_full = udb_cplex_milp( Hp, Jp, G, F, g, D, Cd, N, pCon )
D_bm         = udb_bigM( Hp, Jp, G, F, g, D, Cd, N, pCon )

%D = 34.90   % This is the value for D given in the paper
