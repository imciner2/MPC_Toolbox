function [ c, f, varargout ] = fgm_fp_fraclength( sys, N, Q, R, roundMode )
%FGM_FP_FRACLENGTH Compute the fraction length needed for FGM in fixed-point
%
% This function will compute the minimum number of fraction bits needed
% when implementing the Fast Gradient Method in embedded MPC.
%
% The input roundMode specifies the type of rounding that occurs after an
% arithmetic operation. Valid rounding modes are:
%   * 'nearest' - Round to nearest
%   * 'zero'    - Round to zero
%   * 'ceil'    - Round to +infinity
%   * 'floor'   - Round to -infinity
%
% In the round to nearest and round to zero modes, the horizon length can
% be set to infinity to compute a horizon-independent fraction length.
%
%
% Usage:
%   [ c, f ] = FGM_FP_FRACLENGTH( sys, N, Q, R, roundMode );
%   [ c, f, eta ] = FGM_FP_FRACLENGTH( sys, N, Q, R, roundMode );
%
% Inputs:
%   sys - The discrete-time system being predicted
%   N   - Horizon length
%   Q   - The state weighting matrix
%   R   - The input weighting matrix
%   roundMode - The rounding mode to use
%
% Outputs:
%   c   - Scaling factor for the Hessian matrix
%   f   - Number of fraction bits needed
%   eta - Rounding stability margin for H
%
%
% See also: FGM_FP_INTLENGTH
%
% Created by: Ian McInerney
% Created on: March 9, 2019
% Version: 1.0
% Last Modified: March 9, 2019
%
% Revision History
%   1.0 - Initial release

%% Extract the matrices
sys = ss(sys);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;

% Compute the terminal weight matrix
P = dlyap(A', Q);

%% Create some helper variables
[k, m] = size(B);
I = eye(k);


%% Make sure it is a valid rounding mode
if ( ~strcmp(roundMode, 'nearest') && ~strcmp(roundMode, 'zero') && ~strcmp(roundMode, 'floor') && ~strcmp(roundMode, 'ceil') )
    error('Unsupported rounding mode');
end


%% Compute the scaling factor needed to make the normal Hessian in (0,1)
% Find the maximum and minimum eigenvalue
[ ~, maxE, minE ] = condensed_primal_hessian_cond_lyap( sys, Q, R );
if ( minE <= 0 )
    error('Problem is not strictly convex');
end

if ( maxE >= 1 )
    % Compute the scaling factor
    c = 1.1*maxE;
else
    % No scaling is needed
    c = 1;
end

% Compute the scaled weight matrices
Qhat = Q./(c);
Rhat = R./(c);
Phat = P./(c);


%% Create the matrix symbols needed for the analysis
z = tf('z', sys.Ts);

Pgam = z*sys;                                       % Symbol for the prediction matrix
PH = Pgam'*Qhat*Pgam + Rhat;                        % Symbol for the Hessian
Pp = z*ss(A, B, B'*Phat, zeros(m,m), sys.Ts);       % Symbol for the truncation


%% Compute the rounding stability factor
%Find the pseudospectra at the points 0 and 1
try
    pHinv = inv( -PH );
    sp0 = norm( pHinv, inf );
catch
    sp0 = inf;
end

try
    pHinv = inv( eye(m) - PH );
    sp1 = norm( pHinv, inf );
catch
    sp1 = inf;
end

eta = min(1/sp0, 1/sp1);


%% Compute the fraction length
if ( N == inf )
    if ( strcmp(roundMode, 'nearest') )
        epsFunc = @(f) 2^(-(f+1));
    elseif ( strcmp(roundMode, 'zero') )
        epsFunc = @(f) 2^(-f);
    else
        error('Infinite horizon is only supported with nearest or zero rounding');
    end
    
    % Iterate over the diagonal numbers to compute the norm
    for (k=0:1:60)
        if ( k == 0 )
            F = B'*Phat*(A^k)*B + Rhat;
        else
            F = B'*Phat*(A^k)*B;
        end
        Nhinf(k+1) = norm(F, inf);
    end
    
    f = NaN;
    for (f=0:1:60)
        eps = epsFunc(f);
        
        % Find the diagonal which is just above this round-off error
        k = find( Nhinf >= eps, 1, 'last');
        if ( isempty(k) )
            continue;
        end
        k = k-1;
        
         % Compute the contribution of EG to the error
        EG = abs( eps*m*(2*k - 1) );
        
        if ( EG >= eta )
            % If EG is greater than the margin, this won't work
            % Move onto the next diagonal
            continue;
        end
        
        % EG is less than the margin, so see what the tail error is
        % Compute for the round to nearest
        sub1 = 0;
        for (i=0:1:k-1)
            sub1 = sub1 + A^(i)*z^(-i);
        end

        ET = 2*norm( Pp - B'*Phat*sub1*B, inf);
        
        if ( ET + EG < eta )
            % The two errors are less than the margin, save result and end
            break;
        end 
    end
else
    % Compute the horizon-dependent fraction length
    f = ceil( -log2( eta / (m*N) ) );
    
    % If nearest is used, remove 1 bit from the length
    if ( strcmp( roundMode, 'nearest' ) )
        f = f-1;
    end
end

%% Output eta
varargout{1} = eta;

end

