function [ u, varargout  ] = fgm_quad( H, b, projOp, L, mu, Imax )


%% Save the residuals
if (nargout == 2)
    r = zeros(Imax, 1);
    saveRes = 1;
else
    saveRes = 0;
end


%% Find the size and create an identity matrix
[n, ~] = size(H);
I = eye(n);


%% Compute the step-size
beta = (sqrt(L) - sqrt(mu)) / (sqrt(L) + sqrt(mu));


%% Create the initial condition
z = ones(n,1);
y = z;


%% Iterate the desired number of times
for (i=1:1:(Imax+1))
    % Compute the new step
    t   = y - (1/L)*(H*y + b);
    
    if (saveRes)
        r(i) = norm(t, 2);
    end
    
    % Project the new step onto the feasible space
    z_n = projOp( t );
    
    % Use acceleration to get the next point
    y_n = (1 + beta)*z_n - beta*z;
    
    % Update the variables for the next iteration
    y = y_n;
    z = z_n;
end


%% Save the sequence for output
u = z;

if (saveRes)
    varargout{1} = r;
end

end