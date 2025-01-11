//------------------------------------------------------------------------------
//  Simple RBC Model with Endogenous Labor Supply
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//  Parameters
//------------------------------------------------------------------------------
parameters beta, sigma, psi, alpha, delta, rho, sigma_e;
parameters rho_ss, y_k, k_l;  % Useful constants for steady-state computations

% Parameter values
beta    = 0.97;    % Discount factor
sigma   = 1;       % Risk aversion (log utility in consumption)
psi     = 1;       % Inverse Frisch elasticity of labor supply
alpha   = 0.33;    % Capital share in output
delta   = 0.10;    % Depreciation rate of capital
rho     = 0.95;    % Persistence of technology shock
sigma_e = 0.007;   % Std. dev. of technology shock

% Useful constants (steady-state ratios)
rho_ss  = 1 / beta - 1;               % Steady-state interest rate
y_k     = (rho_ss + delta) / alpha;   % Output-capital ratio
k_l     = y_k^(1 / (alpha - 1));      % Capital-labor ratio


//------------------------------------------------------------------------------
//  Variable declaration block
//------------------------------------------------------------------------------
var c, l, k, y, z;  % Consumption, labor, capital, output, productivity
varexo e;          % Technology shock


//------------------------------------------------------------------------------
//  Model equations
//------------------------------------------------------------------------------
model;

% 1. Euler equation (intertemporal optimization)
c^(-sigma) = beta * (alpha * y(+1) / k + (1 - delta)) * c(+1)^(-sigma);

% 2. Labor supply (intratemporal FOC)
l^psi = (1 - alpha) * y / l * c^(-sigma);

% 3. Budget constraint (resource constraint)
c + k = y + (1 - delta) * k(-1);

% 4. Output equation (Cobb-Douglas production function)
y = exp(z) * k(-1)^alpha * l^(1 - alpha);

% 5. Technology process (AR(1) process)
z = rho * z(-1) + e;

end;


//------------------------------------------------------------------------------
//  Shocks block
//------------------------------------------------------------------------------
shocks;
var e; stderr sigma_e;  % Std. dev. of technology shock
end;


//------------------------------------------------------------------------------
//  Steady-state values (derived from the model)
//------------------------------------------------------------------------------
initval;

% Steady-state values (derived from useful constants)
l = ((1 - alpha) * k_l^alpha)^(1 / psi);  % Labor
k = k_l * l;  % Capital
y = y_k * k;  % Output
c = y - delta * k;  % Consumption
z = 0;  % Technology at its mean value

end;


//------------------------------------------------------------------------------
//  Steady-state computation and stability check
//------------------------------------------------------------------------------
steady;  % Compute the steady state

check;  % Check for stability (eigenvalue check)


//------------------------------------------------------------------------------
//  Stochastic simulation (3rd-order perturbation)
//------------------------------------------------------------------------------
stoch_simul(order=3, irf=0);  % 3rd-order approximation, suppress IRFs