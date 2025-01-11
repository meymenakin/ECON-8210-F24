%-------------------------------------------------------------------------------
% Plot Policy Functions from Dynare Output
% Purpose: Creates plots of policy functions for c, k, and l using 3rd order
%          perturbation solution from Dynare
%-------------------------------------------------------------------------------
clear variables; close all; clc;

addpath /Applications/Dynare/6.2-arm64/matlab

% Run the dynare model
dynare econ8210hw2q3_meakin.mod;

%-------------------------------------------------------------------------------
% Set up grids for state variables
%-------------------------------------------------------------------------------
% Parameters for grid
kssfac = 0.25;     % capital grid range (±25% around steady state)
zsigma = 3;        % productivity grid range (±3 std dev)

% Steady states from Dynare
kss = oo_.steady_state(3);  % capital is third variable
zss = 0;                    % productivity steady state

% Create grids as deviations from steady state
nk = 100;  % number of points for capital grid
nz = 3;   % number of points for productivity grid

% Capital grid
kmin = (1-kssfac)*kss;
kmax = (1+kssfac)*kss;
kg = linspace(kmin, kmax, nk)';
dK = kg - kss;

% Productivity grid
sig_z = sigma_e/sqrt(1-rho^2);  % unconditional std dev of productivity
zmin = -zsigma*sig_z;
zmax = zsigma*sig_z;
zg = linspace(zmin, zmax, nz);
dZ = zg;

%-------------------------------------------------------------------------------
% Compute policy functions
%-------------------------------------------------------------------------------
% Get variable indices in DR structure
pf_names = {'c','l','k'};
pf_idx = zeros(length(pf_names),1);
for i = 1:length(pf_names)
    idx = find(strcmp(M_.endo_names, pf_names{i}));
    pf_idx(i) = find(oo_.dr.order_var == idx);
end

% Initialize arrays for policy functions
[cpf, lpf, kpf] = deal(zeros(nk, nz));

% Compute policy functions for each productivity level
for iz = 1:nz
    cpf(:, iz) = perturbed_pf(dK, dZ(iz), oo_.dr, pf_idx(1), options_.order);
    lpf(:,iz) = perturbed_pf(dK, dZ(iz), oo_.dr, pf_idx(2), options_.order);
    kpf(:,iz) = perturbed_pf(dK, dZ(iz), oo_.dr, pf_idx(3), options_.order);
end

%-------------------------------------------------------------------------------
% Plot policy functions
%-------------------------------------------------------------------------------
% Create figure
figure('Position', [100 100 600 800]);  % Adjusted for 3 rows, 1 column layout
tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'compact');

% Plot consumption policy
nexttile
hold on;
plot(kg, cpf, 'LineWidth', 1);
plot([kss kss], ylim, 'k--', 'LineWidth', 0.5);  % vertical line at steady state
hold off;
title('Consumption Policy Function');
xlabel('Capital Stock');
ylabel('Consumption');
grid on;

% Plot labor policy
nexttile
hold on;
plot(kg, lpf, 'LineWidth', 1);
plot([kss kss], ylim, 'k--', 'LineWidth', 0.5);  % vertical line at steady state
hold off;
title('Labor Policy Function');
xlabel('Capital Stock');
ylabel('Labor Supply');
grid on;

% Plot capital policy
nexttile
hold on;
plot(kg, kpf, 'LineWidth', 1);
plot(kg, kg, 'k--', 'LineWidth', 0.5);  % 45-degree line
plot([kss kss], ylim, 'k--', 'LineWidth', 0.5);  % vertical line at steady state
hold off;
title('Capital Policy Function');
xlabel('Capital Stock');
ylabel('Next Period Capital');
grid on;

% Add a legend showing productivity levels (for all subplots)
legend_entries = cell(1, nz);
for iz = 1:nz
    legend_entries{iz} = sprintf('z = %.2f', zg(iz));
end
legend(legend_entries, 'Location', 'southoutside');  % Positioned below the entire figure

% Overall title
sgtitle('Policy Functions: 3rd Order Perturbation', 'FontSize', 14);
%-------------------------------------------------------------------------------
% Create 3D surface plot for capital policy
%-------------------------------------------------------------------------------
figure('Position', [100 550 800 600]);

% Create meshgrid for surface plot
[K, Z] = meshgrid(kg, zg);

% Create surface plot
surf(K, Z, kpf', 'EdgeColor', 'none');
colormap('turbo');
colorbar;

% Add labels and title
xlabel('Capital Stock (k_t)');
ylabel('Productivity (z_t)');
zlabel('Next Period Capital (k_{t+1})');
title('Capital Policy Function (3D Surface)', 'FontSize', 14);

% Add steady state marker
hold on;
plot3(kss, 0, kss, 'k*', 'MarkerSize', 15);
text(kss, 0, kss, '  Steady State', 'VerticalAlignment', 'bottom');

% Add 45-degree reference line at z=0
iz_ss = find(abs(zg) == min(abs(zg)));  % Find index closest to z=0
plot3(kg, zeros(size(kg)), kg, 'k--', 'LineWidth', 2);

% Adjust view
view(45, 30);
grid on;
hold off;


%-------------------------------------------------------------------------------
% Helper function for policy function computation
%-------------------------------------------------------------------------------
function pf = perturbed_pf(dK, dZ, DR, i, approx_order)
    % Get steady states (in DR order)
    G0 = DR.ys(DR.order_var);
    
    % Get stochastic steady-state correction (if order > 1)
    if approx_order > 1
        GS2 = DR.ghs2(DR.order_var);
    else
        GS2 = zeros(size(G0));
    end
    
    % Get first-order coefficients
    G1 = DR.ghx;
    
    % Constants (including stochastic correction)
    constants = G0(i) + 0.5 * GS2(i);
    
    % First-order terms
    g_k = G1(i,1) * dK;      % Capital effect
    g_z = G1(i,2) * dZ;      % Productivity effect
    order1 = g_k + g_z;
    
    % Second-order terms
    order2 = 0;
    if approx_order > 1
        G2 = DR.ghxx;
        g_k2 = G2(i,1) * dK.^2;
        g_kz = 0.5*sum(G2(i,[2 3])) * dK * dZ;
        g_z2 = G2(i,4) * dZ.^2;
        order2 = 0.5 * (g_k2 + g_kz + g_z2);
    end
    
    % Third-order terms
    order3 = 0;
    if approx_order > 2
        G3 = DR.ghxxx;
        g_k3 = G3(i,1) * dK.^3;
        g_zk2 = (1/3) * sum(G3(i,[2 3 5])) * dK.^2 * dZ;
        g_kz2 = (1/3) * sum(G3(i,[4 6 7])) * dK * dZ.^2;
        g_z3 = G3(i,8) * dZ.^3;
        order3 = (1/6) * (g_k3 + g_zk2 + g_kz2 + g_z3);
    end
    
    % Combine all terms
    pf = constants + order1 + order2 + order3;
end