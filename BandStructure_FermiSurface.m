% =========================================================================
% BandStructure_FermiSurface.m
%
% Computes the low-energy band structure of AB-stacked (Bernal) bilayer
% graphene including trigonal warping via the Slonczewski-Weiss-McClure
% tight-binding model. The script:
%   1. Builds the full 4x4 Hamiltonian at each k-point on a fine grid
%      centered on the K valley, incorporating skew interlayer couplings
%      (gamma3, gamma4) and an in-plane magnetic field via the Peierls
%      substitution (controlled by phi1, phi2).
%   2. Estimates the Fermi energy by equating the enclosed k-space areas
%      of the conduction and valence Fermi pockets (charge neutrality).
%   3. Plots 2D energy contours of both bands at the Fermi level.
%   4. Plots a publication-quality 3D surface of the two middle bands
%      with the Fermi plane and highlighted Fermi contours.
%
% Parameters to vary:
%   v    : displacement field (interlayer potential asymmetry) in eV
%   phi1 : Peierls flux (x-component), proportional to B_parallel * a * l
%   phi2 : Peierls flux (y-component), set to 0 for armchair field direction
%
% Reference: Chakraborty et al., PRB (2025); Seiler et al., Nat. Commun. (2024)
%
% Requirements: MATLAB Parallel Computing Toolbox (parfor loops)
% =========================================================================

clc; clear; close all;

%% --- Tight-Binding Parameters (eV) ---
% Values from Kuzmenko et al. PRB 2009 / Seiler et al. Nat. Commun. 2024
gamma0 =  -2.610;   % Nearest-neighbor intralayer hopping
gamma1 =   0.361;   % Interlayer dimer hopping (B_t - A_b)
gamma3 =  -0.283;   % Skew interlayer hopping (trigonal warping)
gamma4 =   0.138;   % Skew interlayer hopping (electron-hole asymmetry)
Delta  =   0.015;   % On-site energy difference (dimer vs non-dimer sites)

%% --- Lattice Geometry (units of carbon bond length a = 1) ---
% Nearest-neighbor vectors from A to B within a layer
delta1 = [ 1,  sqrt(3)] / 2;
delta2 = [ 1, -sqrt(3)] / 2;
delta3 = [-1,  0      ];

% K and K' valley positions in the Brillouin zone
kd1 = [2*pi/3,  2*pi/(3*sqrt(3))];   % K  valley
kd2 = [2*pi/3, -2*pi/(3*sqrt(3))];   % K' valley

%% --- Momentum Grid (centered on K valley) ---
% Resolution of 2000x2000 captures trigonal-warping fine structure
% (~1/100 of the full BZ, within ~1 meV of the Dirac point)
kx_vals = linspace(kd1(1) - 0.015, kd1(1) + 0.015, 2000);
ky_vals = linspace(kd1(2) - 0.015, kd1(2) + 0.015, 2000);
[kx_grid, ky_grid] = meshgrid(kx_vals, ky_vals);

% Flatten for parallelized k-loop
kx_list = kx_grid(:);
ky_list = ky_grid(:);
num_k   = numel(kx_list);

%% --- External Field Parameters ---
v    = 0.0000;   % Displacement field (interlayer potential asymmetry) in eV
                 % Set v > V_c ~ 0.00063 eV to open a bulk gap
phi1 = 0.0054;   % In-plane magnetic flux (x); phi = B*a*l/Phi_0
phi2 = 0.0000;   % In-plane magnetic flux (y); 0 for armchair field direction

%% --- Function Handles ---
% Lattice structure factor f(k) = sum_{n} exp(-i k.delta_n)
f_scalar = @(kx, ky) ...
    exp(-1i*(kx*delta1(1) + ky*delta1(2))) + ...
    exp(-1i*(kx*delta2(1) + ky*delta2(2))) + ...
    exp(-1i*(kx*delta3(1) + ky*delta3(2)));

% Full 4x4 Hamiltonian; Peierls substitution shifts k -> k +/- phi per layer
% Basis order: [A_top, B_top, A_bot, B_bot]
Hamiltonian = @(kx, ky, v, phi1, phi2) ...
    [-v/2,                            gamma0 * f_scalar(kx+phi1, ky+phi2),  gamma4 * f_scalar(kx, ky),         gamma3 * conj(f_scalar(kx, ky));
     gamma0 * conj(f_scalar(kx+phi1, ky+phi2)), -v/2 + Delta,              gamma1,                             gamma4 * f_scalar(kx, ky);
     gamma4 * conj(f_scalar(kx, ky)),            gamma1,                   v/2 + Delta,                        gamma0 * f_scalar(kx-phi1, ky-phi2);
     gamma3 * f_scalar(kx, ky),                  gamma4*conj(f_scalar(kx,ky)), gamma0*conj(f_scalar(kx-phi1,ky-phi2)), v/2];

% Returns the four sorted eigenvalues at a given k-point
Bands = @(kx, ky, v, phi1, phi2) sort(eig(Hamiltonian(kx, ky, v, phi1, phi2)));

%% --- Band Eigenvalue Calculation (Parallelized) ---
E_cond = zeros(num_k, 1);   % 3rd eigenvalue: bottom of conduction band
E_val  = zeros(num_k, 1);   % 2nd eigenvalue: top of valence band

parfor idx = 1:num_k
    bands       = Bands(kx_list(idx), ky_list(idx), v, phi1, phi2);
    E_cond(idx) = bands(3) * 1000;   % Convert eV -> meV
    E_val(idx)  = bands(2) * 1000;
end

% Reshape results back to grid
E_cond = reshape(E_cond, size(kx_grid));
E_val  = reshape(E_val,  size(kx_grid));

%% --- Fermi Energy Estimation (Charge Neutrality via Pocket Areas) ---
% At charge neutrality the total electron and hole pocket areas in k-space
% must be equal. We scan a range of trial energies and find the crossing.

E_search  = linspace(0.0, 1.5, 1000);   % Trial energy range in meV
area_cond = zeros(size(E_search));
area_val  = zeros(size(E_search));

% Contour area helper: sums polyarea over all disconnected contour loops
computeContourArea = @(C) computeArea(C);

parfor k = 1:length(E_search)
    E_t = E_search(k);
    % Conduction pocket area at trial energy
    C_c          = contourc(kx_vals, ky_vals, E_val,  [E_t, E_t]);
    area_cond(k) = computeArea(C_c);
    % Valence pocket area at trial energy (three satellite pockets)
    C_v          = contourc(kx_vals, ky_vals, E_cond, [E_t, E_t]);
    area_val(k)  = computeArea(C_v);
end

% Interpolate to find exact crossing (E_F where areas are equal)
diff_area   = area_cond - area_val;
E_F         = interp1(diff_area, E_search, 0, 'linear', NaN);

fprintf('Estimated Fermi level: E_F = %.4f meV  (at V = %.2f meV, phi1 = %.4f)\n', ...
        E_F, v*1000, phi1);

%% --- Figure 1: 2D Fermi Surface Contour Plot ---
figure('Color', 'w');

% Energy levels bracketing E_F for a narrow contour window
E_levels = linspace(E_F - 0.002, E_F + 0.002, 20);

contour(kx_grid, ky_grid, E_cond, E_levels, '-', 'LineWidth', 2.0, 'LineColor', 'b');
hold on
contour(kx_grid, ky_grid, E_val,  E_levels, '-', 'LineWidth', 2.0, 'LineColor', 'r');
hold off

xlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Fermi surface contours:  $E_F$ = %.4f meV,  $V$ = %.2f meV', ...
              E_F, v*1000), 'Interpreter', 'latex', 'FontSize', 14);

xlim([kd1(1) - 0.01, kd1(1) + 0.01]);
ylim([kd1(2) - 0.015, kd1(2) + 0.015]);
axis equal
set(gca, 'FontSize', 14, 'LineWidth', 1.5);

%% --- Figure 2: 3D Band Surface with Fermi Plane ---
% Color shading: conduction band in blue-to-white, valence in red-to-white

norm1   = (E_cond + 2) / 6;
norm2   = (E_val  + 2) / 6;
blue    = [0, 0, 1];
red     = [1, 0, 0];
color1  = (1 - norm1) .* reshape(blue, 1,1,3) + norm1 .* ones(size(E_cond,1), size(E_cond,2), 3);
color2  = (1 - norm2) .* reshape(red,  1,1,3) + norm2 .* ones(size(E_val,1),  size(E_val,2),  3);

figure('Color', 'w', 'Position', [100, 100, 700, 600]);

% Conduction band surface
h1       = surf(kx_grid, ky_grid, E_cond, 'EdgeColor', 'none', 'FaceColor', 'interp');
h1.CData = color1;
hold on

% Valence band surface
h2       = surf(kx_grid, ky_grid, E_val, 'EdgeColor', 'none', 'FaceColor', 'interp');
h2.CData = color2;

% Semi-transparent Fermi plane
surf(kx_grid, ky_grid, E_F * ones(size(kx_grid)), ...
     'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Fermi contours highlighted on both surfaces
contour3(kx_grid, ky_grid, E_cond, [E_F, E_F], 'LineColor', blue, 'LineWidth', 5.0);
contour3(kx_grid, ky_grid, E_val,  [E_F, E_F], 'LineColor', red,  'LineWidth', 5.0);

view(-20, 20);
zlim([-1.5, 2.5]);
pbaspect([2 2 1]);
set(gca, 'FontSize', 20, 'FontWeight', 'bold', 'LineWidth', 2, ...
         'Box', 'on', 'TickDir', 'out', 'TickLabelInterpreter', 'latex');

%% --- Export Figure ---
filename = sprintf('BGTW_Contour_V%0.2fmeV_phi1_%0.4f.png', v*1000, phi1);
print(filename, '-dpng', '-r1200');
fprintf('Figure saved as: %s\n', filename);

%% =========================================================================
% Local Functions
% =========================================================================

function total_area = computeArea(C)
% COMPUTEAREA  Sums the polygon areas of all contour loops returned by
%              contourc. Handles multiple disconnected loops correctly.
    total_area = 0;
    i = 1;
    while i < size(C, 2)
        n   = C(2, i);                      % Number of vertices in this loop
        x   = C(1, i+1 : i+n);
        y   = C(2, i+1 : i+n);
        total_area = total_area + polyarea(x, y);
        i   = i + n + 1;
    end
end
