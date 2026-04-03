% =========================================================================
% PhaseDiagram_V_Phi.m
%
% Maps the insulator-metal (IM) phase boundary in the (V, Phi) plane for
% AB-stacked bilayer graphene with trigonal warping. For each value of the
% in-plane magnetic flux phi, the script finds the critical displacement
% field V_c(phi) at which the indirect band gap closes, using a bisection
% search over V. The result is the phase boundary curve that separates the
% gapped insulating phase (V > V_c) from the compensated semimetallic phase
% (V < V_c), as shown in Fig. 2(c) of the main text.
%
% Algorithm:
%   For each phi in phi_vals:
%     1. Precompute structure factors f(k), f(k+phi), f(k-phi) for all
%        k-points on the full grid (vectorized for speed).
%     2. Run a bisection search over V to locate the V_c where the indirect
%        gap E_g = min_k[E3(k) - E2(k)] changes sign.
%     3. Store V_c(phi) in Uc_vals.
%
% The indirect gap is defined as E_g = min_k(E3) - max_k(E2), where E2 and
% E3 are the two middle eigenvalues of the 4x4 Hamiltonian. E_g > 0
% signals an insulator; E_g <= 0 signals a compensated semimetal.
%
% Parameters to vary:
%   phi_vals     : range of in-plane magnetic flux values to sweep
%   v_min/max    : bisection search bounds for V (in eV)
%   N_k          : momentum grid resolution (trade-off: accuracy vs speed)
%
% Reference: Chakraborty et al., PRB (2025)
%
% Requirements: MATLAB Parallel Computing Toolbox (parfor loops)
% =========================================================================

clc; clear; close all;

%% --- Tight-Binding Parameters (eV) ---
gamma0 =  -2.610;   % Nearest-neighbor intralayer hopping
gamma1 =   0.361;   % Interlayer dimer hopping (B_t - A_b)
gamma3 =  -0.283;   % Skew interlayer hopping (trigonal warping)
gamma4 =   0.138;   % Skew interlayer hopping (electron-hole asymmetry)
Delta  =   0.015;   % On-site energy difference (dimer vs non-dimer sites)

%% --- Lattice Geometry ---
delta1 = [ 1,  sqrt(3)] / 2;
delta2 = [ 1, -sqrt(3)] / 2;
delta3 = [-1,  0      ];

%% --- Momentum Grid ---
% Full grid spanning the K-valley region; 500x500 captures the relevant
% low-energy pocket structure without excessive memory overhead
N_k     = 500;
kx_vals = linspace(-2*pi/3,              2*pi/3,           N_k);
ky_vals = linspace(-2*pi/(3*sqrt(3)),    2*pi/(3*sqrt(3)), N_k);
[kx_grid, ky_grid] = meshgrid(kx_vals, ky_vals);

Kx       = kx_grid(:);
Ky       = ky_grid(:);
N_points = numel(Kx);

% Structure factor f(k) = sum_n exp(-i k.delta_n)
calc_f = @(kx, ky) ...
    exp(-1i*(kx*delta1(1) + ky*delta1(2))) + ...
    exp(-1i*(kx*delta2(1) + ky*delta2(2))) + ...
    exp(-1i*(kx*delta3(1) + ky*delta3(2)));

%% --- Simulation Parameters ---
% Sweep range covers both low-field (trigonal-warped) and high-field regimes
phi_vals = linspace(0.02, 1.0, 100);   % In-plane magnetic flux values

% Bisection search bounds for the critical displacement field V_c
v_min_search = 0.0  / 1000;   % Lower bound (eV)
v_max_search = 1000 / 1000;   % Upper bound (eV)
tolerance    = 1e-6;           % Convergence criterion for bisection

fprintf('Starting phase boundary calculation on %d cores...\n', feature('numcores'));

%% --- Main Loop: Bisection Search for V_c at Each Phi (Parallelized) ---
Uc_vals = NaN(size(phi_vals));

parfor p_idx = 1:length(phi_vals)
    phi1 = phi_vals(p_idx);
    phi2 = 0;   % Zero y-flux: field along armchair direction

    % Precompute structure factors for all k-points at this flux value
    % (vectorized: avoids rebuilding inside the bisection loop)
    f_0  = calc_f(Kx, Ky);
    f_p  = calc_f(Kx + phi1, Ky + phi2);   % Top-layer shifted
    f_m  = calc_f(Kx - phi1, Ky - phi2);   % Bottom-layer shifted

    f_0c = conj(f_0);
    f_pc = conj(f_p);
    f_mc = conj(f_m);

    % --- Bisection over V to find V_c where E_g changes sign ---
    low  = v_min_search;
    high = v_max_search;

    % Safety check: if already insulating at the lowest V, record and skip
    gap_low = get_band_gap(low, Delta, gamma0, gamma1, gamma3, gamma4, ...
                           f_0, f_p, f_m, f_0c, f_pc, f_mc, N_points);
    if gap_low > 0
        Uc_vals(p_idx) = low;
    else
        while (high - low) > tolerance
            mid = (low + high) / 2;
            gap = get_band_gap(mid, Delta, gamma0, gamma1, gamma3, gamma4, ...
                               f_0, f_p, f_m, f_0c, f_pc, f_mc, N_points);
            if gap < 0
                low  = mid;   % Still semimetallic: increase V
            else
                high = mid;   % Already insulating: decrease V
            end
        end
        Uc_vals(p_idx) = high;
    end

    fprintf('  phi = %.5f  ->  V_c = %.5f eV  (%.4f meV)\n', ...
            phi1, Uc_vals(p_idx), Uc_vals(p_idx)*1000);
end

%% --- Figure: V_c vs Phi Phase Boundary ---
fig = figure('Position', [100, 100, 800, 600], 'Color', 'w');

plot(phi_vals, Uc_vals * 1000, 'o-', ...
    'LineWidth',      2, ...
    'MarkerSize',     8, ...
    'MarkerFaceColor','b', ...
    'MarkerEdgeColor','k', ...
    'Color',          'b');

xlabel('$\phi$',      'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
ylabel('$V_c$ (meV)', 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('IM Phase Boundary: Critical Displacement Field $V_c$ vs.\ In-Plane Flux $\phi$', ...
      'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');

grid on; grid minor; box on;
set(gca, 'FontSize', 14, 'FontName', 'Times New Roman', ...
         'LineWidth', 1.5, 'TickDir', 'out', 'GridAlpha', 0.15, ...
         'MinorGridAlpha', 0.05, 'Box', 'off');

%% --- Export Figure ---
exportgraphics(fig, 'PhaseDiagram_V_Phi.png', 'ContentType', 'vector', 'Resolution', 1200);
fprintf('Phase diagram saved as: PhaseDiagram_V_Phi.png\n');

%% =========================================================================
% Local Functions
% =========================================================================

function gap = get_band_gap(v, Delta, g0, g1, g3, g4, f0, fp, fm, f0c, fpc, fmc, N)
% GET_BAND_GAP  Computes the indirect band gap E_g = min_k(E3) - max_k(E2)
%               for the full 4x4 bilayer Hamiltonian at a given V and flux.
%
%   Inputs:
%     v     : displacement field (eV)
%     Delta : dimer/non-dimer on-site energy difference (eV)
%     g0,g1,g3,g4 : tight-binding hopping parameters (eV)
%     f0,fp,fm    : precomputed structure factors (N x 1 complex vectors)
%     f0c,fpc,fmc : conjugates of the above
%     N           : number of k-points
%
%   Output:
%     gap   : indirect gap in meV; positive = insulator, negative = semimetal

    min_cond =  inf;   % Tracks minimum of conduction band (E3)
    max_val  = -inf;   % Tracks maximum of valence band (E2)

    v2   = v / 2;
    vDp  =  v2 + Delta;   % Diagonal: A_bot site
    vDm  = -v2 + Delta;   % Diagonal: B_top site

    for k = 1:N
        % Assemble 4x4 Hamiltonian manually (avoids function-call overhead)
        % Basis: [A_top, B_top, A_bot, B_bot]
        H = [-v2,       g0*fp(k),  g4*f0(k),  g3*f0c(k);
             g0*fpc(k), vDm,       g1,         g4*f0(k);
             g4*f0c(k), g1,        vDp,        g0*fm(k);
             g3*f0(k),  g4*f0c(k), g0*fmc(k),  v2      ];

        E = sort(eig(H));   % Sorted eigenvalues: E(1) < E(2) < E(3) < E(4)

        if E(2) > max_val,  max_val  = E(2); end
        if E(3) < min_cond, min_cond = E(3); end
    end

    gap = (min_cond - max_val) * 1000;   % Convert eV -> meV
end
