% ========================================================================
% Exercise 3.5 — Winnowing Process
% Simulations of Mechanical Processes, OVGU Magdeburg, SoSe 2025
% Author: Jaideep Buksagar
% ========================================================================

clear all
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor',   'w');

% =============================== PARAMETERS =============================
g = 9.81;

% Air properties
rho_f = 1.2;
mu_f  = 1.8e-5;

% Jet properties
uf0 = 0.20;
h   = 0.10;

% Release point
x0 = 0.50;  y0 = 0.50;
u0 = 0;     v0 = 0;

% Bin separator (x_c, y_c)
x_bin =  0.55;
y_bin = -0.50;

% Particles
grain.rho_p = 750;     grain.dp = 2.5e-3;
chaff.rho_p =  50;     chaff.dp = 3.25e-3;

% Time settings
dt    = 0.001;
t_max = 5;

% ============================== TASK 1 ============================================================================================================

% Particle trajectories with Euler and RK4

grain_euler = eulerSolver(grain, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt,t_max, x_bin,y_bin, g);
chaff_euler = eulerSolver(chaff, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt,t_max, x_bin,y_bin, g);

grain_rk4   = rk4Solver(grain, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt,t_max, x_bin,y_bin, g);
chaff_rk4   = rk4Solver(chaff, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt,t_max, x_bin,y_bin, g);

fprintf('\n=== TASK 1: Landing positions (RK4) ===\n');
fprintf('Grain final x: %.4f m  -> Bin %d\n', grain_rk4.x(end), classifyBin(grain_rk4.x(end), x_bin));
fprintf('Chaff final x: %.4f m  -> Bin %d\n', chaff_rk4.x(end), classifyBin(chaff_rk4.x(end), x_bin));


% ---- Figure 1: Euler trajectories ------------------------------------------------------------------------------
figure(1); clf;
plot(grain_euler.x, grain_euler.y, 'b-', 'LineWidth', 2); hold on;
plot(chaff_euler.x, chaff_euler.y, 'r-', 'LineWidth', 2);

y_range = [min([grain_euler.y, chaff_euler.y, y_bin]) - 0.05, y0 + 0.05];
x_range = [0, max([grain_euler.x, chaff_euler.x]) + 0.05];

plot([x_bin x_bin], y_range, 'k--', 'LineWidth', 1);
plot(x_range, [y_bin y_bin], 'k:',  'LineWidth', 1);
hold off;

xlabel('x Position (m)'); ylabel('y Position (m)');
title('Particle Trajectories: Euler');
legend('Grain', 'Chaff', 'Bin boundary', 'Ground', 'Location', 'best');
grid on; axis equal;


% ---- Figure 2: RK4 trajectories --------------------------------------------------------------------------------
figure(2); clf;
plot(grain_rk4.x, grain_rk4.y, 'b-', 'LineWidth', 2); hold on;
plot(chaff_rk4.x, chaff_rk4.y, 'r-', 'LineWidth', 2);

y_range = [min([grain_rk4.y, chaff_rk4.y, y_bin]) - 0.05, y0 + 0.05];
x_range = [0, max([grain_rk4.x, chaff_rk4.x]) + 0.05];

plot([x_bin x_bin], y_range, 'k--', 'LineWidth', 1);
plot(x_range, [y_bin y_bin], 'k:',  'LineWidth', 1);
hold off;

xlabel('x Position (m)'); ylabel('y Position (m)');
title('Particle Trajectories: RK4');
legend('Grain', 'Chaff', 'Bin boundary', 'Ground', 'Location', 'best');
grid on; axis equal;

% ============================== TASK 2 =============================================================================================================

% Convergence verification using fixed-time error metric

dt_ref = 1e-6;
ref    = rk4Solver(grain, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt_ref, t_max, x_bin, y_bin, g);
t_check = 0.10;   % well before any particle lands

dt_list = [1e-2, 5e-3, 1e-3, 5e-4, 1e-4];

err_euler = zeros(size(dt_list));
err_rk4   = zeros(size(dt_list));

x_ref_check = interp1(ref.t, ref.x, t_check, 'linear');

for k = 1:length(dt_list)
    dt_k = dt_list(k);
    
    sol_e = eulerSolver(grain, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt_k, t_max, x_bin, y_bin, g);
    sol_r = rk4Solver  (grain, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt_k, t_max, x_bin, y_bin, g);
    
    x_euler_check = interp1(sol_e.t, sol_e.x, t_check, 'linear');
    x_rk4_check   = interp1(sol_r.t, sol_r.x, t_check, 'linear');
    
    err_euler(k) = abs(x_euler_check - x_ref_check);
    err_rk4(k)   = abs(x_rk4_check   - x_ref_check);
end

% Measured slopes (least-squares fit on log-log)
p_euler = polyfit(log(dt_list), log(err_euler), 1);
p_rk4   = polyfit(log(dt_list), log(err_rk4),   1);

fprintf('\n=== TASK 2: Measured convergence orders ===\n');
fprintf('Euler slope: %.3f (expected 1)\n', p_euler(1));
fprintf('RK4   slope: %.3f (expected 4)\n', p_rk4(1));

% Local (per-step) slopes for diagnostic table
fprintf('\n%-10s %-15s %-15s %-12s %-12s\n','dt','err_euler','err_rk4','slope_E','slope_R');
for k = 2:length(dt_list)
    s_e = log(err_euler(k)/err_euler(k-1)) / log(dt_list(k)/dt_list(k-1));
    s_r = log(err_rk4(k)/err_rk4(k-1))     / log(dt_list(k)/dt_list(k-1));
    fprintf('%-10.1e %-15.3e %-15.3e %-12.2f %-12.2f\n', dt_list(k), err_euler(k), err_rk4(k), s_e, s_r);
end

% ---- Figure 3: Convergence study -------------------------------------------------------------------------------------
figure(3); clf;
loglog(dt_list, err_euler, 'b-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog(dt_list, err_rk4,   'r-s', 'LineWidth', 2, 'MarkerSize', 8);

% Reference slope lines (anchored at largest dt)
dt_plot = [dt_list(end), dt_list(1)];
loglog(dt_plot, err_euler(1) * (dt_plot/dt_list(1)).^1, 'b--', 'LineWidth', 1.2);
loglog(dt_plot, err_rk4(1)   * (dt_plot/dt_list(1)).^4, 'r--', 'LineWidth', 1.2);
hold off;

set(gca, 'XDir', 'reverse');   % conventional: large dt on left
xlabel('\Delta t (s)'); ylabel('Error in x-position at t* (m)');
title(sprintf('Convergence Study (t* = %.2f s): Euler order %.2f, RK4 order %.2f', t_check, p_euler(1), p_rk4(1)));
legend('Euler','RK4','Slope 1 reference','Slope 4 reference','Location','southwest');
grid on;

% ============================== TASK 3 ===========================================================================================================
% Minimum jet velocity for separation via bisection

precision   = 1e-6;
uf0_low     = 0.00;
uf0_high    = 0.20;
dt_task3    = 1e-4;
max_iter    = 100;
uf0_history = zeros(1, max_iter);

iter = 0;
while (uf0_high - uf0_low) > precision && iter < max_iter
    iter = iter + 1;
    uf0_mid = (uf0_low + uf0_high) / 2;
    x_land_chaff = landingX(uf0_mid, chaff, rho_f,mu_f, h, x0,y0, u0,v0, dt_task3, t_max, x_bin, y_bin, g);
    uf0_history(iter) = uf0_mid;
    
    if x_land_chaff - x_bin < 0
        uf0_low  = uf0_mid;
    else
        uf0_high = uf0_mid;
    end
end

uf0_min = (uf0_low + uf0_high) / 2;

fprintf('\n=== TASK 3: Minimum jet velocity ===\n');
fprintf('uf0_min = %.6f m/s = %.4f cm/s\n', uf0_min, uf0_min*100);
fprintf('Converged in %d iterations to tolerance %.1e m/s\n', iter, precision);

% ---- Figure 4: Bisection convergence ----
figure(4); clf;
valid_iters = uf0_history(1:iter);
semilogy(1:iter, abs(valid_iters - uf0_min), 'k-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('Iteration'); ylabel('|u_0^{(k)} - u_0^*|  (m/s)');
title('Bisection Convergence: Minimum Jet Velocity');
grid on;

% ============================== TASK 4 ==========================================================================================================
% Non-zero initial velocity: |v_init| = v_terminal, angle in [-108, -72] deg

vt_grain = terminalVelocity(grain, rho_f, mu_f, g);
vt_chaff = terminalVelocity(chaff, rho_f, mu_f, g);

fprintf('\n=== TASK 4: Terminal velocities ===\n');
fprintf('Grain v_t: %.4f m/s\n', vt_grain);
fprintf('Chaff v_t: %.4f m/s\n', vt_chaff);

% --- Compute the four extreme-angle trajectories at default uf0 = 0.20 ------------------------------
% Angle convention: theta measured from +x axis (so -90 is straight down)
% Range: theta in [-108, -72]. Endpoints define worst cases by continuity.

angles_deg = [-72, -108];

[ext_default, paths_default] = extremeAnglesLanding(grain, chaff, vt_grain, vt_chaff, angles_deg, rho_f,mu_f, uf0,h, x0,y0, dt_task3, t_max, x_bin,y_bin, g);

fprintf('\n=== TASK 4: Extreme-angle landing positions at uf0 = %.2f m/s ===\n', uf0);
fprintf('%-8s %-8s %-12s\n', 'Particle','Angle','x_land (m)');
fprintf('%-8s %-8d %-12.4f\n', 'Grain', -72,  ext_default.grain_m72);
fprintf('%-8s %-8d %-12.4f\n', 'Grain', -108, ext_default.grain_m108);
fprintf('%-8s %-8d %-12.4f\n', 'Chaff', -72,  ext_default.chaff_m72);
fprintf('%-8s %-8d %-12.4f\n', 'Chaff', -108, ext_default.chaff_m108);

% ---- Figure 5: Task 4 trajectories at default uf0 -----------------------------------------------------------------
figure(5); clf;
plot(paths_default.grain_m72.x,  paths_default.grain_m72.y,  'b-',  'LineWidth', 1.5); hold on;
plot(paths_default.grain_m108.x, paths_default.grain_m108.y, 'b--', 'LineWidth', 1.5);
plot(paths_default.chaff_m72.x,  paths_default.chaff_m72.y,  'r-',  'LineWidth', 1.5);
plot(paths_default.chaff_m108.x, paths_default.chaff_m108.y, 'r--', 'LineWidth', 1.5);
plot([x_bin x_bin], [-0.6 0.6], 'k--', 'LineWidth', 1);
plot([0 1.5], [y_bin y_bin], 'k:', 'LineWidth', 1);
hold off;
xlabel('x (m)'); ylabel('y (m)');
title('Task 4: Trajectories at extreme angles');
legend('Grain -72°','Grain -108°','Chaff -72°','Chaff -108°','Bin boundary','Ground','Location','best');
grid on; axis equal;


% ---------------------------------------------------------------------------------------------------------------------------------------------------------
% TASK 4 ANALYSIS: ranges of landing positions at default uf0
% ---------------------------------------------------------------------------------------------------------------------------------------------------------
% By continuity in the initial-velocity angle, the grain,s landing x-range across theta in [-108,-72] lies between min and max of its two extremes.
% Same for chaff. Perfect separation requires the grain-range and chaff-range to be disjoint AND a single x_c to lie between them.

grain_xmin = min(ext_default.grain_m72, ext_default.grain_m108);
grain_xmax = max(ext_default.grain_m72, ext_default.grain_m108);
chaff_xmin = min(ext_default.chaff_m72, ext_default.chaff_m108);
chaff_xmax = max(ext_default.chaff_m72, ext_default.chaff_m108);

fprintf('\n=== TASK 4 ranges at uf0 = %.2f m/s ===\n', uf0);
fprintf('Grain x-range: [%.4f, %.4f] m\n', grain_xmin, grain_xmax);
fprintf('Chaff x-range: [%.4f, %.4f] m\n', chaff_xmin, chaff_xmax);

% ---- (a) Can we achieve perfect separation by changing x_c only? ----------------------------------------------------------------------------------------
fprintf('\n--- Task 4(a): vary x_c only ---\n');
if grain_xmax < chaff_xmin
    fprintf('YES: choose x_c in (%.4f, %.4f) m\n', grain_xmax, chaff_xmin);
elseif chaff_xmax < grain_xmin
    fprintf('YES: choose x_c in (%.4f, %.4f) m (chaffs left of grains)\n', chaff_xmax, grain_xmin);
else
    fprintf('NO: grain and chaff ranges overlap at uf0 = %.2f m/s\n', uf0);
end


% ---- (b) Can we achieve perfect separation by changing y_c only? ----------------------------------------------------------------------------------------------------------
% Key insight: changing y_c means classifying at a different vertical level.
% We need to check whether at SOME y (not just y = -0.5), the grain x-coordinate range is fully separated from the chaff x-coordinate range across the entire angle interval.
fprintf('\n--- Task 4(b): vary y_c only ---\n');

% Sample y_c values from y0 down to a deep level
yc_samples = y0 - 0.01 : -0.01 : -0.99;
yc_works   = [];

for j = 1:length(yc_samples)
    
    yc_try = yc_samples(j);
    
    % Recompute extreme landings at this y_c
    [ext_y, ~] = extremeAnglesLanding(grain, chaff, vt_grain, vt_chaff, angles_deg, rho_f,mu_f, uf0,h, x0,y0, dt_task3, t_max, x_bin, yc_try, g);
    
    gmin = min(ext_y.grain_m72, ext_y.grain_m108);
    gmax = max(ext_y.grain_m72, ext_y.grain_m108);
    cmin = min(ext_y.chaff_m72, ext_y.chaff_m108);
    cmax = max(ext_y.chaff_m72, ext_y.chaff_m108);
    
    if gmax < cmin || cmax < gmin
        yc_works = [yc_works, yc_try]; %#ok<AGROW>
    end
end


if isempty(yc_works)
    fprintf('NO: no y_c in tested range gives disjoint grain/chaff x-ranges\n');
else
    fprintf('YES: y_c values that work range from %.3f to %.3f m\n', min(yc_works), max(yc_works));
end

% ---- (c) Can we achieve perfect separation by changing u0 only? --------------------------------------------------------------------------------------------------------------------------------
% u0 constrained to [0.05, 0.40] m/s. Look for u0 where grain and chaff x-ranges are disjoint and the existing x_bin = 0.55 lies between them (since we are NOT allowed to change x_c here).
fprintf('\n--- Task 4(c): vary u_0 only, u_0 in [0.05, 0.40] m/s ---\n');

uf0_samples = 0.05:0.01:0.40;
uf0_works   = [];

for j = 1:length(uf0_samples)
    
    uf0_try = uf0_samples(j);
    
    [ext_u, ~] = extremeAnglesLanding(grain, chaff, vt_grain, vt_chaff, angles_deg, rho_f,mu_f, uf0_try,h, x0,y0, dt_task3, t_max, x_bin,y_bin, g);
    
    gmin = min(ext_u.grain_m72, ext_u.grain_m108);
    gmax = max(ext_u.grain_m72, ext_u.grain_m108);
    cmin = min(ext_u.chaff_m72, ext_u.chaff_m108);
    cmax = max(ext_u.chaff_m72, ext_u.chaff_m108);
    
    % Disjoint AND x_bin = 0.55 in the gap
    if gmax < cmin && gmax < x_bin && x_bin < cmin
        uf0_works = [uf0_works, uf0_try]; %#ok<AGROW>
    end
end


if isempty(uf0_works)
    fprintf('NO: no u_0 in [0.05, 0.40] m/s separates with the fixed x_c = %.2f m\n', x_bin);
else
    fprintf('YES: u_0 values that work range from %.3f to %.3f m/s\n', min(uf0_works), max(uf0_works));
end


% ---- (d) Can we achieve perfect separation by changing x_c AND y_c? -----------------------------------------------------------------------------------
% More degrees of freedom than (b). We re-use (b)s logic but also allow picking ANY x_c in the gap between grain and chaff ranges at that y_c.
% This is feasible whenever there exists a y_c with disjoint ranges.


fprintf('\n--- Task 4(d): vary x_c and y_c ---\n');

if ~isempty(yc_works)
    fprintf('YES: same y_c range as (b) works, with x_c chosen in the gap at each y_c\n');
    fprintf('     Feasible y_c range: %.3f to %.3f m\n', min(yc_works), max(yc_works));
else
    % Need to check more carefully: even if (b) fails at uf0=0.20, perhaps at some y_c the grain and chaff x-ranges become disjoint.
    % Our yc_samples already swept this. So same answer as (b).
    fprintf('NO: no (x_c, y_c) combination achieves perfect separation at uf0 = %.2f m/s\n', uf0);
end

% ---- (e) Can we achieve perfect separation by changing x_c, y_c, AND u_0? ----
% Maximum degrees of freedom. Sweep u_0 in [0.05, 0.40] and for each,
% check if any y_c gives disjoint ranges (then x_c can be placed in the gap).
fprintf('\n--- Task 4(e): vary x_c, y_c, and u_0 ---\n');

uf0_e_works = [];
for j = 1:length(uf0_samples)
    uf0_try = uf0_samples(j);
    
    % Check if any y_c gives disjoint ranges at this u_0
    for k = 1:length(yc_samples)
        yc_try = yc_samples(k);
        [ext_e, ~] = extremeAnglesLanding(grain, chaff, vt_grain, vt_chaff, angles_deg, rho_f,mu_f, uf0_try,h, x0,y0, dt_task3, t_max, x_bin, yc_try, g);
        
        gmin = min(ext_e.grain_m72, ext_e.grain_m108);
        gmax = max(ext_e.grain_m72, ext_e.grain_m108);
        cmin = min(ext_e.chaff_m72, ext_e.chaff_m108);
        cmax = max(ext_e.chaff_m72, ext_e.chaff_m108);
        
        if gmax < cmin || cmax < gmin
            uf0_e_works = [uf0_e_works, uf0_try]; %#ok<AGROW>
            break  % move to next u_0 once we have found a working y_c
        end
    end
end

if isempty(uf0_e_works)
    fprintf('NO: no (x_c, y_c, u_0) combination achieves perfect separation\n');
else
    fprintf('YES: u_0 values with feasible (x_c, y_c) range from %.3f to %.3f m/s\n', ...
        min(uf0_e_works), max(uf0_e_works));
end

fprintf('\n========== ALL TASKS COMPLETE ==========\n');

% ==============================================================================================================================================================================
% ========================= FUNCTIONS ==========================================================================================================================================
% ==============================================================================================================================================================================

% -----------------------------------------------------------------------------------Euler Solver-------------------------------------------------------------------------------

function sol = eulerSolver(particle, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt,t_max, x_bin,y_bin, g) %#ok<INUSD>
    rho_p = particle.rho_p;
    dp    = particle.dp;
    time_steps = round(t_max/dt);
    
    Y = zeros(4, time_steps);
    Y(:,1) = [x0; y0; u0; v0];
    t = zeros(1, time_steps);
    
    final_idx = time_steps;
    for i = 1:time_steps-1
        RHS = winnowingODEsRHS(Y(:,i), rho_f,mu_f, rho_p,dp, uf0,h, g);
        Y(:,i+1) = Y(:,i) + dt*RHS;
        t(i+1)   = t(i) + dt;
        
        if Y(2,i+1) <= y_bin
            % linear interpolation to exact landing
            frac = (Y(2,i) - y_bin) / (Y(2,i) - Y(2,i+1));
            Y(:,i+1) = Y(:,i) + frac * (Y(:,i+1) - Y(:,i));
            t(i+1)   = t(i) + frac * dt;
            final_idx = i+1;
            break
        end
    end
    
    sol.x = Y(1, 1:final_idx);
    sol.y = Y(2, 1:final_idx);
    sol.u = Y(3, 1:final_idx);
    sol.v = Y(4, 1:final_idx);
    sol.t = t(1:final_idx);
end

%---------------------------------------------------------------------------------------RK4 Solver-----------------------------------------------------------------------------

function sol = rk4Solver(particle, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt,t_max, x_bin,y_bin, g) %#ok<INUSD>
    rho_p = particle.rho_p;
    dp    = particle.dp;
    time_steps = round(t_max/dt);
    
    Y = zeros(4, time_steps);
    Y(:,1) = [x0; y0; u0; v0];
    t = zeros(1, time_steps);
    
    final_idx = time_steps;
    for i = 1:time_steps-1
        k1 = winnowingODEsRHS(Y(:,i),               rho_f,mu_f, rho_p,dp, uf0,h, g);
        k2 = winnowingODEsRHS(Y(:,i) + dt*k1/2,     rho_f,mu_f, rho_p,dp, uf0,h, g);
        k3 = winnowingODEsRHS(Y(:,i) + dt*k2/2,     rho_f,mu_f, rho_p,dp, uf0,h, g);
        k4 = winnowingODEsRHS(Y(:,i) + dt*k3,       rho_f,mu_f, rho_p,dp, uf0,h, g);
        
        Y(:,i+1) = Y(:,i) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        t(i+1)   = t(i) + dt;
        
        if Y(2,i+1) <= y_bin
            frac = (Y(2,i) - y_bin) / (Y(2,i) - Y(2,i+1));
            Y(:,i+1) = Y(:,i) + frac * (Y(:,i+1) - Y(:,i));
            t(i+1)   = t(i) + frac * dt;
            final_idx = i+1;
            break
        end
    end
    
    sol.x = Y(1, 1:final_idx);
    sol.y = Y(2, 1:final_idx);
    sol.u = Y(3, 1:final_idx);
    sol.v = Y(4, 1:final_idx);
    sol.t = t(1:final_idx);
end

%-------------------------------------------------------------------------------------RHS Caluculator-------------------------------------------------------------------

function RHS = winnowingODEsRHS(Y, rho_f,mu_f, rho_p,dp, uf0,h, g)
    x = Y(1); y = Y(2);
    u = Y(3); v = Y(4);
    
    Vp = (dp^3) * pi/6;
    mp = rho_p * Vp;
    
    % Fluid velocity (planar jet)
    if x > 1e-9
        uf_x = 6.2 * uf0 * sqrt(h/x) * exp(-50 * (y^2)/(x^2));
    else
        uf_x = 0;
    end
    uf_y = 0;
    
    ur_x   = uf_x - u;
    ur_y   = uf_y - v;
    ur_mag = sqrt(ur_x^2 + ur_y^2);
    
    Re = rho_f * ur_mag * dp / mu_f;
    
    if Re < 1e-12
        Cd = 0;
    elseif Re < 800
        Cd = (24/Re) * (1 + 0.15*Re^0.687);
    else
        Cd = 0.44;
    end
    
    Fd_x = 0.5 * pi * dp^2 * rho_f * Cd * ur_mag * ur_x;
    Fd_y = 0.5 * pi * dp^2 * rho_f * Cd * ur_mag * ur_y;
    
    Fg = -rho_p * Vp * g;
    Fb =  rho_f * Vp * g;
    
    ax = Fd_x / mp;
    ay = (Fg + Fb + Fd_y) / mp;
    
    RHS = [u; v; ax; ay];
end

%-------------------------------------------------------------------------------Landing X coordinate calculator--------------------------------------------------------------

function x_land = landingX(uf0, particle, rho_f,mu_f, h, x0,y0, u0,v0, dt, t_max, x_bin, y_bin, g)
    sol = rk4Solver(particle, rho_f,mu_f, uf0,h, x0,y0, u0,v0, dt, t_max, x_bin, y_bin, g);
    x_land = sol.x(end);
end

%----------------------------------------------------------------------------------Terminal velocity caalculator-------------------------------------------------------------

function vt = terminalVelocity(particle, rho_f, mu_f, g)
    rho_p = particle.rho_p;
    dp    = particle.dp;
    Vp    = pi * dp^3 / 6;
    
    F_net = (rho_p - rho_f) * Vp * g;
    
    vt = F_net / (3 * pi * mu_f * dp);   % Stokes initial guess
    
    Precision = 1e-10;
    max_iter  = 50;
    
    for k = 1:max_iter
        Re = rho_f * vt * dp / mu_f;
        if Re < 800
            Cd      = (24/Re) * (1 + 0.15*Re^0.687);
            dCd_dRe = -(24/Re^2)*(1 + 0.15*Re^0.687) + (24/Re)*(0.15*0.687*Re^(0.687-1));
        else
            Cd      = 0.44;
            dCd_dRe = 0;
        end
        dRe_dvt   = rho_f * dp / mu_f;
        drag      = 0.5 * pi * dp^2 * rho_f * Cd * vt^2;
        f         = drag - F_net;
        ddrag_dvt = 0.5 * pi * dp^2 * rho_f * (dCd_dRe * dRe_dvt * vt^2 + Cd * 2 * vt);
        
        vt_new = vt - f/ddrag_dvt;
        if abs(vt_new - vt) < Precision
            vt = vt_new;
            return
        end
        vt = vt_new;
    end
end

%--------------------------------------------------------------------------------------Bin classifier-----------------------------------------------------------------------

function bin = classifyBin(particle_x, x_bin)
    if particle_x < x_bin
        bin = 1;
    else
        bin = 2;
    end
end

%-----------------------------------------------------------------------------------Extreme angles trajectory---------------------------------------------------------------

function [extremes, paths] = extremeAnglesLanding(grain, chaff, vt_grain, vt_chaff, ...
    angles_deg, rho_f,mu_f, uf0,h, x0,y0, dt, t_max, x_bin,y_bin, g)
% Compute landing x-positions for the 4 (particle, angle) combinations
% used in Task 4. Returns a struct of 4 named landing positions and a
% struct of 4 named full trajectories.

    angles_rad = deg2rad(angles_deg);
    
    % theta_1 = -72 (mostly down, slight +x), theta_2 = -108 (mostly down, slight -x)
    th_72  = angles_rad(angles_deg == -72);
    th_108 = angles_rad(angles_deg == -108);
    
    % --- Grain, -72 degrees ---
    u_init = vt_grain * cos(th_72);
    v_init = vt_grain * sin(th_72);
    sol = rk4Solver(grain, rho_f,mu_f, uf0,h, x0,y0, u_init,v_init, dt, t_max, x_bin,y_bin, g);
    extremes.grain_m72 = sol.x(end);
    paths.grain_m72    = sol;
    
    % --- Grain, -108 degrees ---
    u_init = vt_grain * cos(th_108);
    v_init = vt_grain * sin(th_108);
    sol = rk4Solver(grain, rho_f,mu_f, uf0,h, x0,y0, u_init,v_init, dt, t_max, x_bin,y_bin, g);
    extremes.grain_m108 = sol.x(end);
    paths.grain_m108    = sol;
    
    % --- Chaff, -72 degrees ---
    u_init = vt_chaff * cos(th_72);
    v_init = vt_chaff * sin(th_72);
    sol = rk4Solver(chaff, rho_f,mu_f, uf0,h, x0,y0, u_init,v_init, dt, t_max, x_bin,y_bin, g);
    extremes.chaff_m72 = sol.x(end);
    paths.chaff_m72    = sol;
    
    % --- Chaff, -108 degrees ---
    u_init = vt_chaff * cos(th_108);
    v_init = vt_chaff * sin(th_108);
    sol = rk4Solver(chaff, rho_f,mu_f, uf0,h, x0,y0, u_init,v_init, dt, t_max, x_bin,y_bin, g);
    extremes.chaff_m108 = sol.x(end);
    paths.chaff_m108    = sol;
end

%============================================================================================= END ===========================================================================================
