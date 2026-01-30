%% implicit lgvi 2 with magnet
clear; clc;
% 1. Simulation Parameters
h = 0.1;               
t_end = 8*24*3600;             
steps = t_end / h;
t_span = 0:h:t_end-h;
% 2. Physical Inputs
% Initial velocity: 5 deg/s on all axes
w = [deg2rad(10); deg2rad(5); deg2rad(5)]; 
%J = diag([0.03, 0.03, 0.03]);   
I_xx = 0.00551; 
I_yy = 0.02552; 
I_zz = 0.02565;
% Define the matrix
J = diag([I_xx, I_yy, I_zz]);
invJ = inv(J);
R = eye(3);
M = J * w;         
% --- ADDED: Orbit & Magnet Parameters ---
R_earth = 6378e3;
r_orbit = R_earth + 500e3;      % 500 km Altitude
mean_motion = sqrt(3.986e14 / r_orbit^3); 
mu_body = [0.3; 0; 0];          % Magnetic dipole (Am^2) aligned with Body X
%% Physical constants of hysteresis rods
mu_0 = 4 * pi * 10^(-7);
rod_L = 0.095; % length
rod_r = 0.0005; % radius
V_hyst = pi * rod_r^2 * rod_L; % volume of one rod
% 3. Pre-allocate Data
omega_deg_hist = zeros(steps, 3); 
quat_hist      = zeros(steps, 4); % [qw, qx, qy, qz]
E_err_hist     = zeros(steps, 1);
M_mag_hist     = zeros(steps, 1);
ortho_err_hist = zeros(steps, 1);
% --- UPDATED: Energy Check (Kinetic + Potential) ---
% Need initial B-field for initial potential energy
r_eci_0 = circularPropagatorGerhardt(0); 
B_inertial_0 = getImprovedEarthField(r_eci_0, 0);
U_pot_0 = -dot(mu_body, R' * B_inertial_0);
E0 = (0.5 * M' * (invJ * M)) + U_pot_0;
M0_mag = norm(M);

% --- ADDED: Wall Clock Init ---
fprintf('Starting 8-day LGVI simulation...\n');
tic; 
last_print_time = 0;

% 4. Main Integration Loop
for k = 1:steps
    t = (k-1) * h;
    
    % --- ADDED: Orbital Position & Magnetic Field ---
    r_eci = circularPropagatorGerhardt(t);
    B_inertial = getImprovedEarthField(r_eci, t);
    
    % --- ADDED: Current Torque ---
    B_body = R' * B_inertial;
    
    % --- INSERTED: Hysteresis torque calculation (For Current Step) ---
    rod_axes = eye(3); 
    tau_hyst = zeros(3,1);
    
    % Create dummy state 'y' for your function [q_dummy; w]
    % Note: w must be current angular velocity
    w_curr = invJ * M;
    y = [0;0;0;1; w_curr]; 
    
    for i = 1:3
        rod_dir = rod_axes(:,i);
        % Determine how many rods are on this axis (CSSWE Configuration)
        if i == 1
            num_rods = 0; % CSSWE used a Bar Magnet on X, no hysteresis rods.
        else
            num_rods = 2; % Two rods each for Y and Z axes.
        end
    
        if num_rods > 0
            % This calculates the magnetic flux density of hystersis rods
            % using the Henretty hysteresis model.
            B_rod = calculate_B_hyst(y, B_body, rod_dir);
        
            % The magnetic moment is calculated using flux density and
            % rod volume.
            m_mag = (B_rod * (V_hyst * num_rods)) / mu_0;
            m_vec = m_mag * rod_dir;
        
            % Torque is calculated using T = m x B_external
            tau_hyst = tau_hyst + cross(m_vec, B_body);
        end
    end
    
    % Permanent Magnet Torque: tau = m x B.
    tau_bar = cross(mu_body, B_body);
    
    % Total Torque (Current)
    tau_total = tau_hyst + tau_bar;
    
    % --- UPDATED: Effective Momentum (Predictor) ---
    M_eff = M + (0.5 * h * tau_total);
    
    % --- STEP 1: Exact Newton-Raphson Implicit Solve ---
    f = h * (invJ * M_eff); % Initial Guess using M_eff
    tol = 1e-12;        % Tolerance for machine-precision consistency
    
    for i = 1:7 
        theta = norm(f);
        theta2 = theta^2;
        
        % 1. Compute Coefficients (c1, c2) and Derivatives (dc1/theta, dc2/theta)
        if theta < 1e-4
            % Taylor Series (to prevent division by zero)
            c1 = 1 - theta2/6 + theta2^2/120;
            c2 = 0.5 - theta2/24 + theta2^2/720;
            
            % Analytic derivatives divided by theta: (dc/dtheta) / theta
            dc1_th = -1/3 + theta2/30;  
            dc2_th = -1/12 + theta2/180;
        else
            s = sin(theta); 
            c = cos(theta);
            c1 = s / theta;
            c2 = (1 - c) / theta2;
            
            % Exact analytic derivatives
            dc1_val = (theta*c - s) / theta2;
            dc1_th  = dc1_val / theta;
            
            dc2_val = (theta*s - 2*(1 - c)) / (theta2*theta);
            dc2_th  = dc2_val / theta;
        end
        
        % 2. Helper Variables
        Jf = J * f;
        f_cross_Jf = cross(f, Jf);
        
        % Skew-symmetric matrices
        S_f  = [0, -f(3), f(2); f(3), 0, -f(1); -f(2), f(1), 0];
        S_Jf = [0, -Jf(3), Jf(2); Jf(3), 0, -Jf(1); -Jf(2), Jf(1), 0];
        
        % 3. Calculate Residual G (UPDATED with M_eff)
        G = (c1 * Jf) + (c2 * f_cross_Jf) - (h * M_eff);
        
        % Dynamic exit if converged
        if norm(G) < tol
            break;
        end
        
        % 4. Exact Jacobian Construction
        % Term A: Derivatives related to J*f
        % d(c1*Jf)/df = c1*J + (dc1/theta) * (Jf * f')
        Jac_A = c1 * J + dc1_th * (Jf * f');
        
        % Term B: Derivatives related to f x Jf
        % d(c2 * (f x Jf))/df = c2 * d(f x Jf)/df + d(c2)/df * (f x Jf)
        % Using identity: d(a x b)/da = skew(a)*db/da - skew(b)
        % Here: d(f x Jf)/df = S_f * J - S_Jf
        Jac_B = c2 * (S_f * J - S_Jf) + dc2_th * (f_cross_Jf * f');
        
        Jac = Jac_A + Jac_B;
        
        % 5. Newton Update
        f = f - Jac \ G; 
    end
    
    % --- STEP 2: Compute F (Rodrigues) ---
    theta = norm(f);
    S_f = [0, -f(3), f(2); f(3), 0, -f(1); -f(2), f(1), 0];
    if theta < 1e-6
        c1 = 1 - theta^2/6; c2 = 0.5 - theta^2/24;
    else
        c1 = sin(theta)/theta; c2 = (1-cos(theta))/(theta^2);
    end
    F = eye(3) + c1 * S_f + c2 * (S_f * S_f);
    
    % --- STEP 3 & 4: Updates ---
    R = R * F; 
    
    % --- ADDED: Momentum Update (Trapezoidal Torque) ---
    B_body_next = R' * B_inertial;
    
    % --- INSERTED: Hysteresis torque calculation (For Next Step) ---
    tau_hyst_next = zeros(3,1);
    
    % We use 'w_curr' again here to avoid a second implicit loop (semi-implicit approx)
    % or we could use the predicted velocity: w_next = invJ * (F' * M);
    y = [0;0;0;1; w_curr]; 
    
    for i = 1:3
        rod_dir = rod_axes(:,i);
        if i == 1, num_rods = 0; else, num_rods = 2; end
    
        if num_rods > 0
            B_rod = calculate_B_hyst(y, B_body_next, rod_dir);
            m_mag = (B_rod * (V_hyst * num_rods)) / mu_0;
            m_vec = m_mag * rod_dir;
            tau_hyst_next = tau_hyst_next + cross(m_vec, B_body_next);
        end
    end
    
    tau_bar_next = cross(mu_body, B_body_next);
    tau_total_next = tau_hyst_next + tau_bar_next;
    
    M = F' * M_eff + (0.5 * h * tau_total_next);
    
    % --- 5. DATA LOGGING ---
    omega_rad = invJ * M;
    omega_deg_hist(k, :) = rad2deg(omega_rad)';
    
    % Robust R to Quaternion Conversion
    tr = trace(R);
    if (tr > 0)
        S = sqrt(tr + 1.0) * 2;
        q_raw = [0.25 * S, (R(3,2)-R(2,3))/S, (R(1,3)-R(3,1))/S, (R(2,1)-R(1,2))/S];
    else
        [~, idx] = max(diag(R));
        if idx == 1
            S = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3)) * 2;
            q_raw = [(R(3,2)-R(2,3))/S, 0.25*S, (R(1,2)+R(2,1))/S, (R(1,3)+R(3,1))/S];
        elseif idx == 2
            S = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3)) * 2;
            q_raw = [(R(1,3)-R(3,1))/S, (R(1,2)+R(2,1))/S, 0.25*S, (R(2,3)+R(3,2))/S];
        else
            S = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2)) * 2;
            q_raw = [(R(2,1)-R(1,2))/S, (R(1,3)+R(3,1))/S, (R(2,3)+R(3,2))/S, 0.25*S];
        end
    end
    
    % --- CONTINUITY CORRECTION ---
    % If the dot product with previous q is negative, flip signs to prevent jumps
    if k > 1
        if dot(q_raw, quat_hist(k-1, :)) < 0
            q_raw = -q_raw;
        end
    end
    quat_hist(k, :) = q_raw;
    
    % UPDATED: Energy Check (Total Energy = Kinetic + Potential)
    U_pot = -dot(mu_body, B_body_next);
    E_total = (0.5 * M' * omega_rad) + U_pot;
    
    E_err_hist(k)     = E_total - E0;
    M_mag_hist(k)     = norm(M) - M0_mag;
    ortho_err_hist(k) = norm(R' * R - eye(3), 'fro');

    % --- ADDED: Wall Clock Update Every 10 Seconds ---
    current_time = toc;
    if current_time - last_print_time > 10
        fprintf('Progress: %.1f%%  |  Sim Time: %.2f Days  |  Wall Clock: %.0fs\n', ...
            (k/steps)*100, t/(3600*24), current_time);
        last_print_time = current_time;
    end
end
fprintf('Simulation Complete. Total Wall Time: %.2f seconds.\n', toc);

% 6. Create Tabbed Figure
hFig = figure('Name', 'LGVI Comprehensive Satellite Analysis', 'Color', 'w', ...
              'Units', 'Normalized', 'OuterPosition', [0.1 0.1 0.8 0.8]);
hTabGroup = uitabgroup(hFig);
% --- TAB 1: Motion Data ---
tab1 = uitab(hTabGroup, 'Title', 'Motion Data');
% Subplot 1: Angular Velocity
ax1 = subplot(2,1,1, 'Parent', tab1);
omega_total_mag = sqrt(sum(omega_deg_hist.^2, 2));
plot(ax1, t_span, omega_deg_hist, 'LineWidth', 1.5); 
hold(ax1, 'on');
plot(ax1, t_span, omega_total_mag, 'k', 'LineWidth', 2, 'DisplayName', '|\omega|_{total}');
title(ax1, 'Angular Velocity (Deg/s)'); 
ylabel(ax1, '^\circ/s'); 
grid(ax1, 'on');
legend(ax1, '\omega_x', '\omega_y', '\omega_z', '|\omega|_{total}', 'Location', 'northeastoutside');
% Subplot 2: Smooth Quaternions
ax2 = subplot(2,1,2, 'Parent', tab1);
plot(ax2, t_span, quat_hist, 'LineWidth', 1.5);
title(ax2, 'Smooth Continuous Quaternions (Corrected for Antipodal Flips)'); 
ylabel(ax2, 'Value'); 
xlabel(ax2, 'Time (s)'); 
grid(ax2, 'on');
legend(ax2, 'q_w', 'q_x', 'q_y', 'q_z', 'Location', 'northeastoutside');
% --- TAB 2: Physics Integrity ---
tab2 = uitab(hTabGroup, 'Title', 'Physics Integrity');
ax3 = subplot(2,1,1, 'Parent', tab2);
plot(ax3, t_span, E_err_hist, 'r', 'LineWidth', 1, 'DisplayName', 'Energy Error'); 
hold(ax3, 'on');
plot(ax3, t_span, M_mag_hist, 'b--', 'LineWidth', 1, 'DisplayName', 'Momentum Error');
title(ax3, 'Conservation Deviations (Machine Precision Zoom)'); 
ylabel(ax3, 'Absolute Error'); 
grid(ax3, 'on'); 
% Automatic zoom to show machine epsilon (10^-15 level)
axis(ax3, 'tight'); 
if max(abs(E_err_hist)) < 1e-12
    ylim(ax3, [-2e-15, 2e-15]); 
end
legend(ax3, 'show', 'Location', 'northeastoutside');
ax4 = subplot(2,1,2, 'Parent', tab2);
semilogy(ax4, t_span, ortho_err_hist + 1e-18, 'k', 'LineWidth', 1);
title(ax4, 'Orthogonality Error (Lie Group Preservation)'); 
ylabel(ax4, 'Log Error ||R^T R - I||'); 
xlabel(ax4, 'Time (s)'); 
grid(ax4, 'on');
function B_eci = getImprovedEarthField(r_eci, time_seconds)
    % --- Constants ---
    M_earth = 7.94e22;      
    mu0_4pi = 1e-7;         
    
    % --- Magnetic North Pole Orientation (Tilted Dipole) ---
    % The magnetic pole is tilted ~11.5 degrees from the rotation axis.
    % It rotates with the Earth, so m_hat changes in the ECI frame over 24 hours.
    tilt = deg2rad(11.5); 
    omega_earth = 7.292115e-5; % rad/s
    
    % Rotate the dipole vector around the Earth's Z-axis based on time
    % This simulates the North Pole "wobbling" from the ECI perspective
    m_hat = [sin(tilt)*cos(omega_earth * time_seconds);
             sin(tilt)*sin(omega_earth * time_seconds);
             cos(tilt)];
    
    % --- Distance Calculation ---
    r_mag = norm(r_eci);    
    r_hat = r_eci / r_mag;  
    
    % --- The Dipole Equation ---
    m_dot_r = dot(m_hat, r_hat);
    B_eci = (mu0_4pi * M_earth / r_mag^3) * (3 * r_hat * m_dot_r - m_hat);
end
function position = circularPropagatorGerhardt(t)
    % --- Gerhardt / CSSWE Orbital Parameters ---
    % Earth Radius: 6371 km
    % Altitude: 600 km (Gerhardt uses 600 km for all CSSWE simulations)
    R = (6371 + 600) * 10^3; 
    
    % Inclination: 55 degrees (Gerhardt explicitly assumes 55 deg)
    i = deg2rad(55.0);
    
    % Standard Gravitational Parameter (Earth)
    mu = 3.986004418 * 10^14; 
    
    % 1. Calculate Mean Motion (n)
    % This is the angular velocity of the orbit
    n = sqrt(mu / R^3);
    
    % 2. Calculate Argument of Latitude (u)
    % Gerhardt defines 'u' as the position in the orbit relative to the 
    % ascending node. u = n*t + u0.
    u0 = 0; % Starting at the equator (ascending node)
    u = n * t + u0;
    
    % 3. Calculate Position in ECI Coordinates
    % This rotates the circular orbit by the 55 degree inclination
    x = R * cos(u);
    y = R * sin(u) * cos(i);
    z = R * sin(u) * sin(i);
    
    position = [x; y; z];
end
function B_hyst = calculate_B_hyst(y, B_body, rod_axis)
    %global p
    mu_0 = 4 * pi * 10^(-7);
    %B_s = 0.86; H_c = 1.59; % MEMEsat saturation and coercivity
    B_s = 0.027; H_c = 12; B_r = 0.004; %CCSWE
    p = (1 / H_c) * tan((pi * B_r) / (2 * B_s));
    % 1. Project the magnetic field onto the rod axis
    H_vec = B_body / mu_0;
    H_parallel = dot(H_vec, rod_axis);
    
    % 2. Calculate H_dot (stateless)
    % This is the change in the field as seen by the rotating rod
    w = y(5:7);
    H_dot_vec = -cross(w, H_vec); 
    H_dot_parallel = dot(H_dot_vec, rod_axis);
    
% 3. Revised Switching Logic
    % s must be +1 when H_dot > 0 (increasing) to subtract Hc
    % s must be -1 when H_dot < 0 (decreasing) to add Hc
    Hdot_scale = 1e-6; 
    %Hdot_scale = 1e-12; % added line
    s = tanh(H_dot_parallel / Hdot_scale); % REMOVED the negative sign here
    
    % 4. The Henretty sigmoid equation
    % To create the lag, we use (H - s*Hc)
    B_hyst = (2 / pi) * B_s * atan(p * (H_parallel - s * H_c));
end