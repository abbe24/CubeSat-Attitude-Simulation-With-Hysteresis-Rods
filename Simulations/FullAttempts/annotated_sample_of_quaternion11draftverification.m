%% anotated sample of quaternion11draftverification
% Truman Abbe

% Started 12/16/25
% Truman Abbe | Utah State University | truman.abbe23@gmail.com
% quaternion11draftverification.m

%%% NEAR PERFECT COPY OF GERHARDT SIMULATION


%%%%%%%%%%%%% RUN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial state y0 is given with scalar first quaternion:
% [w, x, y, z, w_x, w_y, w_z]
%y0 = [1 0 0 0 deg2rad(10) deg2rad(5) deg2rad(5)]';
y0 = [1 0 0 0 deg2rad(10) deg2rad(5) deg2rad(5)]';
%B_s = 0.027; H_c = 12; B_r = 0.004; 
%B_s = 0.86; H_c = 1.59; B_r = 0.35; %MEMEsat numbers
B_s = 0.027; H_c = 12; B_r = 0.004; %CCSWE

% define 30 days in seconds
t_end = 1 * 8 * 24 * 3600;%1*24*3600; 
t_span = 0:240:t_end;   % <- must come BEFORE you use it
%t_span = 0:10:5400; % 90 min in seconds

global B_hyst_prev B_hyst_history p 
B_hyst_prev = zeros(3,1);  % [X; Y; Z] rods
B_hyst_history = zeros(3, length(t_span));
p = (1 / H_c) * tan((pi * B_r) / (2 * B_s));
disp("this is p:")
disp(p)
disp(" ")










 

 




fprintf('Starting 30-day simulation... this may take a few minutes.\n');
tic; % starts timer

% tolerances are relaxed to speed up integration
%options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'MaxStep', 0.5);
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-7);%, 'MaxStep', 0.05);

% define 30 days in seconds
%t_end = 30 * 24 * 3600; 

% use a specified output time vector to save memory
% instead of saving every tiny step, save data every 120 seconds
%t_span = 0:120:t_end;

% stiff solver ode15s is used to accommodate magnetic torques  
[T, Y] = ode45(@myODE, t_span, y0, options);
%[T, Y] = ode15s(@myODE, t_span, y0, options);


elapsedTime = toc; % ends timer
fprintf('Simulation complete in %.2f seconds.\n', elapsedTime);
final_state = Y(end, :);
disp(T(end, :)); % time elapsed is displayed for debugging
disp(final_state); % final state is displayed for debugging




%%%%%% Store B_hyst history %%%%%%
for k = 1:length(T)
    % call a modified function to get B_hyst at this time without changing global memory
    B_hyst_history(:,k) = get_B_hyst_for_plot(Y(k,:)');
end



%%%%%% FINAL MISSION DASHBOARD: 14/30-DAY ANALYTICS %%%%%%%%%%%%%%%%%%%%%%%%%%
T_days = T(:) / 86400; % column vector
Y = Y(:,:); % ensures matrix is in correct shape



% downsampling factor for visual clarity (plots every 10th point)
% this prevents "solid block" artifacts in long-duration simulations
%idx = 1:20:size(Y,1);
%idx = 1:size(Y,1);
%idx = 1:60:size(Y,1);   % Plot every 10th point (1.0s) to keep the GUI fast
idx = 1:size(Y,1);     % Plot every point saved (which is 1 per minute)







% --- Hysteresis rod flux history ---
figure;
plot(T_days, B_hyst_history(1,:), 'r', 'LineWidth',1); hold on;
plot(T_days, B_hyst_history(2,:), 'g', 'LineWidth',1);
plot(T_days, B_hyst_history(3,:), 'b', 'LineWidth',1);
xlabel('Time (Days)'); ylabel('B_{hyst} (T)');
legend('Rod X','Rod Y','Rod Z'); title('Hysteresis Rod Flux vs Time');
grid on;











%disp(size(T_days))
%disp(size(Y))
%disp(idx(1:10))

% --- Minimal test to check plotting ---
figure; 
plot(T_days(idx), Y(idx,1), 'r', 'LineWidth', 1); 
title('Minimal Test: First Quaternion'); 
xlabel('Time (Days)'); ylabel('q_w'); 
grid on;

fig = figure('Name', 'Satellite Mission Dashboard: Long-Term Analysis', 'NumberTitle', 'off', 'Color', 'w');

% --- 1. Top Left: Raw Quaternions (Cleaned) ---
subplot(2,2,1);
plot(T_days(idx), Y(idx, 1:4), 'LineWidth', 1);
title('State Vector: Raw Quaternions');
xlabel('Time (Days)'); ylabel('Component Value');
legend('q_w', 'q_x', 'q_y', 'q_z', 'Location', 'northeast');
grid on; xlim([0 T_days(end)]);
ylim([-1.1 1.1]);

% --- 2. Top Right: High-Precision Norm Error ---
subplot(2,2,2);
q_norm = sqrt(sum(Y(:,1:4).^2, 2)); 
error_norm = 1 - q_norm;            
plot(T_days(idx), error_norm(idx), 'Color', [0, 0.447, 0.741]); 
title('High-Precision Norm Error (1 - ||q||)');
xlabel('Time (Days)'); ylabel('Error Magnitude');
grid on; xlim([0 T_days(end)]);

% --- 3. Bottom Left: Angular Velocity (The Physics) ---
subplot(2,2,3);
Y_deg = Y(:, 5:7) * (180/pi); 
w_total_deg = sqrt(sum(Y_deg.^2, 2));

% Plot individual components and the magnitude
plot(T_days, Y_deg, 'LineWidth', 1); hold on;
plot(T_days, w_total_deg, 'k', 'LineWidth', 2); 

% Zero reference line for clarity
yline(0, '-k', 'HandleVisibility', 'off'); 

% Detumble Threshold Line
yline(1, '--r', 'Detumble Threshold', 'LabelVerticalAlignment', 'bottom');

title('Physics: Angular Velocity Decay');
xlabel('Time (Days)'); ylabel('deg/sec');
legend('\omega_x', '\omega_y', '\omega_z', '|\omega|_{total}');
grid on; xlim([0 T_days(end)]);

% --- THE CHANGE IS HERE ---
% Manually set the y-axis to your requested range
ylim([-10, 15]);

% --- 4. Bottom Right: Unit Quaternion Health (Visual Fix) ---
subplot(2,2,4);
% We relax the ylim slightly to avoid the "Solid Red Block" effect
plot(T_days(idx), q_norm(idx), 'r', 'LineWidth', 1.5);
title('Unit Quaternion Health (Stability)');
xlabel('Time (Days)'); ylabel('Magnitude');
grid on; xlim([0 T_days(end)]);

% adjusting ylim to show a stable line rather than a thick band
% if error is 1e-4, this range shows the stability clearly.
ylim([0.999 1.001]); 

% Auto-save the figure as a high-res PNG for your report
drawnow;
saveas(fig, 'Final_30Day_Mission_Dashboard11ver1.png');



% --- NEW FIGURE: MAGNETIC ENVIRONMENT & TORQUE ANALYSIS ---
fig2 = figure('Name', 'Magnetic Environment & Control Torques', 'Color', 'w');

% 1. Earth Magnetic Field in ECI (Verification of getEarthField)
subplot(3,1,1);
B_eci_history = zeros(length(T), 3);
for k = 1:length(T)
    pos = circularPropagator(T(k)); 
    B_eci_history(k,:) = getEarthField(pos);
end
plot(T_days, B_eci_history * 1e6, 'LineWidth', 1); % Convert to microTesla
title('Earth Magnetic Field (ECI Frame)');
ylabel('\muT'); xlabel('Time (Days)');
legend('B_X', 'B_Y', 'B_Z'); grid on;

% 2. Earth Magnetic Field in Body Frame (Verification of R_bi)
% This shows what the satellite actually "sees" as it tumbles
subplot(3,1,2);
B_body_history = zeros(length(T), 3);
for k = 1:length(T)
    q_curr = Y(k, 1:4);
    R_bi = quat2rotm_manual(q_curr);
    B_body_history(k,:) = (R_bi' * B_eci_history(k,:)')';
end
plot(T_days, B_body_history * 1e6, 'LineWidth', 1);
title('Magnetic Field in Satellite Body Frame');
ylabel('\muT'); xlabel('Time (Days)');
legend('B_x', 'B_y', 'B_z'); grid on;

% 3. Component Torques (The "Brakes" vs the "Compass")
% This will show if your hysteresis rods are actually stronger than the tumble
subplot(3,1,3);
% Note: You'll need to extract tau_bar and tau_hyst from a loop 
% similar to the B_hyst_history loop you already have.
% ... (Code to recalculate and plot tau_bar and tau_hyst) ...
title('Control Torques: Hysteresis vs. Bar Magnet');
ylabel('N-m'); xlabel('Time (Days)');
grid on;





% Add this inside your post-processing loop
pointing_error = zeros(length(T), 1);
for k = 1:length(T)
    % 1. Get the field in body frame for this moment
    B_body_vec = B_body_history(k,:)'; 
    % 2. Our bar magnet is on the Body X-axis [1; 0; 0]
    body_x = [1; 0; 0];
    % 3. Calculate angle between them in degrees
    pointing_error(k) = rad2deg(acos(dot(body_x, B_body_vec) / norm(B_body_vec)));
end

% Plotting the Pointing Error
figure;
plot(T_days, pointing_error, 'm', 'LineWidth', 1);
title('Satellite Pointing Error (X-Axis vs B-Field)');
xlabel('Time (Days)'); ylabel('Degrees (°)');
grid on;






%%%%%%%%%%%% ODE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = myODE(t, y)

    %% Quaternion is unpacked.
    q_in = y(1:4);

    %% Hard normalization for quaternion and converted to row vector.
    % Quaternion drift occurs in variable integrators like ODE45 because
    % the solver truncates quaternion values. A stabalizing method is
    % needed to prevent drift. In this function, a two-step stabalizing
    % method is implemented. A hard normalization at the start ensures that
    % a perfect unit quaternion is used for the current timestep's physics
    % calculations.
    q_norm = norm(q_in);
    q = (q_in / q_norm)'; 
    
    %% Angular velocity is unpacked.
    w = y(5:7); 
    
    %% Physical constants of hysteresis rods
    mu_0 = 4 * pi * 10^(-7);
    rod_L = 0.095; % length
    rod_r = 0.0005; % radius
    V_hyst = pi * rod_r^2 * rod_L; % volume of one rod
    
    %% Inertia values for Gerhardt's CSSWE 3U CubeSat
    I_xx = 0.00551; 
    I_yy = 0.02552; 
    I_zz = 0.02565;
    I = diag([I_xx, I_yy, I_zz]);
    
    %% Magnetic flux density of Earth's field is generated and transformed
    %% to body frame.
    % A propagator function models a 600 km altitude orbit with inclination
    % of 51.6 deg. Given an elapsed time, it returns the position along the
    % orbit. An Earth field function is used that is a simplified dipole
    % model meaning it assumes spherical symmetry and zero tilt. A more
    % realistic field function will be implemented soon. A frame
    % transformation function is used that converts from body to inertial
    % frame, requiring transpose (') for inertial to body.
    pos = circularPropagatorGerhardt(t);
    B_eci = getEarthField(pos);
    R_bi = quat2rotm_manual(q); 
    B_body = R_bi' * B_eci;
    
    %% Hysteresis torque calculation.
    rod_axes = eye(3); 
    tau_hyst = zeros(3,1);
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
    
    %% Kinematics: q_dot (Quaternion Rate).
    % This is the quaternion kinematic equation using scalar first
    % convention and hamilton product.
    qdot_raw = 0.5 * quatmultiply(q, [0, w']);
    
    %% Baumgarte Stabilization: Gently nudges the norm back to 1.0.
    % This is the second norm stabilizing technique.
    k_stabilize = 0.1; 
    norm_error = 1 - q_norm; % Used the original norm to calculate error.
    qdot = qdot_raw + k_stabilize * norm_error * q;
    
    %% Permanent Magnet: 0.3 Am^2 aligned with Body X-axis.
    m_bar = [0.3; 0; 0]; 
    
    %% Permanent Magnet Torque: tau = m x B.
    tau_bar = cross(m_bar, B_body);
    
    %% Dynamics: w_dot (Angular Acceleration).
    % External torques tau_hyst and tau_bar are included in calcuation.
    w_col = y(5:7);
    wdot = I \ (tau_hyst + tau_bar - cross(w_col, I * w_col));
    
    %% Repack into state matrix.
    dydt = [qdot(:); wdot(:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%HAMILTON PRODUCT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = quatmultiply(p, q)
    % r = quatmultiply(p, q)
    %   Performs the Hamilton product (p * q) for quaternions.
    %   *** ASSUMES SCALAR-FIRST CONVENTION: q = [qw, qx, qy, qz] ***
    %   
    %   Inputs p and q MUST be 1x4 row vectors.
    %   Output r is a 1x4 row vector.

    % --- Unpack components of the first quaternion (p) ---
    % In your case, p = Omega = [0, wx, wy, wz]
    ps = p(1);   % Scalar part (p_s)
    pv = p(2:4); % Vector part (p_v)

    % --- Unpack components of the second quaternion (q) ---
    % In your case, q = Orientation = [qw, qx, qy, qz]
    qs = q(1);   % Scalar part (q_s)
    qv = q(2:4); % Vector part (q_v)

    % --- Calculate the Product Components ---
    
    % 1. Calculate the NEW Scalar Part (r_s):
    % r_s = p_s*q_s - dot(p_v, q_v)
    r_s = ps * qs - dot(pv, qv);
    
    % 2. Calculate the NEW Vector Part (r_v):
    % r_v = p_s*q_v + q_s*p_v + cross(p_v, q_v)
    r_v = ps * qv + qs * pv + cross(pv, qv);
    
    % Combine product: r = [r_s, r_v] (SCALAR-FIRST order)
    r = [r_s, r_v];
end


%%%%%%%%%%%%%%%%%%%%%QUAT2ROTM FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R_bi = quat2rotm_manual(q)
% QUAT2ROTM_MANUAL Converts a quaternion to a 3x3 rotation matrix.
% This replaces the need for the Navigation or Robotics toolboxes.
%
% Input:  q - A 4x1 or 1x4 vector [w, x, y, z] (Scalar first)
% Output: R_bi - A 3x3 Rotation Matrix

    % 1. Ensure q is a unit quaternion (essential for 30-day stability)
    %q = q / norm(q); this is redundant and taken care of my Baumgarte

    % 2. Extract components
    w = q(1); 
    x = q(2); 
    y = q(3); 
    z = q(4);

    % 3. Calculate the rotation matrix (Hamiltonian convention)
    R_bi = [1 - 2*y^2 - 2*z^2,   2*x*y - 2*w*z,       2*x*z + 2*w*y;
            2*x*y + 2*w*z,       1 - 2*x^2 - 2*z^2,   2*y*z - 2*w*x;
            2*x*z - 2*w*y,       2*y*z + 2*w*x,       1 - 2*x^2 - 2*y^2];
end


%%%%%%%%%%%%%%%%%MAGNETIC FIELD FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B_eci = getEarthField(r_eci)
    % --- Constants ---
    M_earth = 7.94e22;      % Earth's magnetic moment (Amp-m^2)
    mu0_4pi = 1e-7;         % Constant (mu0 / 4*pi)
    
    % --- Distance Calculation ---
    r_mag = norm(r_eci);    % Distance from Earth center (meters)
    r_hat = r_eci / r_mag;  % Unit vector pointing to the satellite
    
    % --- Magnetic North Pole Orientation ---
    % For a standard model, the dipole is aligned near the Z-axis
    % You can refine this later with an 11.5 degree tilt if needed.
    m_hat = [0; 0; 1]; 
    
    % --- The Dipole Equation ---
    % B = (mu0/4pi) * (3*r_hat*(m_dot_r) - m) / r^3
    m_dot_r = dot(m_hat, r_hat);
    B_eci = (mu0_4pi * M_earth / r_mag^3) * (3 * r_hat * m_dot_r - m_hat);
end


%%%%%%%%%%%%%%%%%%%ORBITAL PROPAGATOR FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function position = circularPropagator(t)
    R = 6371 * 10^3 + 400 * 10^3;
    i = deg2rad(51.6);
    mu = 3.985892 * 10^14;
    n = sqrt(mu / R^3);
    theta0 = 0;
    theta = n * t + theta0;
    x = R * cos(theta);
    y = R * sin(theta) * cos(i);
    z = R * sin(theta) * sin(i);

    position = [x; y; z];
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


%%%%%%%%%%%%%%%%%%CALCULATE_B_HYST FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function m_hyst = calculate_m_hyst(y, B_body, V_hyst, B_s, H_c, B_r)
% y has to be a vertical matrix
    
    % mu_0 = 4 * pi * 10^(-7);
    % H = B_body / mu_0;
    % p = (1 / H_c) * tan((pi * B_r) / (2 * B_s));

    % H_dot = -cross(y(5:7), H); % assuming y(5:7) is your [wx; wy; wz]
    % find the sign (s) for each axis (x, y, z)
    % s will be a 3x1 vector of 1, -1, or 0
    
    % first sign flip method
    % s = sign(H_dot);
    % fallback: if sign is 0 (not moving), assume increasing to avoid errors
    % s(s == 0) = 1;

    % Hdot_scale = 0.05;      % A/m/s  (tuned to LEO dynamics)
    % s = tanh(H_dot / Hdot_scale);

    % continuous approximation - better sign flip method
    % eps = 1e-6; % smoothing parameter
    % s = H_dot ./ sqrt(H_dot.^2 + eps^2);

    % Hdot_ref = 5e-4;  % A/m/s, consistent with Gerhardt test rates
    % s = tanh(H_dot / Hdot_ref);

    % B_hyst = (2 / pi) * B_s * atan(p * (H - s * H_c));

    % m_hyst = (B_hyst * V_hyst) / mu_0;

% end

function B_hyst = calculate_B_hyst(y, B_body, rod_axis)
    global p
    mu_0 = 4 * pi * 10^(-7);
    %B_s = 0.86; H_c = 1.59; % MEMEsat saturation and coercivity
    B_s = 0.027; H_c = 12; B_r = 0.004; %CCSWE


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

function B_hyst = get_B_hyst_for_plot(y)
    % We must recalculate B_body for each state point to plot correctly
    pos = circularPropagator(0); % Time is less critical for B_body magnitude
    B_eci = getEarthField(pos);
    R_bi = quat2rotm_manual(y(1:4)'); 
    B_body = R_bi' * B_eci; 

    rod_axes = eye(3);    
    B_hyst = zeros(3,1);
    for i = 1:3
        % Pass exactly 3 arguments to match the function definition
        B_hyst(i) = calculate_B_hyst(y, B_body, rod_axes(:,i));
    end
end 