% Started 12/16/25
% Truman Abbe | Utah State University | truman.abbe23@gmail.com
% quaternion2.m


y0 = [1 0 0 0 1 2 3]';
fprintf('Starting 30-day simulation... this may take a few minutes.\n');
tic; % Start timer
%options = odeset('RelTol', 1e-8, 'AbsTol', 1e-11);
% Simulation starts from 0 deg and 400 km altitude set in circularPropagator
%[T, Y] = ode45(@myODE, [0 3600], y0, options);


% Relaxing tolerances slightly to speed up long-term integration
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 5);

% Define 7 days in seconds
t_end = 30 * 24 * 3600; 

% Use a specified output time vector to save memory
% Instead of saving every tiny step, save data every 60 seconds
t_span = 0:120:t_end;

[T, Y] = ode45(@myODE, t_span, y0, options);


elapsedTime = toc;
fprintf('Simulation complete in %.2f seconds.\n', elapsedTime);
final_state = Y(end, :);
disp(T(end, :));
disp(final_state);





%%%%%% FINAL MISSION DASHBOARD: 14/30-DAY ANALYTICS %%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion: Time to Days
T_days = T / 86400;

% Downsampling factor for visual clarity (plots every 10th point)
% This prevents "solid block" artifacts in long-duration simulations.
idx = 1:20:length(T);

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

plot(T_days, Y_deg, 'LineWidth', 1); hold on;
plot(T_days, w_total_deg, 'k', 'LineWidth', 2); 

% Detumble Threshold Line
yline(1, '--r', 'Detumble Threshold', 'LabelVerticalAlignment', 'bottom');

title('Physics: Angular Velocity Decay');
xlabel('Time (Days)'); ylabel('deg/sec');
legend('\omega_x', '\omega_y', '\omega_z', '|\omega|_{total}');
grid on; xlim([0 T_days(end)]);
ylim([0 max(w_total_deg)*1.1]); % Ensures magnitude is always visible

% --- 4. Bottom Right: Unit Quaternion Health (Visual Fix) ---
subplot(2,2,4);
% We relax the ylim slightly to avoid the "Solid Red Block" effect
plot(T_days(idx), q_norm(idx), 'r', 'LineWidth', 1.5);
title('Unit Quaternion Health (Stability)');
xlabel('Time (Days)'); ylabel('Magnitude');
grid on; xlim([0 T_days(end)]);

% Adjusting ylim to show a stable line rather than a thick band
% If your error is 1e-4, this range shows the stability clearly.
ylim([0.999 1.001]); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auto-save the figure as a high-res PNG for your report
saveas(fig, 'Final_30Day_Mission_Dashboard.png');






%%%%%%%%%%%%ODE FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = myODE(t, y)
    % y is [qw, qx, qy, qz, wx, wy, wz] column vector

    pos = circularPropagator(t);
    B_eci = getEarthField(pos);
    % y(1:4)' is your 1x4 quaternion vector [w, x, y, z]
    R_bi = quat2rotm_manual(y(1:4)');
    B_body = R_bi * B_eci;
    
    V_hyst = pi * (0.0005)^2 * (0.095) * 4;
    B_s = 0.027; H_c = 12; B_r = 0.004; % apparent values
    m_hyst = calculate_m_hyst(y, B_body, V_hyst, B_s, H_c, B_r);

    tau = cross(m_hyst, B_body);
    
    % unpack quaternion and angular velocity into row vectors
    q = y(1:4)';
    w = y(5:7)';
    
    % get pure quaternion row
    wq = [0 w];
    
    % compute qdot
    qdot = 0.5 * quatmultiply(q, wq);

    % compute wdot
    Ix = 0.01; Iy = 0.01; Iz = 0.01;
    I = diag([Ix Iy Iz]);
    
    %tau = [0 0 0];  %example torque
    % --- IMPLEMENT DAMPING HERE ---
    %C = 0.0001;          % Damping constant (smaller = slower decay)
    %tau = -C * w;        % Torque opposes the current angular velocity

    %wdot = (tau - cross(w, I*w')) / I;
    wdot = (I \ (tau' - cross(w', I*w')))';

    % Baumgarte Stabilization
    % This forces the norm back to 1.0 if it drifts
    k = 20; % Gain factor
    norm_error = 1 - (q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2);
    q_dot_corrected = qdot + k * norm_error * q;
    
    % pack derivatives into column for ode45
    %dydt = [q_dot_corrected'; wdot'];
    %dydt = [q_dot_corrected(:); wdot(:)];
    %dydt = [q_dot_corrected(:); wdot(:)];
    % Ensure we only take the exact elements needed (4 for q, 3 for w)
    % This prevents the "13 elements" error
    final_q_dot = q_dot_corrected(1:4); 
    final_w_dot = wdot(1:3);

    % Pack into a single 7x1 COLUMN vector
    dydt = [final_q_dot(:); final_w_dot(:)];
end

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

%%%%%%%%%%%%%%%%%NORM FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y_normalized = normalize_state(y)
    % y_normalized = normalize_state(y)
    %   Normalizes the quaternion portion of the 7-element state vector y.
    %   *** Input y MUST be a COLUMN vector: y = [qw; qx; qy; qz; wx; wy; wz] ***

    % 1. Extract the quaternion components 
    q = y(1:4); % q should be 4x1 (column vector)
    %disp(q);

    % 2. Calculate the magnitude and normalize
    q_norm = norm(q); 
    %disp(q_norm);
    %if q_norm == 0
        % ... warning handling ...
    %    q_normalized = q;
    %else
    q_normalized = q / q_norm; % q_normalized is 4x1 (column vector)
 

    % 3. Extract the angular velocity
    w = y(5:7); 

    % 4. Reassemble the final state vector using vertical concatenation (;)
    % We ensure w is a 3x1 column vector by transposing it, 
    % just in case it was a 1x3 row vector after slicing.
    
    % *** CORRECTED LINE ***
    y_normalized = [q_normalized; w(:)]; 
    % OR: y_normalized = [q_normalized; w']; if w was row vector
    % OR: y_normalized = [q_normalized; w]; if w was column vector.
    
    % The w(:) syntax forces w into a column vector regardless of its current shape.
    % This is the safest way to ensure compatibility for vertical concatenation.
end


%%%%%%%%%%%%%%%%%%%%%QUAT2ROTM FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R_bi = quat2rotm_manual(q)
% QUAT2ROTM_MANUAL Converts a quaternion to a 3x3 rotation matrix.
% This replaces the need for the Navigation or Robotics toolboxes.
%
% Input:  q - A 4x1 or 1x4 vector [w, x, y, z] (Scalar first)
% Output: R_bi - A 3x3 Rotation Matrix

    % 1. Ensure q is a unit quaternion (essential for 30-day stability)
    q = q / norm(q);

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


%%%%%%%%%%%%%%%%%%CALCULATE_B_HYST FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m_hyst = calculate_m_hyst(y, B_body, V_hyst, B_s, H_c, B_r)
% y has to be a vertical matrix
    
    mu_0 = 4 * pi * 10^(-7);
    H = B_body / mu_0;
    p = (1 / H_c) * tan((pi * B_r) / (2 * B_s));

    H_dot = -cross(y(5:7), H); % assuming y(5:7) is your [wx; wy; wz]
    % find the sign (s) for each axis (x, y, z)
    % s will be a 3x1 vector of 1, -1, or 0
    s = sign(H_dot);

    % fallback: if sign is 0 (not moving), assume increasing to avoid errors
    s(s == 0) = 1;

    B_hyst = (2 / pi) * B_s * atan(p * (H - s * H_c));

    m_hyst = (B_hyst * V_hyst) / mu_0;

end