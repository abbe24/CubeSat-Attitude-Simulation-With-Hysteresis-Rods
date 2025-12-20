% Started 12/16/25
% Truman Abbe | Utah State University | truman.abbe23@gmail.com


y0 = [1 0 0 0 1 2 3]';
fprintf('Starting 30-day simulation... this may take a few minutes.\n');
tic; % Start timer
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-11);
[T, Y] = ode45(@myODE, [0 3600], y0, options);
elapsedTime = toc;
fprintf('Simulation complete in %.2f seconds.\n', elapsedTime);
final_state = Y(end, :);
disp(T(end, :));
disp(final_state);

%%%%%%DASHBOARD DISPLAY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- 1. Create a Unified Dashboard Figure ---
figure('Name', 'Satellite Mission Dashboard: 30-Day Simulation', 'NumberTitle', 'off', 'Color', 'w');

% --- 2. Top Left: Raw Quaternions (Subplot 1) ---
subplot(2,2,1);
plot(T, Y(:,1:4), 'LineWidth', 1.5);
title('State Vector: Raw Quaternions');
xlabel('Time (s)'); ylabel('Component Value');
legend('q_w', 'q_x', 'q_y', 'q_z');
grid on;

% --- 3. Top Right: High-Precision Norm Error (Subplot 2) ---
subplot(2,2,2);
q_norm = sqrt(sum(Y(:,1:4).^2, 2)); % Calculate ||q||
error_norm = 1 - q_norm;            % Calculate deviation
plot(T, error_norm, 'Color', [0, 0.447, 0.741]); 
title('High-Precision Norm Error (1 - ||q||)');
xlabel('Time (s)'); ylabel('Error Magnitude');
grid on;

% --- 4. Bottom Left: Angular Velocity in DEGREES (Subplot 3) ---
subplot(2,2,3);
% Conversion: Radians to Degrees (1 rad = ~57.3 degrees)
Y_deg = Y(:, 5:7) * (180/pi); 
% Calculate the total magnitude in degrees
w_total_deg = sqrt(sum(Y_deg.^2, 2));

% Plot individual components (x, y, z)
plot(T, Y_deg, 'LineWidth', 1.2); 
hold on;
% Plot TOTAL magnitude as SOLID BLACK line
plot(T, w_total_deg, 'k', 'LineWidth', 2); 

title('Physics: Angular Velocity (Degrees/sec)');
xlabel('Time (s)'); ylabel('deg/sec');
legend('\omega_x', '\omega_y', '\omega_z', '|\omega|_{total}');
grid on;

% --- 5. Bottom Right: Unit Quaternion Health (Subplot 4) ---
subplot(2,2,4);
plot(T, q_norm, 'r', 'LineWidth', 1.5);
title('Unit Quaternion Health (||q||)');
xlabel('Time (s)'); ylabel('Magnitude');
% Use tight limits to show the stability near 1.0
ylim([0.999999 1.000001]); 
grid on;






%%%%%%%%%%%%ODE FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = myODE(t, y)
    % y is [qw, qx, qy, qz, wx, wy, wz] column vector
    
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
    C = 0.0001;          % Damping constant (smaller = slower decay)
    tau = -C * w;        % Torque opposes the current angular velocity

    %wdot = (tau - cross(w, I*w')) / I;
    wdot = (I \ (tau' - cross(w', I*w')))';

    % Baumgarte Stabilization
    % This forces the norm back to 1.0 if it drifts
    k = 20; % Gain factor
    norm_error = 1 - (q(1)^2 + q(2)^2 + q(3)^2 + q(4)^2);
    q_dot_corrected = qdot + k * norm_error * q;
    
    % pack derivatives into column for ode45
    dydt = [q_dot_corrected'; wdot'];
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
