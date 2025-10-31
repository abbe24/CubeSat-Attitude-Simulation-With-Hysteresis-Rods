%% MAIN SCRIPT: Rigid Body Dynamics Solver
% This section defines initial conditions, parameters, and solves the ODE.

% --- 1. DEFINE INITIAL CONDITIONS ---

% Initial Angles: 0 deg, 0 deg, 0 deg
theta_initial_deg = [10; 15; 20]; 
deg_to_rad = pi/180;
theta_initial_rad = theta_initial_deg * deg_to_rad; 

% Initial Rotational Velocities: Must be non-zero to generate complex motion
omega_initial = [1.0; 1.1; 1.2]; % Example: Initial spin about x and z axes (rad/s)

% Assemble the complete initial state vector Y0 (required by ode45)
% Y0 = [theta1; theta2; theta3; omega_x; omega_y; omega_z]
Y0 = [theta_initial_rad; omega_initial];

% --- 2. DEFINE PARAMETERS ---

% Moments of Inertia [Ixx; Iyy; Izz] (Example values in kg*m^2)
% For a 1U cubesat Ixx = Iyy = Izz. 
% I_params = [0.00222; 0.00222; 0.00222] is for a satellite whose 
% mass is 1.33kg with dimensions 0.1m x 0.1m x 0.1m
I_params = [0.00222; 0.00222; 0.00222]; 

% External Torques [Lx; Ly; Lz] 
% Setting L_params to non-zero values to generate the motion seen in your last plot.
%Effective damping coefficient K:
%K = 0.0000025;
mu_r = 40000;
V = 2.83 * 10^-7;
G = 0.25;
mu_0 = 4 * pi * 10^-7;
rho_eff = 1;
B_mag = 40 * 10^-6;
K = calculate_K(mu_r, V, G, mu_0, rho_eff, B_mag);
disp(K);
%L_params = -K * omega_initial; % Example: Torque causing acceleration and spin

% --- 3. SOLVE THE ODE ---
tspan = [0, 3000]; % Time span for integration (0 to 10 seconds)

% Use an anonymous function to pass the fixed parameters to the dynamics function
odefun = @(t, Y) rigid_body_dynamics(t, Y, I_params, K);

% Call the ODE Solver
[T, Y] = ode45(odefun, tspan, Y0);

disp('Initial State Vector Y0 (radians and rad/s):');
disp(Y0);
disp(['ODE solution complete for a time span of ', num2str(tspan(2)), ' seconds.']);

%% Optional: Plotting Results with Angle Wrapping

% Y is the solution matrix: Y(:,1) is theta1, Y(:,2) is theta2, etc.

% --- 1. Angle Wrapping ---
% Convert all angles from RADIANS (solver output) to DEGREES (plotting)
theta1_deg = Y(:,1) * 180/pi;
theta2_deg = Y(:,2) * 180/pi;
theta3_deg = Y(:,3) * 180/pi;

% Apply the modulus function to theta1 to remove total accumulated rotations.
% This wraps the angle into the range [0, 360) degrees, stopping the "shoot up".
theta1_wrapped = mod(theta1_deg, 360); 

%% Plotting Results with Angular Acceleration

% --- 1. CALCULATE ANGULAR ACCELERATION ---
% The dYdt vector contains [dtheta/dt; domega/dt]. We only need domega/dt (rows 4, 5, 6).
% The output of the ODE function is dYdt, so we call it for every time step.
% Pre-allocate the matrix for efficiency
domega_dt = zeros(size(Y, 1), 3); 

% Retrieve fixed parameters used in the ODE call
I_params = [0.00220; 0.00225; 0.00222]; % [Ixx; Iyy; Izz]
%K = 0.0000025; % Damping coefficient K
%K =

for i = 1:size(Y, 1)
    % Call the dynamics function at time T(i) with state Y(i,:)
    dYdt_i = rigid_body_dynamics(T(i), Y(i, :)', I_params, K); 
    % Store d(omega)/dt components (rows 4, 5, 6 of dYdt)
    domega_dt(i, :) = dYdt_i(4:6)'; 
end

% Extract acceleration components for clarity
domega_x = domega_dt(:, 1);
domega_y = domega_dt(:, 2);
domega_z = domega_dt(:, 3);


%% Plotting Results: All on One Window (4 Subplots)

figure(1);

% Convert angles for plotting (done here to keep the subplot section clean)
theta1_deg = Y(:,1) * 180/pi;
theta2_deg = Y(:,2) * 180/pi;
theta3_deg = Y(:,3) * 180/pi;

% Ensure domega_dt was calculated in the previous step
domega_x = domega_dt(:, 1); 
domega_y = domega_dt(:, 2);
domega_z = domega_dt(:, 3);


% --- Top Plot (1/4): Euler Angles - Zoomed-In Focus on Theta 2 & Theta 3 ---
% This shows the small-angle nutation motion clearly.
subplot(4,1,1); % <--- NOW 4 ROWS TOTAL
plot(T, theta2_deg, 'r', T, theta3_deg, 'g'); 
title('Pitch/Roll Motion (\theta_2 and \theta_3) - Zoomed In');
xlabel('Time (s)');
ylabel('Angle (degrees)');
legend('\theta_2 (Pitch)', '\theta_3 (Roll)');
ylim([-30, 30]); % <--- Tight Y-LIMITS APPLIED HERE

% --- Second Plot (2/4): Euler Angle - Zoomed-Out Yaw (Theta 1) ---
% This shows the large, continuous rotation.
subplot(4,1,2);
plot(T, theta1_deg, 'b'); 
title('Yaw Rotation (\theta_1) - Zoomed Out');
xlabel('Time (s)');
ylabel('Angle (degrees)');
legend('\theta_1 (Yaw)');

% --- Third Plot (3/4): Rotational Velocities ---
subplot(4,1,3);
plot(T, Y(:,4), 'b', T, Y(:,5), 'r', T, Y(:,6), 'g'); 
title('Rotational Velocities Over Time');
xlabel('Time (s)');
ylabel('Velocity (rad/s)');
legend('\omega_x', '\omega_y', '\omega_z');

% --- Bottom Plot (4/4): Angular Acceleration ---
subplot(4,1,4); 
plot(T, domega_x, 'b', T, domega_y, 'r', T, domega_z, 'g');
title('Angular Acceleration Over Time');
xlabel('Time (s)');
ylabel('Acceleration (rad/s^2)');
legend('d\omega_x/dt', 'd\omega_y/dt', 'd\omega_z/dt');

%% =========================================================================
% LOCAL FUNCTION: The ODE System (dY/dt)
% This function is called by ode45. It must take (t, Y, ...) and return dYdt.
% =========================================================================

function dYdt = rigid_body_dynamics(t, Y, I_params, K)
    % ---------------------------------------------------------------------
    % DESCRIPTION: Calculates the time derivatives (dYdt) for the 6-DOF 
    % rigid body system, including dynamic hysteresis damping torque.
    % ---------------------------------------------------------------------
    % === 1. UNPACK STATE VECTOR (Y) ===
    % Y = [theta1; theta2; theta3; omega_x; omega_y; omega_z]
    theta1  = Y(1); % Yaw angle
    theta2  = Y(2); % Pitch angle
    theta3  = Y(3); % Roll angle
    omega_x = Y(4); % Angular velocity component about body x-axis
    omega_y = Y(5); % Angular velocity component about body y-axis
    omega_z = Y(6); % Angular velocity component about body z-axis
    
    % === 2. UNPACK INERTIA PARAMETERS ===
    Ixx = I_params(1); % Moment of Inertia about body x-axis
    Iyy = I_params(2); % Moment of Inertia about body y-axis
    Izz = I_params(3); % Moment of Inertia about body z-axis
    
    % === 3. DYNAMIC TORQUE CALCULATION (THIS IS THE FIX) ===
    % The external torque L (Lx, Ly, Lz) is the Damping Torque L = -K * omega.
    % This ensures the torque always opposes the current motion.
    Lx = -K * omega_x; 
    Ly = -K * omega_y;
    Lz = -K * omega_z;
    
    % === 4. PRE-CALCULATE TRIGONOMETRIC TERMS ===
    cos2 = cos(theta2);
    sin2 = sin(theta2);
    cos3 = cos(theta3);
    sin3 = sin(theta3);
    
    % === 5. KINEMATICS (Euler Angle Derivatives) ===
    
    % Handle singularity (Gimbal Lock) when cos(theta2) is near zero
    if abs(cos2) < 1e-6
        dtheta1dt = NaN; 
        dtheta3dt = NaN;
    else
        % d(theta1)/dt (Yaw rate)
        dtheta1dt = (omega_y * sin3 + omega_z * cos3) / cos2;
        
        % d(theta3)/dt (Roll rate)
        dtheta3dt = omega_x + (sin3 * sin2 / cos2) * omega_y + ...
                            (cos3 * sin2 / cos2) * omega_z;
    end
    
    % d(theta2)/dt (Pitch rate)
    dtheta2dt = omega_y * cos3 - omega_z * sin3;
    
    % === 6. DYNAMICS (Rotational Velocity Derivatives - Euler's Equations) ===
    % The dynamic Lx, Ly, Lz calculated in Step 3 are used here.
    
    % d(omega_x)/dt
    domega_x_dt = (-(Izz - Iyy) * omega_y * omega_z + Lx) / Ixx;
    
    % d(omega_y)/dt
    domega_y_dt = (-(Ixx - Izz) * omega_z * omega_x + Ly) / Iyy;
    
    % d(omega_z)/dt
    domega_z_dt = (-(Iyy - Ixx) * omega_x * omega_y + Lz) / Izz;
    
    % === 7. PACK DERIVATIVES ===
    % The output must be a single 6x1 vector [d(theta)/dt; d(omega)/dt]
    dYdt = [dtheta1dt; dtheta2dt; dtheta3dt; domega_x_dt; domega_y_dt; domega_z_dt];
end

function K = calculate_K(mu_r, V, G, mu_0, rho_eff, B_mag)
    % ---------------------------------------------------------------------
    % DESCRIPTION: Calculates the linear damping coefficient K based on the 
    % rod's physical properties and the LEO magnetic field.
    % ---------------------------------------------------------------------

    % --- STEP 1: Calculate Rod Constant A ---
    % A = (pi^2 * mu_r * V * G) / (2 * mu_0 * rho_eff)
    A = (pi^2 * mu_r * V * G) / (2 * mu_0 * rho_eff);
    
    % --- STEP 2: Calculate Damping Coefficient K ---
    % K = A * |B|^2
    K = A * B_mag^2;
    
    % The K value is the primary output of this function
end