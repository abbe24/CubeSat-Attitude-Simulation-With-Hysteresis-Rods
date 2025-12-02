%% run_PMAC_sim.m
% Passive Magnetic Attitude Control (PMAC) Simulation
% Based on CSSWE (Gerhardt & Palo 2010)

clear; clc;

%% Physical constants
mu0   = 4*pi*1e-7;     % [T·m/A] permeability of free space
Re    = 6378e3;        % [m] Earth radius
alt   = 600e3;         % [m] assumed 600 km orbit
r0    = Re + alt;
muE   = 3.986004418e14;

% Mean motion
n0    = sqrt(muE / r0^3);

%% Orbit / field model parameters
params.Heq  = 18.3;          % [A/m] equatorial field
params.incl = deg2rad(55);   % [rad]
params.n0   = n0;            % [rad/s]
params.u0   = 0;             % [rad]
params.mu0  = mu0;

%% Spacecraft inertias (CSSWE)
params.Ixx = 0.00551;
params.Iyy = 0.02552;
params.Izz = 0.02565;

%% Bar magnet
params.m_bar = [0.3; 0; 0];  % [A·m^2] permanent magnet along body X

%% Hysteresis rods — YOUR rods (1/16" × 5")
D  = 1.5875e-3;   % [m] diameter 1/16"
L  = 0.127;       % [m] length 5"
Vrod = pi*(D/2)^2 * L;

params.Vrod   = Vrod;
params.Hc_app = 12;      % [A/m] apparent coercive force
params.Br_app = 0.004;   % [T]
params.Bs_app = 0.027;   % [T]

% Rod count
params.nRodsY = 1;
params.nRodsZ = 1;

%% Initial conditions
theta0 = deg2rad([0; 0; 0]);         % Euler angles
w0     = deg2rad([10; 5; 5]);        % body rates [deg/s]

x0 = [theta0; w0];

%% Simulation time
tEnd = 10*24*3600;   % 10 days
opts = odeset('RelTol',1e-7,'AbsTol',1e-7);

%% Integrate
[t, x] = ode45(@(t,x) PMAC_ODE(t,x,params), [0 tEnd], x0, opts);

theta = x(:,1:3);
w     = x(:,4:6);

%% Compute beta(t) — angle between body X-axis and B-field
beta = zeros(size(t));
for k = 1:length(t)
    B_inertial = magFieldDipole(t(k), params);
    C_bi = DCM321(theta(k,1), theta(k,2), theta(k,3));
    B_body = C_bi * B_inertial(:);

    ex = [1;0;0];
    beta(k) = acos( dot(ex,B_body) / (norm(B_body)+eps) );
end
beta_deg = rad2deg(beta);

%% Plots
figure; plot(t/86400, rad2deg(w));
xlabel('Days'); ylabel('\omega [deg/s]');
legend('\omega_x','\omega_y','\omega_z'); grid on;

figure; plot(t/86400, beta_deg);
xlabel('Days'); ylabel('\beta [deg]');
grid on; title('Angle Between Body X and Magnetic Field');
%% === Export Attitude to STK ===

% Convert Euler angles (radians) → degrees
yaw_deg   = rad2deg(theta(:,1));
pitch_deg = rad2deg(theta(:,2));
roll_deg  = rad2deg(theta(:,3));

% Build absolute UTCG timestamps (or fake one if t_datetime is missing)
if exist('t_datetime','var')
    % Real STK time available
    absoluteTimes = t_datetime(1) + seconds(t);
else
    % No STK calendar time → use a fake but valid start time
    absoluteTimes = datetime(2025,1,1,0,0,0) + seconds(t);
end

% Convert timestamps into STK-style text
timeStrings = datestr(absoluteTimes, 'dd mmm yyyy HH:MM:SS.FFF');

% Build table for STK
A = table(string(timeStrings), yaw_deg, pitch_deg, roll_deg, ...
          'VariableNames', ["Time_UTCG","Yaw_deg","Pitch_deg","Roll_deg"]);

% Save to CSV
writetable(A, "CubeSat_Attitude_321.csv");
%% === Export STK .a Attitude File ===

Npts = numel(t);  % number of time points

% Use the same absoluteTimes we built for the CSV.
% If it doesn't exist for some reason, fake an epoch.
if exist('absoluteTimes','var')
    epoch = absoluteTimes(1);
else
    epoch = datetime(2025,1,1,0,0,0);
end

% Text version of epoch in STK style
epochStr = datestr(epoch,'dd mmm yyyy HH:MM:SS.FFF000000');

% Time in seconds from epoch (STK EpSec)
t_sec = seconds( (epoch + seconds(t) - epoch) );  % basically just t

fid = fopen('CubeSat_Attitude_321.a','w');

% ----- Header -----
fprintf(fid,'stk.v.12.2\n');
fprintf(fid,'BEGIN Attitude\n');
fprintf(fid,'NumberOfAttitudePoints %d\n', Npts);
fprintf(fid,'ScenarioEpoch %s\n', epochStr);
fprintf(fid,'CentralBody Earth\n');
fprintf(fid,'CoordinateAxes J2000\n');
fprintf(fid,'Sequence 321\n');
fprintf(fid,'TimeFormat EpSec\n\n');
fprintf(fid,'AttitudeTimeEulerAngles\n');

% ----- Data: time(sec)  yaw  pitch  roll (all degrees) -----
for k = 1:Npts
    fprintf(fid,'%.6f %.8f %.8f %.8f\n', ...
        t_sec(k), yaw_deg(k), pitch_deg(k), roll_deg(k));
end

fprintf(fid,'END Attitude\n');
fclose(fid);
