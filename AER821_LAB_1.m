clear; clc; close all;

% Constants
G = 6.6742e-11; % gravitational constant [m^3/kg/s^2]
M_Earth = 5.974e24; % mass of Earth [kg]
M_Luna  = 7.348e22; % mass of Moon [kg]
r_12    = 384400e3; % distance between Earth and Moon [m]

% Gravitational parameters
Mu_Earth = G * M_Earth; % [m^3/s^2]
Mu_Luna  = G * M_Luna;  % [m^3/s^2]

% Distances of Earth and Moon from barycenter
Pi_Earth = M_Earth / (M_Earth + M_Luna);
Pi_Luna  = M_Luna  / (M_Earth + M_Luna);
d_Earth  = -Pi_Luna * r_12; % Earth position offset [m]
d_Luna   =  Pi_Earth * r_12; % Moon position offset [m]

%% Question 1

Omega = sqrt(G * (M_Earth + M_Luna) / r_12^3); % angular velocity [1/s]
T_sys = 2 * pi / Omega; % lunar period [s]

% Spacecraft initial conditions
Re = 6378e3; % Earth radius [m]
h  = 35786e3; % altitude of spacecraft [m] (GEO)

% At this altitude, a satellites orbital period matches Earth s rotation (24 hours)

r_sc_mag = Re + h; % distance from Earth center [m]

% Start position: spacecraft on +x axis
r0x = d_Earth + r_sc_mag; % [m]
r0y = 0; % [m]
r0z = 0; % [m]

% Start velocity: circular orbit speed, in +y
v0x = 0; % [m/s]
v0y = sqrt(Mu_Earth / r_sc_mag); % [m/s]
v0z = 0;% [m/s]

r0 = [r0x; r0y; r0z];  % make position vector [m]
v0 = [v0x; v0y; v0z];  % make velocity vector [m/s]

% Initial state
y0 = [r0; v0];

% Time span: one lunar cycle
tspan = [0 T_sys];

% Solve system
[t, y] = ode113(@(t,y) crtbp_inertial(t, y, Mu_Earth, Mu_Luna, d_Earth, d_Luna, Omega), tspan, y0);

% Earth and Moon positions
R_E = [d_Earth * cos(Omega * t), d_Earth * sin(Omega * t), zeros(size(t))];
R_M = [d_Luna  * cos(Omega * t), d_Luna  * sin(Omega * t), zeros(size(t))];

% Spacecraft relative to Earth (for orbit view)
rel = y(:,1:2) - R_E(:,1:2);

% Plot: Inertial (barycenter frame)
figure;
hold on;
plot(R_E(:,1), R_E(:,2), 'b', 'LineWidth', 2); % Earth
plot(R_M(:,1), R_M(:,2), 'k', 'LineWidth', 2); % Moon
plot(y(:,1), y(:,2), 'r'); % spacecraft
legend('Earth','Moon','Spacecraft');
xlabel('X [m]'); ylabel('Y [m]');
title('CRTBP Inertial Frame - One Lunar Period');
axis equal; grid on;
hold off;

% Run Simulink model
simOut = sim('AER821_LAB_1_S', 'StopTime', num2str(T_sys));

% Get x position
x_sim = squeeze(simOut.x_sc.Data);

% Get y position
y_sim = squeeze(simOut.y_sc.Data);

% Get time
t_sim = simOut.x_sc.Time;

% Plot: Trajectory from Simulink
figure;
hold on;
plot(x_sim, y_sim, 'r');
xlabel('X [m]'); ylabel('Y [m]');
title('Spacecraft trajectory from Simulink');
axis equal; grid on;
hold off;

% Function: Equations of motion
function dydt = crtbp_inertial(t, y, Mu_Earth, Mu_Luna, d_Earth, d_Luna, Omega)

    % Position and velocity
    r = y(1:3); % elements 1–3 = [rx, ry, rz], spacecraft position [m]
    v = y(4:6); % elements 4–6 = [vx, vy, vz], spacecraft velocity [m/s]

    % Earth and Moon positions
    R_Earth = [d_Earth * cos(Omega * t); d_Earth * sin(Omega * t); 0]; % [m]
    R_Luna  = [d_Luna  * cos(Omega * t); d_Luna  * sin(Omega * t); 0]; % [m]

    % Relative positions
    r1 = r - R_Earth; % vector Earth to spacecraft [m]
    r2 = r - R_Luna; % vector Moon to spacecraft [m]

    % Accelerations
    a = -Mu_Earth * r1 / norm(r1)^3 - Mu_Luna * r2 / norm(r2)^3; % acceleration [m/s^2]

    % Derivative
    dydt = [v; a];

end
