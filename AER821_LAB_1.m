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

% Motion and period
Omega = sqrt(G * (M_Earth + M_Luna) / r_12^3); % angular velocity [rad/s]
T_sys = 2 * pi / Omega; % lunar period [s]

% Lagrange Point L1
mu = M_Luna/(M_Earth + M_Luna);
d_from_moon = r_12*(mu/3)^(1/3);
x_guess = d_Luna - d_from_moon;

L1_fun = @(x) (Mu_Earth/(abs(x - d_Earth)^2) - Mu_Luna/(abs(d_Luna - x)^2) - Omega^2 * x);
xL1 = fzero(L1_fun, x_guess);

% Initial position at L1
r0x = xL1; 
r0y = 0;
r0z = 0;

% Initial velocity at L1
v0x = 0;
v0y = Omega * xL1;
v0z = 0;

% State vector
r0 = [r0x; r0y; r0z];
v0 = [v0x; v0y; v0z];
y0 = [r0; v0];

% Time span: one lunar cycle
tspan = [0 T_sys];

% Solve system
[t, y] = ode113(@(t,y) crtbp_inertial(t, y, Mu_Earth, Mu_Luna, d_Earth, d_Luna, Omega), tspan, y0);

% Smooth plot lines
tt = linspace(0, T_sys, 2000).';
R_E = [d_Earth * cos(Omega * tt), d_Earth * sin(Omega * tt)];
R_M = [d_Luna * cos(Omega * tt), d_Luna * sin(Omega * tt)];

% Inertial-frame plot
figure;
hold on;
plot(0,0,'k+','LineWidth',1.5); % barycenter
plot(R_E(:,1), R_E(:,2), 'b','LineWidth',1.5); % Earth path
plot(R_M(:,1), R_M(:,2), 'k','LineWidth',1.5); % Moon path
plot(y(:,1), y(:,2), 'r','LineWidth',1.2); % spacecraft path
plot(d_Earth,0,'bo','MarkerFaceColor','b'); % Earth at t0
plot(d_Luna, 0,'ko','MarkerFaceColor','k'); % Moon at t0
plot(xL1,0,'ro','MarkerFaceColor','r'); % L1 
plot(y(end,1),y(end,2),'rx','LineWidth',2,'MarkerSize',8); % spacecraft end
legend('Barycenter','Earth orbit','Moon orbit','Spacecraft','Earth(t0)','Moon(t0)','L1','Spacecraft end');
xlabel('X [m]'); ylabel('Y [m]');
title('CRTBP – One Lunar Period (Start at L1)');
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
    r = y(1:3); % spacecraft position [m]
    v = y(4:6); % spacecraft velocity [m/s]

    % Earth and Moon positions
    R_Earth = [d_Earth * cos(Omega * t); d_Earth * sin(Omega * t); 0]; 
    R_Luna  = [d_Luna  * cos(Omega * t); d_Luna  * sin(Omega * t); 0]; 

    % Relative positions
    r1 = r - R_Earth; 
    r2 = r - R_Luna;  

    % Accelerations
    a = -Mu_Earth * r1 / norm(r1)^3 - Mu_Luna * r2 / norm(r2)^3; 

    % Derivative
    dydt = [v; a];

end
