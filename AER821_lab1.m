%Values taken from textbook

% Mass values, kg
M_earth = 5.974 * 10^24; 
M_moon = 73.48 * 10^21;
M_total = M_earth + M_moon; %combined mass of earth and moon

% Gravitational Parameters

G = 6.6743 * 10^-20; % Gravitational Constant km^3 / kg*s^2

mu = G*M_total;
mu_earth = G*M_earth;
mu_moon = G*M_moon;

% dimensionless mass ratios, pi
pi_earth = M_earth / M_total;
pi_moon = M_moon / M_total;

r_12 = 3.844*10^5;% distance between earth and moon (from CoM), km
