function demo_two_body_rotating_frame
% Define parameters (use consistent units!)
p.Mu_Earth = 398600.4418;      % km^3/s^2
p.Mu_Luna  = 4902.800066;      % km^3/s^2
R_EM       = 384400;           % km, Earth–Moon distance
p.d_Earth  = -4.6707e+03;                % x-position of Earth in this frame (km)
p.d_Luna   = R_EM+p.d_Earth;             % x-position of Moon in this frame (km)

% If you're in a synodic (rotating) frame, Omega is the mean motion:
p.Omega    = sqrt( (p.Mu_Earth + p.Mu_Luna) / R_EM^3 )  % rad/s

% Initial conditions [x y z vx vy vz] (example: near Earth on x-axis)
X0 = [300000, 0, 0, 0, 0, 0]; % km, km/s (choose your own ICs)

% Integration window
tspan = [0, 2.3570e+06];   

% Solver options
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

% Integrate
[t,X] = ode45(@(t,X) eom(t,X,p), tspan, X0, opts);

% Plot trajectory projection
figure;
hold on;
plot(X(:,1), X(:,2)); axis equal; grid on
plot(X(1,1), X(1,2),'ob');                             % S/C start point
xlabel('x (km)'); ylabel('y (km)'); title('Trajectory in rotating frame');
hold on; plot(p.d_Earth,0,'bo','MarkerFaceColor','b'); % Earth
plot(p.d_Luna,0,'ko','MarkerFaceColor','k');           % Moon
legend('Trajectory','Earth','Moon');
end

function dXdt = eom(~,X,p)
% Unpack state
x = X(1); y = X(2); z = X(3);
vx = X(4); vy = X(5); vz = X(6);

% Distances to Earth and Moon
rE = sqrt( (x - p.d_Earth)^2 + y^2 + z^2 );
rM = sqrt( (x - p.d_Luna )^2 + y^2 + z^2 );

% Accelerations (your formulas)
ax =  2*p.Omega*vy + (p.Omega^2)*x ...
    - (p.Mu_Earth/(rE^3))*(x - p.d_Earth) ...
    - (p.Mu_Luna /(rM^3))*(x - p.d_Luna);

ay = - (p.Mu_Earth/(rE^3))*y ...
     - (p.Mu_Luna /(rM^3))*y ...
     + (p.Omega^2)*y ...
     - 2*p.Omega*vx;

az = - (p.Mu_Earth/(rE^3))*z ...
     - (p.Mu_Luna /(rM^3))*z;

% State derivative
dXdt = [vx; vy; vz; ax; ay; az];
end
