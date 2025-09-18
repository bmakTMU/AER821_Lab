%% Question 3  — Rotating-frame scalar equations (direct integration)

% Initial state in the rotating frame.
% At L1, the equilibrium has zero synodic velocity:
x0 = xL1; y0 = 0; z0 = 0;
xd0 = 0;  yd0 = 0; zd0 = 0;
s0 = [x0; y0; z0; xd0; yd0; zd0];

% Integrate for the same span as Q2 (5 lunar periods is a nice view)
tspan3 = [0 5*T_sys];

% Tighter tolerances help near equilibria
opts = odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',T_sys/1000);

[t3, s3] = ode113(@(t3,s3) crtbp_rotating_scalar( ...
                    t3, s3, Mu_Earth, Mu_Luna, d_Earth, d_Luna, Omega), ...
                  tspan3, s0, opts);

x3 = s3(:,1);  y3 = s3(:,2);

% --- Plot in the rotating frame, Earth and Moon fixed ---
figure; hold on; grid on; axis equal;
plot(d_Earth,0,'bo','MarkerFaceColor','b');           % Earth
plot(d_Luna, 0,'ko','MarkerFaceColor','k');           % Moon
plot(0,0,'k+');                                       % Barycenter
plot(xL1,0,'ro','MarkerFaceColor','r');               % L1 marker
plot(x3, y3, 'g--','LineWidth',1.4);                  % Q3 trajectory

% If Q2 variables exist, overlay for visual agreement
if exist('X_SC_Rot','var') && exist('Y_SC_Rot','var')
    plot(X_SC_Rot, Y_SC_Rot, 'r-','LineWidth',1.1);
    legend('Earth','Moon','Barycenter','L1', ...
           'Q3: Rotating ODE','Q2: Inertial→Rotating','Location','best');
else
    legend('Earth','Moon','Barycenter','L1','Q3: Rotating ODE','Location','best');
end
xlabel('X [m]'); ylabel('Y [m]');
title('CRTBP — Rotating Frame (Direct Scalar Equations)');
hold off;

function ds = crtbp_rotating_scalar(~, s, Mu1, Mu2, x1, x2, Omega)
% State s = [x; y; z; xd; yd; zd] in the rotating (synodic) frame
    x  = s(1);  y  = s(2);  z  = s(3);
    xd = s(4);  yd = s(5);  zd = s(6);

    % Distances to the primaries (fixed at x1=d_Earth and x2=d_Luna)
    r1 = sqrt((x - x1)^2 + y^2 + z^2);
    r2 = sqrt((x - x2)^2 + y^2 + z^2);

    % Scalar CR3BP equations in rotating frame
    xdd =  2*Omega*yd + Omega^2*x  -Mu1*(x - x1)/r1^3 - Mu2*(x - x2)/r2^3;

    ydd = -2*Omega*xd + Omega^2*y - Mu1*y/r1^3       - Mu2*y/r2^3;

    zdd = - Mu1*z/r1^3        - Mu2*z/r2^3;

    ds = [xd; yd; zd; xdd; ydd; zdd];
end