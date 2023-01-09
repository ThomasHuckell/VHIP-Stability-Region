% physical parameters
g           = 9.81;     % gravitational accelereation [m/s/s]
z0          = 0.95;        % resting CoM height [m]
m           = 70;       % mass [Kg]

tauMax      = 100;     % maximum flywheel torque [Nm]
thetaMax    = 0.5;       % maximum flywheel angular rotation [rad]
J           = 8;      % rotational inertia of the Flywheel

zmin        = 0.9;
zmax        = 1;
DDz_max     = 3;

delta       = 0.2;     % base of support size
% controller parameters
Kp  = 3;       % proportional gain
Kd  = 1;       % deriviative gain

% Store variables in structure
modelParam.gravity          = g;
modelParam.restHeight       = z0;
modelParam.mass             = m;
modelParam.supportSize      = delta;
modelParam.pGain            = Kp;
modelParam.dGain            = Kd;
modelParam.inertia          = J;
modelParam.theta            = thetaMax;
modelParam.zbound           = [zmin zmax];
modelParam.delta            = delta;
modelParam.DDz              = DDz_max;

omega = sqrt(g/z0);

% maximum time for the bang bang controller s.t theta <= theta_max
T = sqrt(J*thetaMax/tauMax);