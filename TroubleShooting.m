syms x(t) omega m z0 tau delta T x_0 xdot_0 z DDz g
Dx = diff(x);

ode1 = diff(x,t,2) == omega^2*(x-delta) - 1/(z0*m)*tau;

cond1 = x(0) == x_0;
cond2 = Dx(0) == xdot_0;

conds = [cond1 cond2];

xSol(t) = dsolve(ode1,conds);

xSol = simplify(xSol);

DxSol = diff(xSol,t);

DxSol = simplify(DxSol);

ode2 = diff(x,t,2) == omega^2*(x-delta) + 1/(z0*m)*tau;

cond1 = x(0) == xSol(T);
cond2 = Dx(0) == DxSol(T);

conds = [cond1 cond2];

xFin(t) = dsolve(ode2,conds);
xFin(t) = simplify(xFin);

DxFin(t) = diff(xFin,t);


CP_end = matlabFunction(simplify(xFin(T) + DxFin(T)/omega - delta));

%%
Dz = diff(z);

ode3 = diff(x,t,2) == g/z*(x-delta);
ode4 = diff(z,t,2) == DDz;

cond1 = x(0) == x_0;
cond2 = Dx(0) == xdot_0;
cond3 = z(0) == z0;
cond4 = Dz(0) == 0;


conds = [cond1 cond2;
         cond3 cond4];



%%
% physical parameters
g           = 9.81;     % gravitational accelereation [m/s/s]
z0          = 0.9;        % resting CoM height [m]
m           = 70;       % mass [Kg]

delta       = 0.2;
tauMax      = 100;     % maximum flywheel torque [Nm]
thetaMax    = 0.5;       % maximum flywheel angular rotation [rad]
J           = 8;      % rotational inertia of the Flywheel

omega = sqrt(g/z0);
T = sqrt(J*thetaMax/tauMax);



CP_end(T,delta,m,omega,tauMax,x_0,xdot_0,z0)
