%% Symbolic derivation for the boundry of the region of stability of the 
% LIP + Flywheel model
% Thomas Huckell Jan 9 2021


% define symbolic variables 

% x0 : initial CoM position
% xdot0 : initial CoM Velocity
% w : omega is the natural frequency of the system w = sqrt(g/z)
% d : is the maximum CoP location i.e the length of the base of support
% tau : is maximuim torque of the flywheel
% m : is the mass of the LIP
% z : is the constant height the LIP is constrained to
% T : is the time used for the bang bang controller
% s : is the symbolic frequency domain variable
% t : is the symbolic time domain variable 

syms x0 xdot0 w d tau m z T s t 


%% Sovling the Zero Input Response

% Laplace transform of the dynamics with tau = 0
X_zir(s) = (s^2*x0 + s*xdot0 - w^2*d)/(s^3-s*w^2);

X_zir(s) = partfrac(X_zir,s);

x_zir(t) = ilaplace(X_zir(s));

Dx_zir(t) = diff(x_zir(t),t);

% Manually simplified
x_zirM(t) = d + ( x0 - d)*cosh(w*t) + xdot0/w*sinh(w*t); %  Confirmed x_zir-x_zirM = 0 
Dx_zirM(t) = w*(x0 - d)*sinh(w*t) + xdot0*cosh(w*t);     %  Confirmed x_zir-x_zirM = 0 

simplifyCheck1 = simplify(x_zir - x_zirM);
simplifyCheck2 = simplify(Dx_zir - Dx_zirM);

if simplifyCheck1 == 0 && simplifyCheck2 == 0
    disp('Zero Input Response simplification valid')
end

%% Solving Zero State Response

% Laplace transform of the dynamics with x0 = 0 and xdot0 = 0
X_zsr(s) = (-d*w^2 -1/(m*z)*tau  +1/(m*z)*2*tau*exp(-T*s)  -1/(m*z)*tau*exp(-2*s*T))/(s^3-s*w^2);

assume(T > 0);

x_zsr(t) = ilaplace(X_zsr(s));

Dx_zsr(t) = diff(x_zsr(t),t);

% % Manually simplified
x_zsrM(t) = -tau/(m*w^2*z)*((cosh(w*t)-1)*heaviside(t) - 2*(cosh(w*(t-T))- 1)*heaviside(t-T) ...
                 + (cosh(w*(t-2*T)) - 1)*heaviside(t-2*T))+ d*(1-cosh(t*w))
             
Dx_zsrM(t) = -tau/(m*w*z)*(sinh(w*t)*heaviside(t) - 2*sinh(w*(t-T))*heaviside(t-T) ...
                + sinh(w*(t-2*T))*heaviside(t-2*T)) - d*w*sinh(t*w)



x(t) = x_zirM + x_zsrM;

Dx(t) = Dx_zirM + Dx_zsrM;

xfin = simplify(x(2*T))

Dxfin = simplify(Dx(2*T))

% xfin = rewrite(x(2*T),'sinhcosh')
% Dxfin = rewrite(Dx(2*T),'sinhcosh')

cp = simplify(xfin + Dxfin/w);



