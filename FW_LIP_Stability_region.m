%% Stability region of existing models.
% Thomas Huckell January 22 2021


%% Variable decleration
close all
clc
clear

addpath("C:\Users\tomhu\OneDrive - Queen's University\Documents\fileExchange\arrow3")
ltxFMT = {'Interpreter','latex'};
spec = plotSpec();


% physical parameters
g           = 9.81;     % gravitational accelereation [m/s/s]
z0          = 0.9;        % resting CoM height [m]
m           = 70;       % mass [Kg]

tauMax      = 100;     % maximum flywheel torque [Nm]
thetaMax    = 0.5;       % maximum flywheel angular rotation [rad]
J           = 8;      % rotational inertia of the Flywheel

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


omega = sqrt(g/z0);

% maximum time for the bang bang controller s.t theta <= theta_max
T = sqrt(J*thetaMax/tauMax);


%% The maximum recoverable capture points uncomment one of the different cp_max

%Region from Stephens Humanoid push recovery

%cp_max = delta + tauMax/(m*g)*(exp(omega*T) - 1)^2;


% Region using stephens general form but pratts leading coefficient on the
% torque term

%cp_max = delta + tauMax/(m*g)*(exp(-2*omega*T)*(exp(omega*T) - 1)^2);

% pratts region

cp_max = delta + tauMax/(m*g)*(1 - exp(-omega*T))^2; %exp(omega*2*T)-2*exp(T*omega)+1)/exp(omega*2*T)



% Region that I derived

%cp_max = (delta -2*delta*(1-exp(2*T*omega)) + tauMax/(omega^2*m*z0)*(1-exp(T*omega))^2)/exp(2*T*omega);

% troubleShoot

%cp_max = delta  + tauMax/(m*g)*(exp(-omega*T)-1)^2-delta*(exp(-2*omega*T)-1);

%% LIP Simulations

% Initialize grid of inital conditions for to simulate LIP balancing
% trajecotries

xdotArray = linspace(-.95,.95,15);
xArray    = linspace(-0.19,0.19,15);

% Define boundry of LIP stability region 
t    = linspace(-2,2,20);

LIP_UB = (-t + delta)*omega;
LIP_LB = (-t - delta)*omega;
eigen = (-t)*omega;

% Plot the trajectory of LIP model with PD ankle control
fig1 = figure('Position',[100 100 750 300]);
xlim([-0.2 0.2])
ylim([-1 1])
daspect([0.1,1,1])

hold all
c = 1;




for i = 1:length(xdotArray)
    for j = 1:length(xArray)
        
        % Set intitial conditions
        X0 = [xArray(j); xdotArray(i)];
        t0 = 0;
        
        % Simulate LIP model from inital conditions
        [t_LIP, X_LIP] = ode45( @(t,X) LIP(t,X,modelParam) , [t0,3], X0);
        
        
        % Plot the phase of trajectory
        if c == 1
            
            plot(X_LIP(:,1),X_LIP(:,2),'Color','#A2142F','HandleVisibility','off');
            plot(X_LIP(1,1),X_LIP(1,2),'ko','MarkerSize',2);
            plot(t,LIP_UB,'Color',spec.LIP{2},'LineWidth',2);
            plot(t,LIP_LB,'Color',spec.LIP{2},'LineWidth',2,'HandleVisibility','off');
            plot(t,eigen,'--k','LineWidth',2)

            
            if abs(X0(1)+X0(2)/omega)<= delta
                addArrow(X_LIP(:,1),X_LIP(:,2),[0.1,0.6],1.9,'_b');
            else
                addArrow(X_LIP(:,1),X_LIP(:,2),[0.1,0.6],1.9,'_r')
            end
        else
            
            if abs(X0(1)+X0(2)/omega)<= delta
                plot(X_LIP(:,1),X_LIP(:,2),'Color','#000097','HandleVisibility','off');
                addArrow(X_LIP(:,1),X_LIP(:,2),[0.1,0.6],1.9,'_b');
            else
                plot(X_LIP(:,1),X_LIP(:,2),'Color','#A2142F','HandleVisibility','off');
                addArrow(X_LIP(:,1),X_LIP(:,2),[0.1,0.6],1.9,'_r')
            end
            plot(X_LIP(1,1),X_LIP(1,2),'ko','MarkerSize',2,'HandleVisibility','off');
            
        end
        
        c = c + 1;
    end
end




% color fill regions on the plot

v1 = [[t,fliplr(t)]',[LIP_UB,fliplr(LIP_LB)]'];
f1 = 1:1:length([t,fliplr(t)]);

patch('vertices',v1,'faces',f1,'FaceColor',spec.LIP{2},'FaceAlpha',0.15,'edgeColor','none')




x1 = [t(1) t(end) t(end)];
y1 = [LIP_UB(1) LIP_UB(end) 1];

x2 = [t(1) t(1) t(end)];
y2 = [-1 LIP_LB(1) LIP_LB(end)];



patch(x1,y1,'r','FaceAlpha',0.15,'edgeColor','none')
patch(x2,y2,'r','FaceAlpha',0.15,'edgeColor','none')

xline(0,':');
yline(0,':');

hold off


xlabel('$x$ [m]',ltxFMT{:})
ylabel('$\dot{x}$ [m/s]',ltxFMT{:})
legend('Initial Condition','$|x+\frac{\dot{x}}{\omega}| = \delta$','$x+\frac{\dot{x}}{\omega}=0$','Location','northwest','Interpreter','latex')

exportgraphics(fig1,'figures/LIP_stability.png','Resolution',300)

%% Flywheel Simulations


% LIP+FW region of stability
FW_UB = (-t + cp_max)*omega;
FW_LB = (-t - cp_max)*omega;


% Set the initial velocity for the FW simulations. These simulations will
% have the intial position be zero.
xdot1 = .5*delta*omega;
xdot2 = delta*omega;
xdot3 = (0.1*delta*omega + 0.9*cp_max*omega);
xdot4 = cp_max*omega+0.01;

xdotArray = [xdot1 xdot2 xdot3 xdot4];

colorArray = {'_a','^a','p','f';
              '#2E6815', '#93BB6D', '#BFC01B', '#B32121'};

% plot the trajectories of LIP+FW starting at the set inital velocities

fig2 = figure('Position',[100 100 750 300]);
daspect([0.1,1,1]);
xlim([-0.2 0.20])
ylim([-1 1])

fig3 = figure('Position',[100 100 750 450]);


for i = 1:length(xdotArray)
    
    % Set the inital conditions
    X0 = [0;0; xdotArray(i);0];
    
    % simulate the LIP + FW
    [t_FW, X_FW , U_FW] = simFW(tauMax,X0,modelParam,cp_max);
    
    figure(fig2)
    hold all
    % plot the trajectory
    if i == 1
        
        plot(X_FW(:,1),X_FW(:,3),'Color',colorArray{2,i},'LineWidth',1.5,'HandleVisibility','off');
        plot(X_FW(1,1),X_FW(1,3),'ko','MarkerSize',2);
        plot(X_FW(end,1),X_FW(end,3),'kx')
        plot(t,eigen,'--k')
        % plot the different regions of stability
        plot(t,LIP_UB,'Color',spec.LIP{2},'HandleVisibility','on','LineWidth',2)
        plot(t,FW_UB,'Color',spec.LIPPFW{2},'HandleVisibility','on','LineWidth',2)
        plot(t,FW_LB,'Color',spec.LIPPFW{2},'HandleVisibility','off','LineWidth',2)
        plot(t,LIP_LB,'Color',spec.LIP{2},'HandleVisibility','off','LineWidth',2)
        plot(t,eigen,'--k','HandleVisibility','off')
        
        addArrow(X_FW(:,1),X_FW(:,3),[0.3,0.8],1.9,colorArray{1,i});
        
        
    else
        
        plot(X_FW(:,1),X_FW(:,3),'Color',colorArray{2,i},'HandleVisibility','off','LineWidth',1.5);
        plot(X_FW(1,1),X_FW(1,3),'ko','MarkerSize',2,'HandleVisibility','off');
        plot(X_FW(end,1),X_FW(end,3),'kx','HandleVisibility','off');
        
        addArrow(X_FW(:,1),X_FW(:,3),[0.3,0.8],1.9,colorArray{1,i});
    end
    
    hold off
    
    figure(fig3);
    subplot(3,1,1)
    hold on
    plot(t_FW,X_FW(:,2),'Color',colorArray{2,i},'LineWidth',1.5)
    subplot(3,1,2)
    hold on
    plot(t_FW,X_FW(:,4),'Color',colorArray{2,i},'LineWidth',1.5)
    subplot(3,1,3)
    hold on
    plot(U_FW(:,1),U_FW(:,2),'Color',colorArray{2,i},'LineWidth',1.5);
    hold off
    c = c + 1;
    
end

figure(fig2);
hold all




% color the regions

patch('vertices',v1,'faces',f1,'FaceColor',spec.LIP{2},'FaceAlpha',0.15,'edgeColor','none')

v2 = [[t,fliplr(t)]',[LIP_UB,fliplr(FW_UB)]'];
f2 = f1;

v3 = [[t,fliplr(t)]',[LIP_LB,fliplr(FW_LB)]'];
f3 = f2;

patch('vertices',v2,'faces',f2,'FaceColor',spec.LIPPFW{2},'FaceAlpha',0.15,'edgeColor','none')
patch('vertices',v3,'faces',f3,'FaceColor',spec.LIPPFW{2},'FaceAlpha',0.15,'edgeColor','none')

x1 = [t(1) t(end) t(end)];
y1 = [FW_UB(1) FW_UB(end) 1];

x2 = [t(1) t(1) t(end)];
y2 = [-1 FW_LB(1) FW_LB(end)];

patch(x1,y1,'r','FaceAlpha',0.15,'edgeColor','none')
patch(x2,y2,'r','FaceAlpha',0.15,'edgeColor','none')

xline(0,':');
yline(0,':');
hold off

xlabel('$x$ [m]',ltxFMT{:})
ylabel('$\dot{x}$[m/s]',ltxFMT{:})
legend('Initial state','Final state','$x+\frac{\dot{x}}{\omega} = 0$','$|x+\frac{\dot{x}}{\omega}| = \delta$','$|x+\frac{\dot{x}}{\omega}| = \delta + \frac{\tau}{mg}(e^{- \omega T}-1)^2$','Location','northwest',ltxFMT{:})



figure(fig3)
subplot(3,1,1)
ylabel('$\theta$ [rad]',ltxFMT{:})
yline(0,':');
yline(thetaMax,':','$\theta_{max}$','LabelHorizontalAlignment','center','interpreter','latex',ltxFMT{:});
ylim([0, .8])
subplot(3,1,2)
ylabel('$\dot{\theta}$ [rad/s]',ltxFMT{:})
subplot(3,1,3)
ylabel('$\tau$ [Nm]',ltxFMT{:})
xlabel('time [s]',ltxFMT{:})
ylim([-150,150])

exportgraphics(fig2,'figures/LIPPFW_stability.png','Resolution',300)
exportgraphics(fig3,'figures/LIPPFW_Maximum.png','Resolution',300)


%% Functions


function [t_FW, X_FW, U_bb] = simFW(tau,X0,modelParam,cp_max)
%   simFW simulates the LIP+FW dynamics given a maximum FW torque: tau [NM}
%   initial condition X0 = [x0,xdot0]^T, model parameters, and the maximum
%   recoverable capture point: cp_max such that the FW control action will
%   be proportional to initial capture point


% extract model parameters
g   = modelParam.gravity;
z0  = modelParam.restHeight;
m   = modelParam.mass;
J   = modelParam.inertia;
Kp  = modelParam.pGain;
Kd  = modelParam.dGain;
delta = modelParam.supportSize;
thetaMax = modelParam.theta;

omega = sqrt(g/z0);

% maximum time for the bang bang controler
Tmax = sqrt(J*thetaMax/tau);

% calculate the initial capture point
cp = X0(1) + X0(3)/omega;

% Set the bang-bang control time to be proportional to the maximum
% recoverable capture point
Tbb = min(1,abs(cp/cp_max))*Tmax;

% Simulate the FW dynamics
[t_bb1, X_bb1] = ode45( @(t,X) FlyWheelLIP(t,X,tau,modelParam,delta) , [0,Tbb], X0);
[t_bb2, X_bb2] = ode45( @(t,X) FlyWheelLIP(t,X,-tau,modelParam,delta) , [Tbb,2*Tbb], X_bb1(end,:));

% store data
t_FW = [t_bb1; t_bb2(2:end)];
X_FW = [X_bb1; X_bb2(2:end,:)];
U_bb = [0,tau; Tbb, tau; Tbb, -tau; 2*Tbb, -tau];
end


function Dy = LIP(t,y,modelParam)
% LIP is function used in to solve the second order
% dynamics of the LIP model with PD ankle control. The state variable is
% y = [x; x']. The function takes (t: time , y: state,
% modelParam: a structure containing the model parameters

% extract model Parameters
g   = modelParam.gravity;
z0  = modelParam.restHeight;

Kp  = modelParam.pGain;
Kd  = modelParam.dGain;

delta = modelParam.supportSize;


% Determine Desired CoP location
xCop = Kp.*y(1) + Kd.*y(2);

xCop = max(-delta, min(xCop,delta));

Dy = zeros(2,1); % initialize output

Dy(1) = y(2);
Dy(2) = g/z0*(y(1) - xCop);


end

function Dy = FlyWheelLIP(t,y,tau,modelParam,CoP)
% FlyWHeelLIP is function used in to solve the second order
% system consisting of the LIP model with a flywheel. The state variable is
% y = [x; theta; x'; theta']. The function takes (t: time , y: state,
% tau: torque, modelParam: a structure containing the model parameters,
% pushParam: a structure that contains the impulsive push parameters)

% Set model Parameters
g   = modelParam.gravity;
z0  = modelParam.restHeight;
m   = modelParam.mass;
J   = modelParam.inertia;

Kp  = modelParam.pGain;
Kd  = modelParam.dGain;

delta = modelParam.supportSize;


% Determine Desired CoP location
if(nargin < 5)
    xCop = Kp.*y(1) + Kd.*y(3);
    
    xCop = max(-delta, min(xCop,delta));
else
    xCoP = CoP;
end

Dy = zeros(4,1);    % initalize output

Dy(1) = y(3);                                                   % x'
Dy(2) = y(4);                                                   % theta'
Dy(3) = g/z0*(y(1)-delta) - tau/(m*z0);                          % x''
Dy(4) = tau/J;                                                  % theta''

end

function addArrow(x,xdot,range,nArr,color)

w = 0.85;
h = 0.6;
gcf();

N = length(x);




for i = floor(range(1)*N):floor((floor((range(2)-range(1))*N))/nArr):floor(range(2)*N)
    
    p1 = [x(i) xdot(i)];
    p2 = [x(i+1) xdot(i+1)]; 
    
    
    
    if (abs(x(i) + xdot(i)/sqrt(9.81/0.9))> 0.4) || (p1(1) == p2(1) && p1(2) == p2(2))
        return
    end
    
    
    arrow3(p1,p2,color,w,h)
    
end
end


