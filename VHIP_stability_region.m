%% Numerically determination of the region of stability of 

clc
close all
clear 
%%
addpath("C:\Users\tomhu\OneDrive - Queen's University\Documents\fileExchange\arrow3")
ltxFMT = {'Interpreter','latex'};

spec = plotSpec();

modelParameters

omega = sqrt(g/z0);

% Store variables in structure
modelParam.gravity          = g;
modelParam.restHeight       = z0;
modelParam.mass             = m;

T  = sqrt(J*thetaMax/tauMax);


%% point foot capture region

dz_min = z0 - zmin;
dz_max = zmax - z0;

xDot = linspace(-0.6,0.6,100);

xCP_LIP = 1/omega*xDot;

xCP_LIPPFWub = xCP_LIP + tauMax/(m*g)*(1-exp(-T*omega))^2;
xCP_LIPPFWlb = xCP_LIP - tauMax/(m*g)*(1-exp(-T*omega))^2;

xCP_VHIPzmin = (sqrt(2*dz_min/g)+ (zmin)/(sqrt(g*zmin)+sqrt(2*g*dz_min)))*xDot;
xCP_VHIPzmax = (z0*(sqrt(2*dz_max)+sqrt(zmax)))/(sqrt(g)*(z0 + 2*dz_max + sqrt(2*zmax*dz_max)))*xDot;

fig1 = figure('Position',[100 100 750 300]);


v1 = [[xCP_VHIPzmin';flip(xCP_VHIPzmax')],[xDot';flip(xDot')]];
f1 = 1:1:2*length(xDot);

v2 = [ [xCP_LIPPFWub'; flip([xCP_VHIPzmax(1:50),xCP_VHIPzmin(51:end)]')],[xDot';flip(xDot')]];
f2 = 1:1:2*length(xDot);

v3 = [ [xCP_LIPPFWlb'; flip([xCP_VHIPzmin(1:50),xCP_VHIPzmax(51:end)]')],[xDot';flip(xDot')]];
f3 = 1:1:2*length(xDot);

patch('vertices',v2,'faces',f2,'FaceColor',spec.hip{3},'FaceAlpha',0.15)
hold on
patch('vertices',v3,'faces',f3,'FaceColor',spec.hip{3},'FaceAlpha',0.15,'handlevisibility','off')
hold on
patch('vertices',v1,'faces',f1,'FaceColor',spec.toe{3},'FaceAlpha',0.15)
hold on
plot(xCP_LIP,xDot,'--','lineWidth',1.5)


xlim([-0.2,0.2])
ylim([-0.5,0.5])

xlabel('$\xi$ [m]',spec.ltxFMT{:});
ylabel('$\dot{x}_0$ [m/s]',spec.ltxFMT{:});

legend('$\xi_{\textrm{LIPPFW}}$','$\xi_{\textrm{VHIP}}$','$\xi_{\textrm{LIP}}$',spec.ltxFMT{:},'location','northwest')

exportgraphics(fig1,'figures\capturePositionInitVel.png','resolution',300);

%%
X0Array = genTestPoints(omega,delta);

colorArray = hot(25);
%% vertial acceleration
t = linspace(-1,1);

zBounds = [zmin zmax];

DDz_array = [ 1, 3 , 5 , 9.81];

label = '$\ddot{z}_{bb}$ =';
fig2 = figure('Position',[100 100 750 300]);
hold all
c = 1;
for DDz = DDz_array 
    
    [stableCon, unstableCon] = stabilityCheck(delta,DDz,zBounds,modelParam);
    
    bnds = boundary(stableCon);
    
    VHIPbnds = stableCon(bnds,:);
    
    if DDz ==3
        StabilityBounds = VHIPbnds;
    end
    
    
    plot(VHIPbnds(:,1),VHIPbnds(:,2),'Color',colorArray(4*c,:),'lineWidth',1.5);
 
    c = c +1;
end

lipUB = (-t + delta)*omega;
lipLB = (-t - delta)*omega;

vhipUB = (-t + delta)*sqrt(g/zmax);
vhipLB = (-t + delta)*sqrt(g/zmax);


flyWheelUB = (-t + delta + tauMax/(m*g)*(exp(-omega*T)*(exp(omega*T) - 1)^2))*omega;
flyWheelLB = (-t - delta - tauMax/(m*g)*(exp(-omega*T)*(exp(omega*T) - 1)^2))*omega;

plot(t,lipUB,'Color',spec.LIP{2},'HandleVisibility','off')
plot(t,lipLB,'Color',spec.LIP{2},'HandleVisibility','off')
plot(t,flyWheelUB,'Color',spec.LIPPFW{2},'HandleVisibility','off')
plot(t,flyWheelLB,'Color',spec.LIPPFW{2},'HandleVisibility','off')

patch('vertices',[[t,fliplr(t)]',[lipUB,fliplr(lipLB)]'],'faces',1:1:2*length(t),'FaceColor',spec.LIP{2},'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','on')
patch('vertices',[[t,fliplr(t)]',[lipUB,fliplr(flyWheelUB)]'],'faces',1:1:2*length(t),'FaceColor',spec.LIPPFW{2},'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','on')
patch('vertices',[[t,fliplr(t)]',[lipLB,fliplr(flyWheelLB)]'],'faces',1:1:2*length(t),'FaceColor',spec.LIPPFW{2},'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off')

x1 = [t(1) t(end) t(end)];
y1 = [flyWheelUB(1) flyWheelUB(end) 1];

x2 = [t(1) t(1) t(end)];
y2 = [-1 flyWheelLB(1) flyWheelLB(end)];

patch(x1,y1,'r','FaceAlpha',0.1,'EdgeColor','none')
patch(x2,y2,'r','FaceAlpha',0.1,'EdgeColor','none')

xline(0,'k:');
yline(0,'k:');

hold off
xlim([-0.2 0.2])
ylim([-1 1])


xlabel('$x_0$ [m]','interpreter','latex')
ylabel('$\dot{x}$ [m/s]','interpreter','latex')
legend( [label , ' ' , num2str(DDz_array(1)), 'm/s$^2$'], [label , ' ' , num2str(DDz_array(2)), 'm/s$^2$'] ,...
        [label , ' ' , num2str(DDz_array(3)), 'm/s$^2$'], [label , ' ' , num2str(DDz_array(4)), 'm/s$^2$'],...
       'LIP capture region','LIPPFW capture region','interpreter','latex','Location','northwest')

 exportgraphics(fig2,'figures/zddot_sensitivity.png','Resolution',300)

%% zBounds variations

zMax_array = [0.96 0.975 1 1.1 ];

label = '$z_{\max}-z_0$ =';

fig3 = figure('Position',[100 100 750 300]);
hold all
c = 1;
for zMax = zMax_array 

    zBounds(2) = zMax;
    
    [stableCon, unstableCon] = stabilityCheck(delta,DDz_max,zBounds,modelParam);
    
    bnds = boundary(stableCon);
    
    VHIPbnds = stableCon(bnds,:);
    
    
    plot(VHIPbnds(:,1),VHIPbnds(:,2),'Color',colorArray(4*c,:),'lineWidth',1.5);
    c = c+1;
end


plot(t,lipUB,'Color',spec.LIP{2},'HandleVisibility','off')
plot(t,lipLB,'Color',spec.LIP{2},'HandleVisibility','off')
plot(t,flyWheelUB,'Color',spec.LIPPFW{2},'HandleVisibility','off')
plot(t,flyWheelLB,'Color',spec.LIPPFW{2},'HandleVisibility','off')

patch('vertices',[[t,fliplr(t)]',[lipUB,fliplr(lipLB)]'],'faces',1:1:2*length(t),'FaceColor',spec.LIP{2},'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','on')
patch('vertices',[[t,fliplr(t)]',[lipUB,fliplr(flyWheelUB)]'],'faces',1:1:2*length(t),'FaceColor',spec.LIPPFW{2},'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','on')
patch('vertices',[[t,fliplr(t)]',[lipLB,fliplr(flyWheelLB)]'],'faces',1:1:2*length(t),'FaceColor',spec.LIPPFW{2},'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off')

x1 = [t(1) t(end) t(end)];
y1 = [flyWheelUB(1) flyWheelUB(end) 1];

x2 = [t(1) t(1) t(end)];
y2 = [-1 flyWheelLB(1) flyWheelLB(end)];

patch(x1,y1,'r','FaceAlpha',0.1,'EdgeColor','none')
patch(x2,y2,'r','FaceAlpha',0.1,'EdgeColor','none')

xline(0,'k:');
yline(0,'k:');

hold off
xlim([-0.20 0.20])
ylim([-1 1])

%title('VHIP Sensitivity $z_{max}$ variation','interpreter','latex')
xlabel('$x_0$ [m]','interpreter','latex')
ylabel('$\dot{x}$ [m/s]','interpreter','latex')
legend([label , ' ' , num2str(zMax_array(1)-z0), 'm'], [label , ' ' , num2str(zMax_array(2)-z0), 'm'] ,...
       [label , ' ' , num2str(zMax_array(3)-z0), 'm'], [label , ' ' , num2str(zMax_array(4)-z0), 'm'],...
       'LIP capture region','LIPPFW capture region','interpreter','latex','Location','northwest')
   
exportgraphics(fig3,'figures/VHIP_zmaxSensitivity.png');
    
%% Simulation
x0 = 0;
xDotArray = [delta*omega/4*3, delta*omega, (delta*omega*0.2 + 0.6639*0.8) , 0.67];

fig4 = figure('Position',[100 100 750 300]);
xlim([-0.02 0.16])
ylim([delta*omega/6 0.7])
daspect([0.08,0.6,1])

fig5 = figure('Position',[100 100 750 450]);

colorArray = {'_a','^a','p','f';
              '#2E6815', '#93BB6D', '#BFC01B', '#B32121'};
c = 1;    
for xdot0 = xDotArray
    
    X0 = [x0,xdot0];
    
    [T,Q,U] = simVHIP(X0,modelParam,0.6639);
    
    figure(fig4)
    hold all
    % plot the trajectory
    if c == 1
       
        plot(Q(:,1),Q(:,3),'Color',colorArray{2,c},'HandleVisibility','off','LineWidth',1.5);
         plot(Q(1,1),Q(1,3),'ko','MarkerSize',4);
        plot(Q(end,1),Q(end,3),'kx')
        
        % plot the different regions of stability
        plot(t,lipUB,'Color',spec.LIP{2},'HandleVisibility','on','LineWidth',2)
        plot(t,vhipUB,'--','Color',spec.LIP{2},'HandleVisibility','on','LineWidth',2)
        plot(StabilityBounds(:,1),StabilityBounds(:,2),'Color',spec.VHIP{2},'lineWidth',2);
        %plot(t,flyWheelUB,'Color','#0072BD','HandleVisibility','on','LineWidth',1.5)
        %plot(t,lipLB,'Color','#0072BD','HandleVisibility','off','LineWidth',1.5)
        %plot(t,flyWheelLB,'Color','#77AC30','HandleVisibility','off','LineWidth',1.5)
        
        
        addArrow(Q(:,1),Q(:,3),[0.3,0.8],1.9,colorArray{1,c});
        
        
    else
        
        plot(Q(:,1),Q(:,3),'Color',colorArray{2,c},'HandleVisibility','off','LineWidth',1.5);
        plot(Q(1,1),Q(1,3),'ko','MarkerSize',4,'HandleVisibility','off');
        plot(Q(end,1),Q(end,3),'kx','HandleVisibility','off');
        
        addArrow(Q(:,1),Q(:,3),[0.3,0.8],1.9,colorArray{1,c});
    end
    
    hold off
    
    figure(fig5);
    subplot(3,1,1)
    hold on
    plot(T,Q(:,2),'Color',colorArray{2,c},'LineWidth',1.5)
    subplot(3,1,2)
    hold on
    plot(T,Q(:,4),'Color',colorArray{2,c},'LineWidth',1.5)
    subplot(3,1,3)
    hold on
    plot(U(:,1),U(:,2),'Color',colorArray{2,c},'LineWidth',1.5);
    hold off
    c = c + 1;
    
end

figure(fig4);
hold all

patch('vertices',[[t,fliplr(t)]',[lipUB,fliplr(lipLB)]'],'faces',1:1:2*length(t),'FaceColor',spec.LIP{2},'FaceAlpha',0.15,'EdgeColor','none')
patch('Faces',[1:1:length([StabilityBounds(236:303,1)',t])],'Vertices',[[StabilityBounds(236:303,1)',t]',[StabilityBounds(236:303,2)',lipUB]'],'FaceColor',spec.VHIP{2},'FaceAlpha',0.15,'EdgeColor','none')
patch([StabilityBounds(236:303,1)',-1,1],[StabilityBounds(236:303,2)',10,10],'r','FaceAlpha',0.15,'EdgeColor','none')
xline(0,':');
yline(0,':');


xlabel('$x$ [m]',ltxFMT{:})
ylabel('$\dot{x}$[m/s]',ltxFMT{:})
legend('Initial state','Final State','$|x+\frac{\dot{x}}{\omega(z_0)}| = \delta$','$|x+\frac{\dot{x}}{\omega(z_{\max})}| = \delta$',...
       'VHIP Stability bound','Location','northeast',ltxFMT{:})



figure(fig5)
subplot(3,1,1)
ylabel('$z$ [m]',ltxFMT{:})
yline(zmax,':','$z_{\max}$','LabelHorizontalAlignment','center',ltxFMT{:});
yline(z0,':','$z_{0}$','LabelHorizontalAlignment','right',ltxFMT{:});
ylim([0.9, 1.05])
subplot(3,1,2)
ylabel('$\dot{z}$ [m/s]',ltxFMT{:})
subplot(3,1,3)
ylabel('$\ddot{z}$ [m/s$^2$]',ltxFMT{:})
xlabel('time [s]',ltxFMT{:})
ylim([-4,4])

exportgraphics(fig4,'figures/VHIP_stability.png','Resolution',300)
exportgraphics(fig5,'figures/VHIP_Maximum.png','Resolution',300)

% %% Support Area Variation
% 
% zbounds(2) = 1.1;
% 
% delta_array = [ 0.1 0.15 0.2 0.3 ];
% 
% label = '$\delta_{sup}^\pm$ =';
% 
% figure(3)
% hold all
% for sup = delta_array 
%     
%     [stableCon, unstableCon] = stabilityCheck(sup,DDz_max,zBounds,modelParam);
%     
%     bnds = boundary(stableCon);
%     
%     VHIPbnds = stableCon(bnds,:);
%     
%     plot(VHIPbnds(:,1),VHIPbnds(:,2));
% 
% end
% 
% capturepointUBs = zeros(size(delta_array,2),size(t,2));
% 
% capturepointLBs = capturepointUBs;
% 
% 
%     
% for i = 1:size(delta_array,2)
%     
%     capturepointUBs(i,:) = (-t + delta_array(i))*omega;
%     capturepointLBs(i,:) = (-t - delta_array(i))*omega;
%     
%     if i == 1
%     plot(t,capturepointUBs(i,:),'--r')
%     plot(t,capturepointLBs(i,:),'--r','HandleVisibility','off')
%     
%     else
%     plot(t,capturepointUBs(i,:),'--r','HandleVisibility','off')
%     plot(t,capturepointLBs(i,:),'--r','HandleVisibility','off')    
%     end
%     
% end
% 
% hold off
%                
% xlim([-0.4 0.4])
% ylim([-1 1])
% title('VHIP Sensitivity with support size variation','interpreter','latex')
% xlabel('Initial CoM position [$m$]','interpreter','latex')
% ylabel('Initial CoM Velocity [$\frac{m}{s}$]','interpreter','latex')
% legend([label, ' ' , num2str(delta_array(1)), '$m$'],[label, ' ' , num2str(delta_array(2)), '$m$'],...
%        [label, ' ' , num2str(delta_array(3)), '$m$'],[label, ' ' , num2str(delta_array(4)), '$m$'],...
%        'LIP stability bounds','interpreter','latex')
% 
% 
% %% 
% supDataTab = readtable('VHIP_sup.csv');
% supData = table2array(supDataTab);
% 
% 
% figure(4) 
% plot(supData(:,1),supData(:,4),'kx')
% xlabel('support size (m)')
% ylabel('Initial Capture Point $\xi = x_0 + \frac{\dot{x}_0}{\omega}$','interpreter','latex')
% title('Capture point which VHIP falls to LIP stability region')
% 
% figure(5)
% plot(supData(:,1),supData(:,5),'kx');
% xlabel('support size (m)')
% ylabel('capture point as percentage of the support size $\frac{\xi}{\delta}$','interpreter','latex')
% title('Capure point as percentage of supportsize Where VHIP falls to LIP stability')
% ylim([0,1])


%% functions

function [t,Q,U] = simVHIP(X0,modelParam,maxCP)
    
    z0 = modelParam.restHeight;
    zmax = modelParam.zbound(2);
    delta = modelParam.delta;
    DDz   = modelParam.DDz;
    
    Tmax = sqrt(2*(zmax-z0)/(2*DDz));
    
    Tbb = min(Tmax,X0(2)/maxCP*Tmax);
    
    Q0 = [X0(1),z0,X0(2),0];
    
    
    [t1,Q1] = ode45(@(t,q) VHIP(t,q,delta,DDz),[0 Tbb],Q0);
    [t2,Q2] = ode45(@(t,q) VHIP(t,q,delta,-DDz),[Tbb,2*Tbb],Q1(end,:));
      
    Q = [Q1;Q2(2:end,:)];
    U = [0,DDz;Tbb,DDz;Tbb,-DDz;2*Tbb,-DDz];
    t = [t1;t2(2:end)];

end







function [stableCon, unstableCon] = stabilityCheck(delta,DDz_max,zBounds,modelParam)





g   = modelParam.gravity;
z0  = modelParam.restHeight;

omega = sqrt(g/z0);

X0Array = genTestPoints(omega,delta);


i_max = size(X0Array,1);

stable_count    = 1;
unstable_count  = 1;

stableCon   = zeros(i_max,2);
unstableCon = zeros(i_max,2);



for i = 1:size(X0Array,1)
    
        
        x0    = X0Array(i,1);
        xdot0 = X0Array(i,2);
        
        Q0 = [x0;z0];
        Qdot0 = [xdot0;0];
        
        orbit = 0.5*Qdot0(1)^2 - g/(2*z0)*Q0(1)^2;

        if orbit >= 0
            t1 = sqrt(2*(zBounds(2)-z0)/(2*DDz_max));
            t2 = 2*t1;
            
            q0 = [Q0;Qdot0];
            [~,q1] = ode45(@(t,q) VHIP(t,q,delta,DDz_max),[0 t1],q0);
            [~,q2] = ode45(@(t,q) VHIP(t,q,delta,-DDz_max),[t1 t2],q1(end,:));
            
            omega = sqrt(g/zBounds(2));

        else
            t1 = sqrt(2*(z0-zBounds(1))/(2*DDz_max));
            t2 = 2*t1;
            
            q0 = [Q0;Qdot0];
            [~,q1] = ode45(@(t,q) VHIP(t,q,delta,-DDz_max),[0 t1],q0);
            [~,q2] = ode45(@(t,q) VHIP(t,q,delta,DDz_max),[t1 t2],q1(end,:));
            
             omega = sqrt(g/zBounds(1));

        end


         xf = q2(end,1);
         xdotf = q2(end,3);
         
         cp = xf + xdotf/omega;
         
         if abs(cp) > delta
             
             unstableCon(unstable_count,:) =  [x0 xdot0];
             
             unstable_count = unstable_count + 1;
         else
             stableCon(stable_count,:) =  [x0 xdot0];
             
             stable_count = stable_count + 1;
         end
            
        
        
   
end

stableCon = stableCon(1:stable_count,:);
unstableCon = unstableCon(1:unstable_count,:);

end
    
function Xarray = genTestPoints(omega,delta)


N = 100;
n = 100;


x = linspace(-2*delta,2*delta,N);

xdotUB = (delta - x)*omega;
xdotLB = (-delta - x)*omega;

p = linspace(-0.3,0.3,n);

c = 1;
Xarray = zeros(N*n*2,2);


    for i = 1:length(x)
        for j= 1:length(p) 
            
        Xarray(c,:) = [x(i),xdotUB(i)+p(j)];
        c = c + 1;
        Xarray(c,:) = [x(i),xdotLB(i)+p(j)];
        c = c + 1;   
        end
    end

end

function [Dy,xCop] = VHIP(t,y,delta,zc)

g = 9.81;

Kp = 3;
Kd = 1;

xCop = Kp.*y(1) + Kd.*y(3);

xCop = max(-delta, min(xCop,delta));

Dy = zeros(4,1);

Dy(1) = y(3);
Dy(2) = y(4);
Dy(3) = (y(1)-xCop)./y(2).*(g+zc);
Dy(4) = zc;

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