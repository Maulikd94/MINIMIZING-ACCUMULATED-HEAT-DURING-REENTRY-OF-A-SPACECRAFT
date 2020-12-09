%%                      3 States
clc;
clear all;
close all;
tic;

%% Global Constants
global R B c1 c2 c3 Sm g0 rho FT
R = 209;                % Radius of Earth in 100,000 ft
B = 4.26;               % Beta
c1 = 1.174;             % Aerodynamic constant
c2 = 0.9;               % Aerodynamic constant
c3 = 0.6;               % Aerodynamic constant
Sm = 53200;             % Surface to mass ratio of spacecraft
g0 = 3.2172e-4;
rho = 2.704e-3;

%% Teminal Conditions

% States
% Initial Conditions
V0 = 0.36;              % Velocity
G00 = -8.1*pi/180 ;     % Flight path angle
Z0 = 4/R;             % Normalized height
% Final Conditions
Vf = 0.27;
Gf = 0;
Zf = 2.5/R;

% Guesses of State & Co-state at Initial time
LV0 = -0.5;             % Costate for veolcity
LG0 = -0.5 ;               % Costate for flight path angle
LZ0 = -2.2;                 % Costate for normalized height

% Guesses of State & Co-state at Final time
LVf = -4;               % Costate for veolcity
LGf =  3.6;               % Costate for flight path angle
LZf = -30;               % Costate for normalized height

% Guess for Final Time
FT = 225;               % Obtained from text book

I0 = [V0 G00 Z0 LV0 LG0 LZ0]'; %Initial condition for States & CoStates
If = [Vf Gf Zf LVf LGf LZf]'; %Initial condition for States & Costates

%% User Parameters
% Time span for Integrations
tvec1 = linspace(0,0.5*FT,51);
tvec2 = linspace(FT,0.5*FT,51);

MaxNiter = 10;
Niter = 1;
K = 1e-5;    % Update factor
A = Jacob;

%% Main Loop

while Niter<=MaxNiter
    Option1 = odeset('RelTol', 1e-7, 'AbsTol',1e-9);
    %   Integration in Forward in initial time to mid-point time
    [t1,state1] = ode45(@(time,state)State_COstate(time,state),tvec1,I0,Option1);
    
    vt1 = state1(end,1);gt1 = state1(end,2);zt1 = state1(end,3);
    lmvt1 = state1(end,4);lmgt1 = state1(end,5);lmzt1 = state1(end,6);
    
    
    %   Integration in Backward in mid-point time to final time
    [t2,state2] = ode45(@(time,state)State_COstate(time,state),tvec2,If,Option1);
    
    vt2 = state2(end,1);gt2 = state2(end,2);zt2 = state2(end,3);
    lmvt2 = state2(end,4);lmgt2 = state2(end,5);lmzt2 = state2(end,6);
    
    % Error in the Variables
    timediff=t2(end)-t1(end);
    err=state1(end,:)-state2(end,:); 
    IdenMat= eye(6);
    
    % Plot the Variables
    for i= 1:51
        u11(i) = atan((6*state1(i,5))/(9*state1(i,1)*state1(i,4)));
        u12(i) = atan((6*state2(i,5))/(9*state2(i,1)*state2(i,4)));
    end
    u1 = u11*1;
    u2 = u12*1;
    du = u1(end)-u2(end);
    
    % Termination Criteria
    if timediff~=0
        disp('Integration Times are not matching up')
        disp(timediff)
        Savedata(I2)
        break;
    end
    if max(isnan(state1))
        disp('Forward Integration has NaN')
        Savedata(I2)
        break;
    end
    if max(isnan(state2))
        disp('Backward Integration has NaN')
        Savedata(I2)
        break;
    end
    if norm(err) < 1e-5
        toc;
        disp('control convergence to the allowed tolerance');
        Savedata(I2)
        break;
    end
    
    
    %  Calculation of the Jacobian for Forward Integration
    V = I0(1);G = I0(2);Z = I0(3);
    LV = I0(4);LG = I0(5);LZ = I0(6);
    B1=eval(A);
    
    %   Integration in Forward in initial time to mid-point time
    [~,phi1] = ode113(@Transition,tvec1,IdenMat(:),Option1,B1);
    Phi11 = reshape(phi1(end,:),[6,6]);
    
    %  Calculation of the Jacobian for Backward Integration
    V = If(1);G = If(2);Z = If(3);
    LV = If(4);LG = If(5);LZ = If(6);
    B2=eval(A);
    
    %   Integration in Backward in mid-point time to final time
    [~,phi2] = ode113(@Transition,tvec2,IdenMat(:),Option1,B2);
    Phi12 = reshape(phi2(end,:),[6,6]);
    
    NC1 = Correction(err,Phi11,Phi12);
    NC = K.*NC1;
  
    disp(num2str(Niter))
    disp(du)
    disp('Error Norm')
    disp(norm(err(1)))
%     K = K*1.05;
    
    %%%%%%%%%%%%%%%%%%%%%%%% UPDATE SECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correction section
    LV0 = LV0 + NC(1);LG0 = LG0 + NC(2);LZ0 = LZ0 + NC(3);
    LVf = LVf + NC(4);LGf = LGf + NC(5);LZf = LZf + NC(6);
    
    L0 = [LV0 LG0 LZ0]';
    Lf = [LVf LGf LZf]';
    
    % Update in Parameters
    Niter=Niter+1;
    
    % Update in Initial Guesses
        I0 = [V0 G00 Z0 LV0 LG0 LZ0]'; %Initial condition for States & CoStates
    If = [Vf Gf Zf LVf LGf LZf]'; %Initial condition for States & Costates
    
%     plots(t1,state1,t2,state2,Niter,u1,u2)
    
end
figure(1)
plot(t1,state1(:,1));grid
axis([0 FT 0.99*min(min(state1(:,1)),min(state2(:,1))) 1.01*max(max(state1(:,1)),max(state2(:,1)))]);
hold on
plot(t2,state2(:,1))
hold off
xlabel('Time (s)')
ylabel('Velocity')
print('-f1','Velocirty1','-dpng')
 
figure(2)
plot(t1,(180/pi)*state1(:,2)); grid
axis([0 FT 1.1*(180/pi)*min(min(state1(:,2)),min(state2(:,2))) 1.1*(180/pi)*max(max(state1(:,2)),max(state2(:,2)))]);
hold on
plot(t2,(180/pi)*state2(:,2))
hold off
xlabel('Time (s)')
ylabel('Gamma')
 print('-f2','Gamma1','-dpng')

figure(3)
plot(t1,state1(:,3)); grid
axis([0 FT 0.95*min(min(state1(:,3)),min(state2(:,3))) 1.05*max(max(state1(:,3)),max(state2(:,3)))]);
hold on
plot(t2,state2(:,3))
hold off
xlabel('Time (s)')
ylabel('Zeta')
print('-f3','Zeta1','-dpng')

figure(4)
plot(t1,u1); grid
axis([0 FT 1.1*min(min(u1),min(u2)) 1.1*max(max(u1),max(u2))]);
hold on
plot(t2,u2)
hold off
xlabel('Time (s)')
ylabel('Control U')
print('-f4','Control1','-dpng')

figure(5)
plot(t1,state1(:,4));grid
axis([0 FT 0.99*min(min(state1(:,4)),min(state2(:,4))) 1.01*max(max(state1(:,4)),max(state2(:,4)))]);
hold on
plot(t2,state2(:,4))
hold off
xlabel('Time (s)')
ylabel('Lambda Velocity')
print('-f5','Lvelocity1','-dpng')

figure(6)
plot(t1,state1(:,5));grid
axis([0 FT 0.99*min(min(state1(:,5)),min(state2(:,5))) 1.01*max(max(state1(:,5)),max(state2(:,5)))]);
hold on
plot(t2,state2(:,5))
hold off
xlabel('Time (s)')
ylabel('Lambda Gamma')
print('-f6','Lgamma1','-dpng')

figure(7)
plot(t1,state1(:,6));grid
axis([0 FT 0.99*min(min(state1(:,6)),min(state2(:,6))) 1.01*max(max(state1(:,6)),max(state2(:,6)))]);
hold on
plot(t2,state2(:,6))
hold off
xlabel('Time (s)')
ylabel('Lambda Zeta')
print('-f7','Lzeta1','-dpng')

%% PLOT FUNCTION
function plots(t1,state1,t2,state2,Niter,u1,u2)
global R FT
subplot(2,4,1)
plot(t1,state1(:,1));grid
axis([0 FT 0.99*min(min(state1(:,1)),min(state2(:,1))) 1.01*max(max(state1(:,1)),max(state2(:,1)))]);
hold on
plot(t2,state2(:,1))
hold off
xlabel('Time (s)')
ylabel('Velocity')

subplot(2,4,2)
plot(t1,(180/pi)*state1(:,2)); grid
axis([0 FT 1.1*(180/pi)*min(min(state1(:,2)),min(state2(:,2))) 1.1*(180/pi)*max(max(state1(:,2)),max(state2(:,2)))]);
hold on
plot(t2,(180/pi)*state2(:,2))
hold off
xlabel('Time (s)')
ylabel('Gamma')
title(['Iteration No. ',num2str(Niter)])

subplot(2,4,3)
plot(t1,state1(:,3)); grid
axis([0 FT 0.95*min(min(state1(:,3)),min(state2(:,3))) 1.05*max(max(state1(:,3)),max(state2(:,3)))]);
hold on
plot(t2,state2(:,3))
hold off
xlabel('Time (s)')
ylabel('Zeta')

subplot(2,4,4)
plot(t1,u1); grid
axis([0 FT 1.1*min(min(u1),min(u2)) 1.1*max(max(u1),max(u2))]);
hold on
plot(t2,u2)
hold off
xlabel('Time (s)')
ylabel('Control U')

subplot(2,4,5)
plot(t1,state1(:,4));grid
axis([0 FT 0.99*min(min(state1(:,4)),min(state2(:,4))) 1.01*max(max(state1(:,4)),max(state2(:,4)))]);
hold on
plot(t2,state2(:,4))
hold off
xlabel('Time (s)')
ylabel('Lambda Velocity')

subplot(2,4,6)
plot(t1,state1(:,5));grid
axis([0 FT 0.99*min(min(state1(:,5)),min(state2(:,5))) 1.01*max(max(state1(:,5)),max(state2(:,5)))]);
hold on
plot(t2,state2(:,5))
hold off
xlabel('Time (s)')
ylabel('Lambda Gamma')

subplot(2,4,7)
plot(t1,state1(:,6));grid
axis([0 FT 0.99*min(min(state1(:,6)),min(state2(:,6))) 1.01*max(max(state1(:,6)),max(state2(:,6)))]);
hold on
plot(t2,state2(:,6))
hold off
xlabel('Time (s)')
ylabel('Lambda Zeta')

pause(0.1)
end

%% State and Costate Equations

function SC = State_COstate(time,state)
global R B c1 c2 c3 Sm g0 rho FT

V = state(1);
G = state(2);
Z = state(3);
LV = state(4);
LG = state(5);
LZ = state(6);
% FT = state(7);

den = sqrt((c3*LG)^2+(c2*LV*V)^2);     % Obtained from dH/dU = 0
cosu = -c2*LV*V/den;
sinu = -c3*LG/den;
CD = c1 - c2*cosu;
CL = c3*sinu;

Vdot = (-0.5*Sm*rho*exp(-B*R*Z)*V^2*CD -(g0*sin(G))/(1+Z)^2);
Gdot = (0.5*Sm*rho*exp(-B*R*Z)*V*CL + (V*cos(G))/(R*(1+Z)) - g0*cos(G)/(V*(1+Z)^2));
Zdot = (V*sin(G)/R);

% Hamiltonian
% H = 10*V^3*sqrt(rho*exp(-B*R*Z)) + LV*Vdot + LG*Gdot + LZ*Zdot;

% Costate Equations
LVdot = LV*(143.8528*V*exp(-890.34*Z)*((0.81*LV*V)/(0.36*LG^2 + 0.81*LV^2*V^2)^(1/2) + 1.174) + 71.9264*V^2*exp(-890.34*Z)*((0.81*LV)/(0.36*LG^2 + 0.81*LV^2*V^2)^(1/2) - (0.6561*LV^3*V^2)/(0.36*LG^2 + 0.81*LV^2*V^2)^(3/2))) - 0.0047846889952153110047846889952153*LZ*sin(G) - 30.0*V^2*(0.002704*exp(-890.34*Z))^(1/2) - 1.0*LG*(cos(G)/(209.0*Z + 209.0) + (0.00032172000000000000293279289742543*cos(G))/(V^2*(Z + 1.0)^2) - (25.893504*LG*exp(-890.34*Z))/(0.36*LG^2 + 0.81*LV^2*V^2)^(1/2) + (20.97373824*LG*LV^2*V^2*exp(-890.34*Z))/(0.36*LG^2 + 0.81*LV^2*V^2)^(3/2));
LGdot = LG*((V*sin(G))/(209*Z + 209) - (5934686503393837*sin(G))/(18446744073709551616*V*(Z + 1)^2)) + (5934686503393837*LV*cos(G))/(18446744073709551616*(Z + 1)^2) - (LZ*V*cos(G))/209;
LZdot = (12.0373968*V^3*exp(-890.34*Z))/(0.002704*exp(-890.34*Z))^(1/2) - 1.0*LG*((0.00064344000000000000586558579485086*cos(G))/(V*(Z + 1.0)^3) - (209.0*V*cos(G))/(209.0*Z + 209.0)^2 + (23054.02235136*LG*V*exp(-890.34*Z))/(0.36*LG^2 + 0.81*LV^2*V^2)^(1/2)) - 1.0*LV*((0.00064344000000000000586558579485086*sin(G))/(Z + 1.0)^3 + 64038.950976*V^2*exp(-890.34*Z)*((0.81*LV*V)/(0.36*LG^2 + 0.81*LV^2*V^2)^(1/2) + 1.174));
% FTdot = 0;

SC = [Vdot;Gdot;Zdot;LVdot;LGdot;LZdot];%FTdot];
end

%%
function dydot = Transition(time,phiC,B)
phiM = reshape(phiC,[6,6]);
dydotM = B*phiM;
dydot = dydotM(:);
time;

end

%%                                 Solving for the Corrections

function New_SC = Correction(err,Phi11,Phi12)

CLV0 = sym('CLV0','real');
CLG0 = sym('CLG0','real');
CLZ0 = sym('CLZ0','real');
CLVf = sym('CLVf','real');
CLGf = sym('CLGf','real');
CLZf = sym('CLZf','real');
% CFT = sym('CFT','real');

C1 = [0 0 0 CLV0 CLG0 CLZ0]';
C2 = [0 0 0 CLVf CLGf CLZf]';
S1 = Phi11*C1;
S2 = Phi12*C2;
S = S1+S2;
err = reshape(err,[6,1]);

for i = 1:6
    eqn(i) = S(i,:) == err(i);
end

solu = solve([eqn(1),eqn(2),eqn(3),eqn(4),eqn(5),eqn(6)],[CLV0,CLG0,CLZ0,CLVf,CLGf,CLZf]);

CLV0_sol = solu.CLV0;
CLG0_sol = solu.CLG0;
CLZ0_sol = solu.CLZ0;
CLVf_sol = solu.CLVf;
CLGf_sol = solu.CLGf;
CLZf_sol = solu.CLZf;
% CFT_sol = solu.CFT;

% vpa(eqn')
Sol = [CLV0_sol;CLG0_sol;CLZ0_sol;CLVf_sol;CLGf_sol;CLZf_sol];%CFT_sol];
New_SC = double(Sol);
end

%% Jacobian

function JB = Jacob

global R B c1 c2 c3 Sm g0 rho FT
V = sym('V','real'); G = sym('G','real'); Z = sym('Z','real');
LV = sym('LV','real'); LG = sym('LG','real'); LZ = sym('LZ','real');
% FT = sym('FT','real');

% H = 10*V^3*sqrt(rho*exp(-B*R*Z)) + LV*Vdot + LG*Gdot + LZ*Zdot

den = sqrt((c3*LG)^2+(c2*LV*V)^2);     % Obtained from dH/dU = 0
cosu = -c2*LV*V/den;
sinu = -c3*LG/den;

CD = c1 - c2*cosu;
CL = c3*sinu;

% State Equations
Vdot = (-0.5*Sm*rho*exp(-B*R*Z)*V^2*CD -(g0*sin(G))/(1+Z)^2);
Gdot = (0.5*Sm*rho*exp(-B*R*Z)*V*CL + (V*cos(G))/(R*(1+Z)) - g0*cos(G)/(V*(1+Z)^2));
Zdot = (V*sin(G)/R);

% Hamiltonian
H = 10*V^3*sqrt(rho*exp(-B*R*Z)) + LV*Vdot + LG*Gdot + LZ*Zdot;

% Costate Equations
LVdot = -diff(H,V);
LGdot = -diff(H,G);
LZdot = -diff(H,Z);

% Constituent differential equations
JV = vpa(jacobian(Vdot,[V,G,Z,LV,LG,LZ]));
JG = vpa(jacobian(Gdot,[V,G,Z,LV,LG,LZ]));
JZ = vpa(jacobian(Zdot,[V,G,Z,LV,LG,LZ]));

JLV = vpa(jacobian(LVdot,[V,G,Z,LV,LG,LZ]));
JLG = vpa(jacobian(LGdot,[V,G,Z,LV,LG,LZ]));
JLZ = vpa(jacobian(LZdot,[V,G,Z,LV,LG,LZ]));

% dH = vpa(jacobian(H,[V,G,Z,LV,LG,LZ]));

% Jacobian
JB = [JV;JG;JZ;JLV;JLG;JLZ];%dH];

end
