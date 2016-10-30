
clear all

addpath /home/riccardo/Matlab/MIQP
% function [xmin, fmin, flag, Extendedflag] = ...
%           miqp(H, f, A, b, Aeq, beq, vartype, lb, ub, x0, Options)




options = [];
options.integtol = 1e-6;
options.solver   = 'quadprog';

% [xsol,fsol,flag,ef] = miqp(Q, b, C, d, [], [], ivar, vlb, vub, x0, options);

% ===================================================
%   goal
% ===================================================

TsHuman = 1e-4;
T1 = [0:TsHuman:10];
k1 = linspace(1,2,length(T1));

f1 = 0.225; % [Hz]
Y1 = k1.*sin(2*pi*f1*T1);

Y2 = Y1(end)*ones(1,length(T1));

f3 = 0.325;
T3 = [0:TsHuman:10];
k3 = linspace(2,0.5,length(T3));
Y3 = k3.*sin(2*pi*f3*T3+pi/2);

Y = [Y1 Y2 Y3];
T = [0:length(Y)-1]*TsHuman;

humanTrajectory.time = T';
humanTrajectory.signals.values = Y';
humanTrajectory.signals.dimensions = 1;

% figure
% plot(humanTrajectory.time,humanTrajectory.signals.values)

xE = 1.8;
fH = 20;

KPh_pos = 30;
KIh_pos = 20;
KDh_pos = 5;

KPh_force = 0.1;
KIh_force = 0;
KDh_force = 0;

kP_slave = 10;% 20;
kD_slave = 5;% 5;
kI_slave = 1;% 5;

% ===================================================
%   plant data
% ===================================================

Ts = 1e-3;

% master + operator (SISO case)
Mm = 1;% 3.5;
Mh = 0.35; %2;
Dm = 1;%2.2;
Dh = 5; %1;
Km = 0;
Kh = 5; %60;
Mmh = Mm + Mh;
Dmh = Dm + Dh; 
Kmh = Km + Kh;
ssAmh = zeros(2,2);
ssAmh(1,2) = 1;
ssAmh(2,1) = -inv(Mmh)*Kmh;
ssAmh(2,2) = -inv(Mmh)*Dmh;
ssBmh = zeros(2,1);
ssBmh(2,1) = inv(Mmh);
ssCmh = eye(2,2);
ssDmh = zeros(2,1);

% slave + environment (SISO case)
Ms = 1; %3.5;
Me = 0;
Ds = 1; %1.6;
De = 10;%10; %10
Ks = 0;
Ke = 50;%1000;

ssAs = zeros(2,2);
ssAs(1,2) = 1;
ssAs(2,1) = -inv(Ms)*Ks;
ssAs(2,2) = -inv(Ms)*Ds;
ssBs = zeros(2,1);
ssBs(2,1) = inv(Ms);
ssCs = eye(2,2);
ssDs = zeros(2,1);


s = tf('s');
Es = 1/(Me*s^2+De*s+Ke);
[ssAe,ssBe,ssCe,ssDe] = ssdata(Es);

% NCS_MPC_QoS_plant(Ms,Ds,Ks,Me,De,Ke);



duration_time = T(end);

return

% sim('NCS_MPC_QoS')
% sim('NCS_MPC_QoS_v2')
% sim('NCS_MPC_QoS_v3')
sim('teleoperation_MPC_QoS')


