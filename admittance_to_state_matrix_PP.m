% mu = 0.01;
% Kf = 1;
% Ms = 0.61;
% Bs = 11;
% Kv = 87.8;
% Kp = 4000;
% 
% s = tf('s');
% 
% num = mu * Kp;
% den = Ms*s^2 + (Bs + Kv)*s + Kp;
% 
% Ys = num / den
% 
% 
% Mm = 0.64;
% Bm = 0.64;
% 
% s = tf('s');
% 
% num = 1;
% den = Mm*s^2 + (Bm + Kv)*s + Kp;
% 
% Ym = num / den
% 
% s = tf('s');
% 
% Y = [Ym -Kf*Ym; -Ym*Ys ((Mm*s^2+Bm*s+mu*Kf*Kp)/(Ms*s^2+(Bs+Kv)*s+Kp))*Ym]
% 
% [A, B, C, D] = ssdata(ss(Y))


mu = 0.01;
Kf = 1;
Ms = 0.61;
Bs = 11;
Kv = 87.8;
Kp = 4000;
Mm = 0.64;
Bm = 0.64;

s = tf('s');

Gm = Mm*s^2 + Bm*s + Kp + Kv*s;
Qm = -Kp -Kv*s;
Gs = Ms*s^2 + Bs*s + Kp + Kv*s;
Qs = -Kp -Kv*s;

%Fs = Xm*Gm + Xs*Qm
%Fm = Xs*Gs + Xm*Qs

det_z = Gm*Gs - Qm*Qs;
Z = [Gs Qm; Qs Gm];
Y = 1/det_z * [Gs -Qm; -Qs Gm];
Ys = s * Y;

[A, B, C, D] = ssdata(ss(Y))
