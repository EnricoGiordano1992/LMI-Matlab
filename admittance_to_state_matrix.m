mu = 0.01;
Kf = 1;
Ms = 0.61;
Bs = 11;
Kv = 87.8;
Kp = 4000;

s = tf('s');

num = mu * Kp;
den = Ms*s^2 + (Bs + Kv)*s + Kp;

Ys = num / den


Mm = 0.64;
Bm = 0.64;

s = tf('s');

num = 1;
den = (0.64*s + 0.64);

Ym = num / den

s = tf('s');

Y = [Ym -Kf*Ym; -Ym*Ys ((Mm*s^2+Bm*s+mu*Kf*Kp)/(Ms*s^2+(Bs+Kv)*s+Kp))*Ym]

[A, B, C, D] = ssdata(ss(Y))
