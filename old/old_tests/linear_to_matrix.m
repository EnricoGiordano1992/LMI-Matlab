mu = 0.01;
Ms = 0.61;
Bs = 11;
Kv = 87.8;
Kp = 4000;

s = tf('s');

num = mu * Kp;
den = Ms*s^2 + (Bs + Kv)*s + Kp;

H = num / den

discrete_sys = c2d(H, 0.001);

disp('Ys(s) Model:')

ss(discrete_sys)

%{

  a = 
            x1       x2
   x1    1.844  -0.8505
   x2        1        0
 
  b = 
             u1
   x1  0.007812
   x2         0
 
  c = 
             x1        x2
   y1  0.003977  0.003768
 
  d = 
       u1
   y1   0
 
Sample time: 0.001 seconds

%}


Mm = 0.64;
Bm = 0.64;

s = tf('s');

num = 1;
den = (0.64*s + 0.64);

H = num / den;

discrete_sys = c2d(H, 0.001);

disp('Ym(s) Model:')

[A, B, C, D] = ssdata(ss(H))


%{

  a = 
          x1
   x1  0.999
 
  b = 
            u1
   x1  0.03125
 
  c = 
            x1
   y1  0.04998
 
  d = 
       u1
   y1   0
 
Sample time: 0.001 seconds


%}

