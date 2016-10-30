Ms = 0.64;
Bs = 12;
Mm = 0.64;
Bm = 12;

Ke = 0%0.01;
Me = 0%1;
De = 0%0.01; %damping
Xe = 9;

K = [
  -20.2018   22.4536  -19.0278   20.1187
  -30.5351   55.1758  -30.9200   52.9287
]


A_master = [0 1; 0 -Bm/Mm];
B_master = [0; 1/Mm];
C_master = [1 0; 0 1];
D_master = [0; 0];

A_slave = [0 1; 0 -Bs/Ms];
B_slave = [0; 1/Ms];
C_slave = [1 0; 0 1];
D_slave = [0; 0];