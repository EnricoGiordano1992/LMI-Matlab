Ms = 0.61;
Bs = 11;
Kv = 40;
Kp = 40;
Mm = 0.64;
Bm = 12;

K = [-0.5155 -1.4188 -0.5155 0; -0.5155 -0.0003 -0.5155 -1.5375];

A_master = [0 1; 0 -Bm/Mm];
B_master = [0; 1/Mm];
C_master = [0 1];
D_master = 0;

A_slave = [0 1; 0 -Bs/Ms];
B_slave = [0; 1/Ms];
C_slave = [0 1];
D_slave = 0;