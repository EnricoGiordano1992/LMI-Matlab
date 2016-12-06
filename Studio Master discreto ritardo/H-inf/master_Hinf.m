Mm = 0.64;
Bm = 12;

A = [0 1; 0 -Bm/Mm];
B = [0; 1/Mm];
C = [1 0; 0 1];
D = [0; 0];
B0 = [0; 1/Mm];
D0 = [0];

Ts = 0.001;
mu = 1;

%z = tf('z', Ts);

csys = ss(A,B,C,D);
%dsys = z^(-4)*c2d(csys, Ts);

dsys = c2d(csys, Ts);

n_delay = 20;
phi = dsys.a;
gamma1 = dsys.b;
gamma2 = zeros(size(dsys.b,1), size(dsys.b,2));

Az1 = [phi gamma1 gamma2 zeros(size(dsys.b,1), n_delay-1)];
Az21 = [zeros(n_delay+1,size(dsys.a,2)+size(dsys.b,2))];
Az22 = [eye(n_delay+1,n_delay)];

A_z = [Az1; Az21 Az22];

B_z=[zeros(size(A_z,1)-1,1); 1]
C_z=[dsys.c zeros(size(dsys.c,1),n_delay+1)]


A = dsys.a;
B = dsys.b;
C = dsys.c;
D = dsys.d;

Ad = A_z;
Bd = B_z;
Cd = C_z;
Dd = dsys.d;

n = size(Ad,1);
p = size(Bd,1);
q = size(Bd,2);
r = size(Dd,1);
s = size(Dd,2);

C0 = [20 1 zeros(1,n_delay+1)];


P = sdpvar(n,n, 'symmetric');
L = sdpvar(q,n);
%X = sdpvar(n,n);

b11 = P;
%b12 = Ad*X+Bd*L;
b12 = Ad*P+Bd*L;
b13 = Bd;
b14 = zeros(n,q);
b21 = b12';
%b22 = X+X'-P;
b22 = P;
b23 = zeros(n, q);
%b24 = X'*Cd'+L'*Dd';
%b24 = X'*C0'+L'*D0';
b24 = P*C0'+L'*D0';
b31 = b13';
b32 = b23';
b33 = eye(q, q);
b34 = D0';
b41 = b14';
b42 = b24';
b43 = b34';
b44 = mu*eye(q,q);

lmi1 = [b11 b12 b13 b14];
lmi2 = [b21 b22 b23 b24];
lmi3 = [b31 b32 b33 b34];
lmi4 = [b41 b42 b43 b44];

LMI = [lmi1; lmi2; lmi3; lmi4];

F = [P >= 0, LMI >=0];

diagnostics = solvesdp(F);
disp(diagnostics.problem);
if diagnostics.problem == 0
 disp('Feasible')
 L_s = value(L);
 P_s = value(P);
 
 K = L_s * inv(P_s)
 K = K(1:2)

elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end

