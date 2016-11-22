Mm = 0.64;
Bm = 12;

A = [0 1; 0 -Bm/Mm];
B = [0; 1/Mm];
C = [1 0; 0 1];
D = [0; 0];
C0 = [1 10];
B0 = [0; 1/Mm];
D0 = [0.1];

n = size(A,1);
p = size(B,1);
q = size(B,2);
r = size(D,1);
s = size(D,2);


Ts = 0.001;

mu = 1;

csys = ss(A,B,C,D);
dsys = c2d(csys, Ts);

Ad = dsys.a;
Bd = dsys.b;
Cd = dsys.c;
Dd = dsys.d;

P = sdpvar(n,n, 'symmetric');
M = sdpvar(q,n);
%X = sdpvar(n,n);

b11 = P;
%b12 = Ad*X+Bd*L;
b12 = Ad*P+Bd*M;
b13 = Bd;
b14 = zeros(n,q);
b21 = b12';
%b22 = X+X'-P;
b22 = P;
b23 = zeros(n, q);
%b24 = X'*Cd'+L'*Dd';
b24 = P*C0';
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
 M_s = value(M);
 P_s = value(P);
 
 K = M_s * inv(P_s)

elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end

