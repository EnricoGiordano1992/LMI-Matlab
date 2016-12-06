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

csys = ss(A,B,C,D);

dsys = c2d(csys, Ts);


A = dsys.a;
B = dsys.b;
C = dsys.c;
D = dsys.d;

Cz = C;
Dzu = D;
Bu = B;
Bw = B;

n = size(A,1);
p = size(Bu,1);
q = size(Bu,2);
r = size(Dzu,1);
s = size(Dzu,2);

C0 = [1 10];


P = sdpvar(n,n, 'symmetric');
W = sdpvar(n,n, 'symmetric');
L = sdpvar(q,n);
X = sdpvar(n,n);

LMI1 = [W Cz*X+Dzu*L; W' X+X'-P];

b11 = P;
b12 = A*X+Bu*L;
b13 = Bw;
b21 = b12';
b22 = X+X'-P;
b23 = zeros(n,q);
b31 = b13';
b32 = b23';
b33 = eye(q);

lmi2_1 = [b11 b12 b13];
lmi2_2 = [b21 b22 b23];
lmi2_3 = [b31 b32 b33];

LMI2 = [lmi2_1; lmi2_2; lmi2_3];

F = [trace(W) <= mu, LMI1 >=0, LMI2 >= 0];

diagnostics = solvesdp(F);
disp(diagnostics.problem);
if diagnostics.problem == 0
 disp('Feasible')
 L_s = value(L);
 X_s = value(X);
 
 K = L_s * inv(X_s)

elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end

