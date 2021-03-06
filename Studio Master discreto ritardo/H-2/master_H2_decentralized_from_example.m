a1 = [0.8189 0.0863 0.0900 0.0813];
a2 = [0.2524 1.0033 0.0313 0.2004];
a3 = [-0.0545 0.0102 0.7901 -0.2580];
a4 = [-0.1918 -0.1034 0.1602 0.8604];

bu1 = [0.0045 0.0044];
bu2 = [0.1001 0.0100];
bu3 = [0.0003 -0.0136];
bu4 = [-0.0051 0.0936];

bw1 = [0.0953 0 0];
bw2 = [0.0145 0 0];
bw3 = [0.0862 0 0];
bw4 = [-0.0011 0 0];

cz1 = [1 0 -1 0];
cz2 = [0 0 0 0];
cz3 = [0 0 0 0];

dzu1 = [0 0];
dzu2 = [1 0];
dzu3 = [0 1];

A = [a1; a2; a3; a4];
Cz = [cz1; cz2; cz3];
Dzu = [dzu1; dzu2; dzu3];
Bu = [bu1; bu2; bu3; bu4];
Bw = [bw1; bw2; bw3; bw4];


n = size(A,1);
p = size(Bu,1);
q = size(Bu,2);
r = size(Dzu,1);
s = size(Dzu,2);
t = size(Cz, 1);

P = sdpvar(n,n, 'symmetric');
W = sdpvar(t,t, 'symmetric');
%L = sdpvar(q,n);
X1 = sdpvar(2,2);
X2 = sdpvar(2,2);
X = [X1 zeros(2); zeros(2) X2];
%X = sdpvar(n,n);
L1 = sdpvar(1,2);
L2 = sdpvar(1,2);
L = [L1 zeros(1,2); zeros(1,2) L2];



lmi1_11 = W;
lmi1_12 = Cz*X+Dzu*L;
lmi1_21 = lmi1_12';
lmi1_22 = X+X'-P;

LMI1 = [lmi1_11 lmi1_12; lmi1_21 lmi1_22];

b11 = P;
b12 = A*X+Bu*L;
b13 = Bw;
b21 = b12';
b22 = X+X'-P;
b23 = zeros(n,size(Bw,2));
b31 = b13';
b32 = b23';
b33 = eye(size(Bw,2));

lmi2_1 = [b11 b12 b13];
lmi2_2 = [b21 b22 b23];
lmi2_3 = [b31 b32 b33];

LMI2 = [lmi2_1; lmi2_2; lmi2_3];

mu = 0.41;

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

