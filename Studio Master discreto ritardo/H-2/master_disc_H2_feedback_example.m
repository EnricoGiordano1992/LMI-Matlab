Mm = 0.64;
Bm = 12;

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

cy1 = [1 0 0 0];
cy2 = [0 0 1 0];

dyw1 = [0 1 0];
dyw2 = [0 0 1];

A = [a1; a2; a3; a4];
Cz = [cz1; cz2; cz3];
Dzu = [dzu1; dzu2; dzu3];
Bu = [bu1; bu2; bu3; bu4];
Bw = [bw1; bw2; bw3; bw4];
Cy = [cy1; cy2];
Dyw = [dyw1; dyw2];


n = size(A,1);
p = size(Bu,1);
q = size(Bu,2);
r = size(Dyw,1);
s = size(Dyw,2);
t = size(Cy,1);
u = size(Cy,2);
v = size(Cz,1);

P = sdpvar(n,n, 'symmetric');
H = sdpvar(n,n, 'symmetric');
W = sdpvar(v,v, 'symmetric');
L = sdpvar(q,n); 
F = sdpvar(n,t); 
Q = sdpvar(n,n); 
R = sdpvar(q,t); 
S = sdpvar(n,n); 
J = sdpvar(n,n); 
X = sdpvar(n,n); 
Y = sdpvar(n,n); 

lmi11 = [W Cz*X+Dzu*L Cz+Dzu*R*Cy];
lmi12 = [(Cz*X+Dzu*L)' X+X'-P eye(size(S))+S'-J];
lmi13 = [(Cz+Dzu*R*Cy)' (eye(size(S))+S'-J)' Y+Y'-H];

LMI1 = [lmi11; lmi12; lmi13];

b11 = P;
b12 = J;
b13 = A*X+Bu*L;
b14 = A+Bu*R*Cy;
b15 = Bw+Bu*R*Dyw;

b21 = b12';
b22 = H;
b23 = Q;
b24 = Y*A+F*Cy;
b25 = Y*Bw+F*Dyw;

b31 = b13';
b32 = b23';
b33 = X+X'-P;
b34 = eye(size(S))+S'-J;
b35 = zeros(size(S,1), size(b15,2));

b41 = b14';
b42 = b24';
b43 = b34';
b44 = Y+Y'-H;
b45 = zeros(size(b44,1), size(b15,2));

b51 = b15';
b52 = b25';
b53 = b35';
b54 = b45';
b55 = eye(size(b54,1), size(b15,2));


lmi21 = [b11 b12 b13 b14 b15];
lmi22 = [b21 b22 b23 b24 b25];
lmi23 = [b31 b32 b33 b34 b35];
lmi24 = [b41 b42 b43 b44 b45];
lmi25 = [b51 b52 b53 b54 b55];

LMI2 = [lmi21; lmi22; lmi23; lmi24; lmi25];

res = [X >= 0, Y >= 0, S >= 0, R >= 0, LMI1 >=0 LMI2 >= 0, trace(W) <= mu, Dzw+Dzu*R*Dyw==0];

diagnostics = solvesdp(res);
disp(diagnostics.problem);
if diagnostics.problem == 0
 disp('Feasible')

 [U,V] = qr(-value(Y)*value(X)+eye(size(value(Y)*value(X))));
 m12 = -inv(V)*value(Y)*Bu;
 mat1 = [inv(V) m12; zeros(size(Bu,2),n) eye(size(Bu,2))];
 mat2 = [value(Q)-value(Y)*A*value(X) value(F); value(L) value(R)];
 mat3 = [inv(value(U)) zeros(n,size(Cy,1)); -Cy*value(X)*inv(value(U)) eye(size(Cy,1))];
 K = mat1*mat2*mat3;

 %state space realization of the controller
 Ac = K(1:n,1:n);
 Bc = K(1:n,n+1:end);
 Cc = K(n+1:end,1:n);
 Dc = K(end,end);

elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end
