Ms = 0.61;
Bs = 11;
Kv = 40;
Kp = 40;
Mm = 0.64;
Bm = 12;

s = tf('s');

A = [0 1 0 0; 0 (-Bm/Mm) 0 0; 0 0 0 1; 0 0 0 -(Bs/Ms)];
B = [0 0; 1/Mm 0; 0 0; 0 1/Ms];
B0 = [0 0; 1/Mm 0; 0 0; 0 1/Ms];
C0 = [1 0 -1 0; 0 1 0 -1];
C = [0 1 0 0; 0 0 0 1];
D = [0 0; 0 0];

n = size(A,1)
p = size(B,2)
P = sdpvar(size(A,1),size(A,1), 'symmetric');

F = [P >= 0, [[eye(n,n) zeros(n,p);A B; C D]'*[zeros(n,n) P zeros(n,p); P zeros(n,n) zeros(n,p); zeros(p,n) zeros(p,n) eye(p,p)]*[eye(n,n) zeros(n,p);A B; C D]] >= 0];
%F = [P >= 0, [A'*P+P*A P*B-C'; B'*P-C, -D'-D] <=0]
diagnostics = solvesdp(F)
disp(diagnostics.problem)
if diagnostics.problem == 0
 disp('Feasible')
 solution = value(P)
elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end
