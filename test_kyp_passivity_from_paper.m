n = size(A,1)
p = size(B,2)
P = sdpvar(size(A,1),size(A,1), 'symmetric');

%F = [P >= 0, [[eye(n,n) zeros(n,p);A B; C D]'*[zeros(n,n) P zeros(n,p); P zeros(n,n) zeros(n,p); zeros(p,n) zeros(p,n) eye(p,p)]*[eye(n,n) zeros(n,p);A B; C D]] >= 0];
F = [P >= 0, [A'*P+P*A P*B-C'; B'*P-C, -D'-D] <=0]
diagnostics = solvesdp(F)
disp(diagnostics.problem)
if diagnostics.problem == 0
 disp('Feasible')
elseif diagnostics.problem == 1
 disp('Infeasible')
else
 disp('Something else happened')
end
