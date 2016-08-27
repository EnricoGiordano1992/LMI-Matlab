n = size(A,1);

t = sdpvar(1,1);
P = sdpvar(n,n);

F = [kyp(A,B,P,blkdiag(C'*C,-t)) <= 0]

sol = optimize(F,t);
