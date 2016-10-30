clc%, clear all
%clear('yalmip')
opt = sdpsettings;
opt = sdpsettings(opt,'verbose',0,'solver','sedumi');
s=tf('s'); n1=1; n2=2; n=2;
%for i=1:1:n1
%    Fi1(i,1)=1/(s+3)^(i-1);
%end
 %Environment
%for i=1:1:n2
%    Fi2(i,1)=1/(s+5)^(i-1);
%end

Fi1 = Ym
Fi2 = Ys

[A1,B1,C1,D1]=ssdata(Fi1);
n1 = size(A1,1)
K1 = sdpvar(n1,n1,'symmetric');
P1=sdpvar(size(A1,1),size(A1,2));
H1=[A1'*P1+P1*A1 P1*B1; B1'*P1 zeros(size(P1*B1,2),size(B1'*P1,1))];
M1=[C1*K1*C1' C1*K1*D1 ; D1'*K1*C1' D1'*K1*D1];
F=(H1+M1) >= 0; % Fi1*K1*Fi1 >= 0

[A2,B2,C2,D2]=ssdata(Fi2);
n2 = size(A2,1)
K2 = sdpvar(n2,n2,'symmetric');
P2=sdpvar(size(A2,1),size(A2,2));
H2=[A2'*P2+P2*A2 P2*B2; B2'*P2 zeros(size(P2*B2,2),size(B2'*P2,1))];
T2=[C2*K2*C2' C2*K2*D2 ; D2'*K2*C2' D2'*K2*D2];
F=[F,(H2+T2) >= 0]; % + Fi2*K2*Fi2 >= 0

K=[ zeros(size(blkdiag(K1,K2))) blkdiag(K1,K2); ...
blkdiag(K1,K2) zeros(size(blkdiag(K1,K2)))];
% Definining Nominal plant G. Parameters are taken from
% \cite{Willaert2011}
%Mm=0.64; Bm=3.4; Ms=0.61; Bs=11;
%Ym=1/(s*Mm+Bm); Ys=1/(s*Ms+Bs);
% Define Controllers
%Cm=1; Cs=50;
%G=(1/(s+Ym*Cm+Ys*Cs))*[Ym*(s+Ys*Cs) Ym*Cm*Ys ; Ys*Cs*Ym Ys*(s+Ym*Cm)];
%Gp=blkdiag(Fi1,Fi2,Fi1,Fi2)*[-G;eye(size(G,2))];
%[A,B,C,D]=ssdata(Gp);

P=sdpvar(size(A,1),size(A,2),'symmetric');
H=[A'*P+P*A P*B; B'*P zeros(size(P*B,2),size(B'*P,1))];
T=[C*K*C' C*K*D ; D'*K*C' D'*K*D];
F=[F,(H+T) >= -0.000001*eye(size(H))];
solution=solvesdp(F,[],opt)
report=solution.problem;
if report == 0
disp('Stable')
elseif report == 1
disp('Infeasible')
else
display('Hmm, something went wrong!');
yalmiperror(solution.problem)
end