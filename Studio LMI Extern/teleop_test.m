Ms = 0.64;
Bs = 12;
Mm = 0.64;
Bm = 12;

Ke = 0;%0.01;
Me = 0;%1;
De = 0;%0.01; %damping
Xe = 9;

A_master = [0 1; 0 -Bm/Mm];
B_master = [0; 1/Mm];
C_master = [1 0; 0 1];
D_master = [0; 0];

A_slave = [0 1; 0 -Bs/Ms];
B_slave = [0; 1/Ms];
C_slave = [1 0; 0 1];
D_slave = [0; 0];

A = [0 1 0 0; 0 (-Bm/Mm) 0 0; 0 0 0 1; 0 0 0 -(Bs/Ms)];
B = [0 0; 1/Mm 0; 0 0; 0 1/Ms];
B0 = [0 0; 1/Mm 0; 0 0; 0 1/Ms];
C = [0 1 0 0; 0 0 0 1];
C0 = [0 1 0 -1];

n = size(A,1);
q = size(B,2);
r = size(C0,1);

E = zeros(r,q);
D0 = E;



disp('H2:');
output = h2_lmi_c(A,B,C)
if output.h2 > 0
    K = output.h2;
else
    disp('infeasible');
end


%disp('Hinf:');
%output = hinf_norm_c_yal(A,B,C,D)
%if output.feas > 0
%    K = output.P;
%else
%    disp('infeasible');
%end
