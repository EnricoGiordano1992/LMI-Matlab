Mm = 0.64;
Bm = 12;
A = [0 1; 0 -Bm/Mm];
B = [0; 1/Mm];
C = [1 0; 0 1];
D = [0; 0];


disp('H2:');
output = h2_lmi_c(A,B,C)
if output.h2 > 0
    K = output.P;
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
