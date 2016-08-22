%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test inequality 28 of paper   (KYP)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   P > 0
%
%   [I 0; A B; C D]'*[0 P 0; P 0 0; 0 0 R]*[I 0; A B; C D] < 0
% writed in:
%   [A'P+PA PB-C'; B'P-C; -D'-D] < 0 
%

a1 = [1.865 -0.8709;1 0];
b1 = [0.0009766; 0];
c1 = [0.000764 0.0007296];
d1 = 0;

setlmis([]);
p = lmivar(1,[2 1]);
lmiterm([1 1 1 p],1,a1,'s');    %A'P + PA   [1]
lmiterm([1 1 2 p],1,b1);        %PB         [2]
lmiterm([1 1 2 1],-1,c1');      %PB - C'    [2]
% NB: third element must be omitted because it's specular
%lmiterm([1 2 1 p],b1',1);       %B'P        [3]
%lmiterm([1 2 1 1],-1,c1);       %B'P -C     [3]
lmiterm([1 2 2 1],-1,d1');      %-D'        [4]
lmiterm([1 2 2 1],-1,d1);       %-D'-D      [4]

lmiterm([-2 1 1 p],1,1)

lmis = getlmis;

[tmin,xfeas] = feasp(lmis)

P_Ys = dec2mat(lmis,xfeas,p);
                                                     