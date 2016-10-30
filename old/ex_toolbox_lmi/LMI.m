A = 1/2*eye(2);
C = [1 1];
B = C';

setlmis([]) 
X=lmivar(1,[2 1]) 
S=lmivar(1,[1 0])

% 1st LMI 
lmiterm([1 1 1 X],1,A,'s') 
lmiterm([1 1 1 S],C',C) 
lmiterm([1 1 2 X],1,B) 
lmiterm([1 2 2 S],-1,1)

% 2nd LMI 
lmiterm([-2 1 1 X],1,1)

% 3rd LMI
lmiterm([-3 1 1 S],1,1) 
lmiterm([3 1 1 0],1)

lmis = getlmis;

[tmin,xfeas] = feasp(lmis)

P = dec2mat(lmis,xfeas,X)