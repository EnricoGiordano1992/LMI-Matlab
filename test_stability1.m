%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test with lyap discrete time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = [1.844 -0.8505;1 0];

setlmis([]);
p = lmivar(1,[2 1]);


lmiterm([1 1 1 p],a1,a1');     % LMI #1 
lmiterm([1 1 1 p],-1,1);     % LMI #1 

lmiterm([-2 1 1 p],1,1);         % LMI #4: P 

lmis = getlmis;

[tmin,xfeas] = feasp(lmis)

P_Ys = dec2mat(lmis,xfeas,p);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a1 = [0.999];

setlmis([]);
p = lmivar(1,[1 1]);


lmiterm([1 1 1 p],a1,a1');     % LMI #1 
lmiterm([1 1 1 p],-1,1);     % LMI #1 

lmiterm([-2 1 1 p],1,1);         % LMI #4: P 

lmis = getlmis;

[tmin,xfeas] = feasp(lmis)

P_Ym = dec2mat(lmis,xfeas,p);
