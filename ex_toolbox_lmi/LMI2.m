a1 = [1 2;1 -3]

setlmis([]) 
p = lmivar(1,[2 1])

lmiterm([1 1 1 p],1,a1,'s')     % LMI #1 
lmiterm([-2 1 1 p],1,1)         % LMI #4: P 

lmis = getlmis

[tmin,xfeas] = feasp(lmis)

P = dec2mat(lmis,xfeas,p)
