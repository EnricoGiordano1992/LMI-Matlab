function X = lyapunov(ua, ta, ub, tb, C)
%LYAPUNOV  Solve continuous-time Lyapunov equations.
%
%   X = LYAPUNOV(ua,ta,ub,tb,C) 
%   solves the special form of the Lyapunov matrix equation:
%
%           A*X + X*A' = -C
%
%   where A = ua*ta*ua'. The following lines of code will generate 
%   ua, ta, ub, tb from A
%
%   [ua,ta] = schur(A); 
%   [ua,ta] = rsf2csf(ua,ta);
%   %Schur decomposition of A' can be calculated from that of A.
%   j = ma:-1:1;
%   ub = ua(:,j);
%   tb = ta(j,j)';
%   %Check all combinations of ta(i,i)+tb(j,j) for zero
%   p1 = diag(ta).'; % Use .' instead of ' in case A and B are not real
%   p1 = p1(ones(mb,1),:);
%   p2 = diag(tb);
%   p2 = p2(:,ones(ma,1));
%   sum = abs(p1) + abs(p2);
%   if any(any(sum == 0)) | any(any(abs(p1 + p2) < 1000*eps*sum))
%     error('Solution does not exist or is not unique.');
%   end

%
%   See also  LYAP.

%       A. Hansson 2002-08-29

[ma,na] = size(ta);
[mb,nb] = size(tb);
[mc,nc] = size(C);
  
% Transform C
ucu = -ua'*C*ub;

% Solve for first column of transformed solution
y = zeros(ma,mb);
ema = eye(ma);
y(:,1) = (ta+ema*tb(1,1))\ucu(:,1);

% Solve for remaining columns of transformed solution
for k=2:mb,
   km1 = 1:(k-1);
   y(:,k) = (ta+ema*tb(k,k))\(ucu(:,k)-y(:,km1)*tb(km1,k));
end

% Find untransformed solution 
X = ua*y*ub';

X = real(X);

X = (X+X')/2;

% end lyapunov