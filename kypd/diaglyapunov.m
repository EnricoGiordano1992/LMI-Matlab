function X = diaglyapunov(ua,uainv,ta,C)
%DIAGLYAPUNOV  Solve continuous-time Lyapunov equations.
%
%   X = DIAGLYAPUNOV(ua,uainv,ta,C) 
%   solves the special form of the Lyapunov matrix equation:
%
%           A*X + X*A' = -C
%
%   where A = ua*ta*uainv. The following lines of code will generate 
%   ua, ta from A
%
%   [ua,ta] = eig(A); 
%   uainv = inv(ua);  

%
%   See also  LYAP.

%       A. Hansson 2002-08-30

[ma,na] = size(ta);
[mc,nc] = size(C);
  
% Transform C
ucu = uainv*C*uainv';

dummy = diag(ta)*ones(1,ma);
y = -ucu./(dummy + dummy');

% Find untransformed solution 
X = ua*y*ua';

X = real(X);

X = (X+X')/2;

% end diaglyapunov