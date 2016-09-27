% F = ellipfilt (n, xi, damp, wF)
%
% computes an elliptic-type body bending filter of nth order
% Inputs
%   n     order of the filter
%   xi    (> 1), defines the ripple, large values give small ripple
%             small values > 1, high ripple but a sharp change at wF
%   damp  damping in dB of filter
%   wF   filter cut-off frequency

% Author: Anders Helmesson, Dept of EE, Linköpings universitet
% email: andersh@isy.liu.se
% http://users.isy.liu.se/rt/andersh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = ellipfilt (n, xi, damp, wF)

if nargin < 3,
  disp ('usage:  F = ellipfilt (n, xi, damp, wF)');
  return;
end

if nargin < 4,
% default values
  wF = 1;
end;

k = 10^(damp/20);
b = real (poly (ellipratz (n, xi, k)));
ww = wF.^(0:length(b)-1);
F = ss (tf ((b.*ww)/sqrt(k), b(end:-1:1).*ww));

% balance the system
F = balreal (F);

% perform a similarity transformation to obtain c positive
% this form is unique if all Hankel singular values are distinct
F = ss2ss (F, diag (sign (F.c)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% local function related to elliptic rational functions

function z = ellipratz (n, xi, k)
[sn, cn, dn] = ellipj (ellipke (1/xi^2)*(2*(1:n)-1)/n, 1/xi^2);
z = (cn./dn)/sqrt(xi);
if nargin > 2,
  b = poly (z);
  z = j * roots (b + j*sqrt(k)*b(end:-1:1));
  i = find (real(z) > 0);
  z(i) = -conj(z(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function k = ellipke (x)
% y = [1, sqrt(1-x)]
% while abs(diff(y)) > eps, y = [0.5*sum(y), sqrt(prod(y))]; end
% k = pi/sum(y);
