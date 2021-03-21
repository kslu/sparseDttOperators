function y=pgfeval(B,c,x,mev)
% PGFEVAL evaluates a polynomial graph filter given a shift B, PGF
% coefficients c, input signal x, and center graph frequency mev
%
% y=pgfeval(B,c,x,mev)
% 
% Input arguments
%   B: shift operator (Laplacian or polynomial of L)
%   c: PGF coefficients
%   x: input signal
%   mev: center graph frequency
% 
% Output arguments
%   y: output signal
%
% 20200801
c=c(:);
n=size(B,1);
k=numel(c);

% y0 = c_k * x
% y1 = (c_k * L + c_k-1) * x
%    = L * y0 + c_k-1 * x
% y2 = (c_k * L^2 + c_k-1 * L + c_k-2) * x
%    = L * y1 + c_k-2 * x
%    ...
% yk = (c_k * L^k + ... + c_1 * L + c_0) * x
%    = L * y_k-1 + c_0 * x
y=c(k)*x;
B0=B-mev*eye(n);
for i=k-1:-1:1
    y=B0*y+c(i)*x;
end