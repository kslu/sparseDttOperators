function y=mpgfeval(Bs,ps,coeffs,x)
% MPGFEVAL returns output of MPGF 
% 
% H=mpgfeval(Bs,ps,xc)
% 
% Input arguments:
%   Bs: nxnxk tensor of k graph operators
%   ps: power list
%   coeffs: filter coefficients
%   x: input signal
%
% Output argument:
%   y: graph filter matrix
% 
% Example:
%    % low-pass filter
%    h0=[2,2,1,1,1,1,0,0]';
%    % generate DCT-II operators
%    [Bs,Vs]=dttoperators(8,'dct2');
%    % find eigenvalues and power list for degree-2 polynomial terms
%    [Ys,powers]=mvpolyexpand(Vs,2);
%    % solve optimal approximate filter with degree 2
%    [coeffs,hb,cb,xc,erb]=mpgfexhaust(Ys,h0,3);
%    % compute the output of x=1:8;
%    y=mpgfeval(Bs,powers,coeffs,x);
%
% KS Lu
% 20200817

cb=find(coeffs);
xc=coeffs(cb);
y=0*x;
for b=1:numel(cb)
    pcidx=find(ps(:,cb(b)));
    yc=x;
    for ii=1:numel(pcidx)
        yc = Bs(:,:,pcidx(ii))^ps(pcidx(ii),cb(b))*yc;
    end
    y=y+xc(b)*yc;
end