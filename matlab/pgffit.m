function [c,hfit,er,mevs]=pgffit(B,h,k,wf,use_minimax,mevs)
% PGFFIT fits a polynomial graph filter given a shift B
%
% [c,hfit,er,mevs]=pgffit(B,h,k,wf,use_minimax)
% or [c,hfit,er,mevs]=pgffit(ev,h,k,wf,use_minimax)
% 
% Input arguments
%   B: shift operator (Laplacian or polynomial of L)
%   (or) ev: eigenvalues of B
%   h: desired frequency response
%   k: polynomial order
%   wf: frequency weighting
% 
% Output arguments
%   c: filter coefficients
%   hfit: approximate frequency response
%   er: error
%   mevs: central graph frequency, (min(ev)+max(ev))/2
%
% 20200427

h=h(:);
n=numel(h);

if size(B,1)==1 || size(B,2)==1
    evs=B;
else
    evs=eig(B);
end

if nargin<5 || isempty(use_minimax)
    use_minimax=0;
end

if use_minimax
    options = optimoptions('linprog','Display','none');
end

% frequency weighting
if nargin<4 || isempty(wf)
    wf=ones(n,1);
%     wf=exp(0.2*(n:-1:1))';
end
%wf=max(wf,0.01);
id_w0=find(wf==0);

% centering the eigenvalues to improve stability
if nargin<6 || isempty(mevs)
    mevs=2;
%     mevs=(min(evs)+max(evs))/2;
end
cevs=evs-mevs;

if use_minimax==1
    Lamb=(cevs*ones(1,k+1)).^(ones(n,1)*(0:k));
    f=[zeros(1,k+1), 1]';
    
    % constraints
    A=[Lamb, -ones(n,1)./(wf+eps); -Lamb, -ones(n,1)./(wf+eps)];
    b=[h;-h];
    
    % remove constraints associated to the transition band
    % (zero entries in wf)
    A([id_w0, n+id_w0],:)=[];
    b([id_w0, n+id_w0],:)=[];
    
    ge=linprog(f,A,b,[],[],[],[],options);
    c=ge(1:end-1);
    hfit=Lamb*c;
else
    Lamb=(cevs*ones(1,k+1)).^(ones(n,1)*(0:k));
    c=(Lamb'*diag(wf+eps)*Lamb)\Lamb'*diag(wf)*h;
    %c=linsolve(Lamb,h);
    hfit=Lamb*c;
end

er=wf'*(hfit-h).^2;

%% Filter matrix construction (deprecated)
% % A. direct polynomial construction 
% H=0*B;
% for i=0:k
%     H=H+c(i+1)*B^i;
% end
% % B. polynomial construction with one order at a time (faster)
% H=c(k+1)*eye(n);
% for i=k:-1:1
%     H=H*(B-mevs*eye(n))+c(i)*eye(n);
% end