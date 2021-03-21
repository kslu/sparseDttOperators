function [x,ha,c,xc,er]=mpgffit(P,h,k,wf,use_minimax,l0_method)
% MPGFFIT computes a multivariate polynomial graph filter (MPGF) 
% in terms of coefficients x that approximately minimizes ||P*x - h||^2
% 
% [x,ha,c,xc,er]=mpgfomp(P,h,k,wf,use_minimax,l0_method)
% 
% Input arguments
%   P: matrix of eigenvalues, each column of which contains the eigenvalues
%      of the associated polynomial term
%   h: desired frequency response
%   k: number of multivariate polynomial terms to use (not applicable if
%      l0_method='leastsquares')
%   wf: frequency weights
%   use_minimax: 1 for minimax criterion 
%                0 for least squares fitting (default)
%   l0_method: 'exhaust' for exhaust
%              'omp' for orthogonal matching pursuit (default)
%              'leastsquares' for least squares solution
% 
% Output arguments
%   x: solution
%   ha: A*x, approximation of h
%   c: indices of nonzero coefficients in x
%   xc: nonzero coefficients in x
%   er: error
% 
% KS Lu
% 20200715

[n,m]=size(P);

if nargin<6 || isempty(l0_method)
    l0_method='omp';
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
id_w0=find(wf==0);

switch l0_method
    case 'omp'
        
        tol=1e-8;
        res=h;
        
        idx_rest=1:m;
        idx_picked=zeros(1,k);

        % normalize columns of A (needed for correlation comparison)
        Pn=normc(P);

        for i=1:k
            % find the atom with maximal correlation to the residue
%             [i,idx_rest]
            corrs=abs(Pn(:,idx_rest)'*res);
            [cmax,idmax]=max(corrs);
            if cmax<tol
                idx_picked=idx_picked(1:i-1);
                break;
            end
            idx_picked(i)=idx_rest(idmax);
            idx_rest(idmax)=[];  % remove the entry that is chosen

            % sub-matrix of Pi
            subm_pi=P(:,idx_picked(1:i));
            
            % fit weights
            if use_minimax
                
                % constraints
                f=[zeros(1,i), 1]';
                A=[  subm_pi, -ones(n,1)./(wf+eps);
                    -subm_pi, -ones(n,1)./(wf+eps)];
                b=[h;-h];
                
                % remove constraints associated to the transition band 
                % (zero entries in wf)
                A([id_w0, n+id_w0],:)=[];
                b([id_w0, n+id_w0],:)=[];

                xc=linprog(f,A,b,[],[],[],[],options);
                xc=xc(1:end-1);
                res=h-subm_pi*xc;
            else
                xc=(subm_pi'*diag(wf+eps)*subm_pi)\subm_pi'*diag(wf)*h;
                res=h-subm_pi*xc;
            end
        end

        x=zeros(m,1);
        x(idx_picked)=xc;
        c=idx_picked;
        ha=h-res;

    case 'exhaust'
        cur_comb=1:k;
        opterr=1e10;
        
        while ~isempty(cur_comb)
            
            % sub-matrix of Pi
            subm_pi=P(:,cur_comb);
            % skip if subm_pi'*subm_pi is singular
            if min(svd(subm_pi'*subm_pi))<1e-4
                cur_comb=nextcomb(m,cur_comb);
                continue;
            end
                
            if use_minimax

                % constraints
                f=[zeros(1,k), 1]';
                A=[  subm_pi, -ones(n,1)./(wf+eps);
                    -subm_pi, -ones(n,1)./(wf+eps)];
                b=[h;-h];
                
                % remove constraints associated to the transition band 
                % (zero entries in wf)
                A([id_w0, n+id_w0],:)=[];
                b([id_w0, n+id_w0],:)=[];

                xa=linprog(f,A,b,[],[],[],[],options);
                xa=xa(1:end-1);
                
            else
                xa=(subm_pi'*diag(wf+eps)*subm_pi)\subm_pi'*diag(wf)*h;
            end
            
            if norm(subm_pi*xa-h)<opterr
                xc=xa;
                opterr=norm(subm_pi*xa-h);
                x=zeros(m,1); x(cur_comb)=xa;
                ha=P*x;
                c=cur_comb;
            end
            cur_comb=nextcomb(m,cur_comb);
        end

    case 'leastsquares'
        if use_minimax
            % constraints
            f=[zeros(1,m), 1]';
            A=[  P, -ones(n,1)./(wf+eps);
                -P, -ones(n,1)./(wf+eps)];
            b=[h;-h];

            % remove constraints associated to the transition band 
            % (zero entries in wf)
            A([id_w0, n+id_w0],:)=[];
            b([id_w0, n+id_w0],:)=[];

            x=linprog(f,A,b,[],[],[],[],options);
            x=x(1:end-1);
        else
            x=(P'*diag(wf+eps)*P)\P'*diag(wf)*h;
        end
        ha=P*x;
        c=find(x~=0);
        xc=x(c);
end

er=wf'*(ha-h).^2;