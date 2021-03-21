function [Bs,Vs]=dttoperators(n,dttstr,use_norm)
% DTTOPERATORS returns sparse operators and their associated eigenvalues 
% given DTT type and transform length
% 
% [Bs,Vs]=dttoperators(n,dttstr,use_norm)
% 
% Bs: B operators
% Vs: a matrix, where each column contains the eigenvalues of Bs(:,:,1)
% n: transform length
% dttstr: string that indicates DTT type, 'dct1' to 'dct8' and 'dst1' to
%     'dst8'
% use_norm: set to 1 to normalize the scale, i.e., some Bs(:,:,1) and 
%     Bs(:,:,end) will be eye(n) or -eye(n) or their flipped versions 
%     instead of the scaled versions when use_norm=1
% 
% 20200413
if nargin<3 || isempty(use_norm)
    use_norm=0;
end

%% Operators
% number of operators
switch dttstr
    case 'dst1'
        m=n+2;
    case {'dct2', 'dct8', 'dst2', 'dst5', 'dst6', 'dst7'}
        m=n+1;
    case {'dct3', 'dct4', 'dct5', 'dct6', 'dct7', 'dst3', 'dst4', 'dst8'}
        m=n;
end

% BC
switch dttstr
    case {'dct2', 'dct4', 'dct6', 'dct8'}
        type1=1;
        a=1;
    case {'dst2', 'dst4', 'dst6', 'dst8'}
        type1=1;
        a=-1;
    case {'dct1', 'dct3', 'dct5', 'dct7'}
        type1=2;
        a=1;
    case {'dst1', 'dst3', 'dst5', 'dst7'}
        type1=3;
        a=-1;
end
switch dttstr
    case {'dct2', 'dct5', 'dst4', 'dst7'}
        type2=1;
        b=1;
    case {'dct4', 'dct7', 'dst2', 'dst5'}
        type2=1;
        b=-1;
    case {'dct1', 'dct6', 'dst3', 'dst8'}
        type2=2;
        b=1;
    case {'dct3', 'dct8', 'dst1', 'dst6'}
        type2=3;
        b=-1;
end

% generate rectangles
Bs=zeros(n,n,m);
Bs(:,:,1)=2*eye(n);
for i=1:m-1
    Bs(:,:,i+1)=dttrectmtx(n,i,type1,type2,a,b);
end

%% eigenvalues
ls=0:m-1;
switch dttstr
    case {'dct1', 'dct2', 'dct5', 'dct6'}
        js=0:n-1;        % j-1
    case {'dst1', 'dst2', 'dst5', 'dst6'}
        js=1:n;          % j
    case {'dct3', 'dct4', 'dct7', 'dct8', 'dst3', 'dst4', 'dst7', 'dst8'}
        js=0.5:1:n-0.5;  % j-1/2
end
switch dttstr
    case {'dct2', 'dct3', 'dct4', 'dst2', 'dst3', 'dst4'}
        nn=n;
    case 'dct1'
        nn=n-1;
    case 'dst1'
        nn=n+1;
    case {'dct5', 'dct6', 'dct7', 'dst8'}
        nn=n-0.5;
    case {'dst5', 'dst6', 'dst7', 'dct8'}
        nn=n+0.5;
end
Vs=2*cos(js'*ls*pi/nn);

%% normalization
if use_norm
    Bs(:,:,1)=Bs(:,:,1)/2;
    Vs(:,1)=Vs(:,1)/2;
    switch dttstr
        case {'dct1', 'dct2'}
            Bs(:,:,end)=Bs(:,:,end)/2;
            Vs(:,end)=Vs(:,end)/2;
        case {'dst1', 'dst2'}
            Bs(:,:,end)=-Bs(:,:,end)/2;
            Vs(:,end)=-Vs(:,end)/2;
    end
end