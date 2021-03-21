function M=dttrectmtx(n,i,type1,type2,a,b,c)
% DTTRECTMTX returns a rectangular matrix as DCT/DST sparse shift
% operators
% 
% M=dttrectmtx(n,i,type1,type2,a,b,c)
% i: index for diagonal, ranging from 0 to n
% type1: upper-left side of rectangle (use i = 3 for example)
%     pattern with type=1:  0 0 1 1 0 0 0 0
%                           0 1 0 0 1 0 0 0
%                           1 0 0 0 0 1 0 0
%                           1 0 0 0 0 0 1 0
%                           0 1 0 0 0 0 0 1
%                           0 0 1 0 0 0 0 1
%                           0 0 0 1 0 0 1 0
%                           0 0 0 0 1 1 0 0
%     pattern with type=2:  0 0 0 c 0 0 0 0
%                           0 0 1 0 1 0 0 0
%                           0 1 0 0 0 1 0 0
%                           c 0 0 0 0 0 1 0
%                           0 1 0 0 0 0 0 1
%                           0 0 1 0 0 0 0 1
%                           0 0 0 1 0 0 1 0
%                           0 0 0 0 1 1 0 0, where c=sqrt(2) by default
%     pattern with type=3:  0 a 0 1 0 0 0 0
%                           a 0 0 0 1 0 0 0
%                           0 0 0 0 0 1 0 0
%                           1 0 0 0 0 0 1 0
%                           0 1 0 0 0 0 0 1
%                           0 0 1 0 0 0 0 1
%                           0 0 0 1 0 0 1 0
%                           0 0 0 0 1 1 0 0
% type2: lower-right side of rectangle (similar to type1)
% a: value of the top-left side (default: 1)
% b: value of the bottom-right side (default: 1)

if nargin<7 || isempty(c)
    c=sqrt(2);
end
if nargin<6 || isempty(b)
    b=1;
end
if nargin<5 || isempty(a)
    a=1;
end

if type1==3 && type2==3 && i==n+1
    M=zeros(n);
else
    M=diag(ones(n-i,1),i)+diag(ones(n-i,1),-i);
end

switch type1
    case 1
        M=M+fliplr(a*diag(ones(i,1),n-i));
    case 2
        M=M+fliplr(a*diag(ones(i+1,1),n-i-1));
        M(1,i+1)=c;
        M(i+1,1)=c;
    case 3
        M=M+fliplr(a*diag(ones(i-1,1),n-i+1));
end

switch type2
    case 1
        M=M+fliplr(b*diag(ones(i,1),-n+i));
    case 2
        M=M+fliplr(b*diag(ones(i+1,1),-n+i+1));
        M(n-i,n)=c;
        M(n,n-i)=c;
    case 3
        M=M+fliplr(b*diag(ones(i-1,1),-n+i-1));
end

% special case
if type1==2 && type2==2 && i==n-1
    M=2*fliplr(eye(n));
end