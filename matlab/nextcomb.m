function nextc=nextcomb(n,c)
% NEXTCOMB find next combination in terms of lexicographic order, 
%     given n and the current combination
% 
% nextc=nextcomb(n,c)
% 
% n: n for "n choose k"
% c:    current combination. It must have size k and is sorted in increasing order
% 
% KS Lu
% 20170401
%
k=numel(c);
if all(c==n-k+1:n)
    %fprintf('The input array is the last combination\n');
    nextc=[];
    return;
end

nextc=c;
m=find(c<n-k+1:n,1,'last');
nextc(m)=c(m)+1;
if m<k
    nextc(m+1:k)=c(m)+2:c(m)+k-m+1;
end