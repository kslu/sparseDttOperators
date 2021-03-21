clear;
close all;

n=24;
tau=2;

[Bds,eig_bs]=dttoperators(n,'dct2');
eig_l=2-eig_bs(:,2);
eig_ls=(2-eig_bs(:,2)).^(0:n);
%%
% h=tau./(tau+eig_l);
%h=exp(-0.5*(eig_l-2).^2);
h=[ones(8,1);zeros(16,1)];
w1=[ones(7,1)*2;0*ones(2,1);ones(15,1)];
bd1=(eig_l(max(find(w1==2)))+eig_l(min(find(w1==0))))/2; % passband boundary
bd2=(eig_l(min(find(w1==1)))+eig_l(max(find(w1==0))))/2; % stopband boundary

k=4;
%%
[cls_u,hls_u,erls_u,mls_u]=pgffit(eig_l,h,k,[],0);
[cls_w,hls_w,erls_w,mls_w]=pgffit(eig_l,h,k,w1,0);
[cm_u,hm_u,erm_u,mm_u]=pgffit(eig_l,h,k,[],1);
[cm_w,hm_w,erm_w,mm_w]=pgffit(eig_l,h,k,w1,1);

figure;
hold on; grid on;
plot(eig_l, h, '-*', 'linewidth', 3, 'markersize', 10);
plot(eig_l, hls_u, '-.o', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hls_w, '-.s', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hm_u, ':+', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hm_w, ':x', 'linewidth', 2.5, 'markersize', 10);
plot([bd1, bd1], [-0.5, 2], '--k');
plot([bd2, bd2], [-0.5, 2], '--k');

legend('Desired', 'LS, unweighted', 'LS, weighted', 'Minimax, unweighted',...
    'Minimax, weighted');

xlabel('Graph frequency');
ylabel('Frequency response');
set(gca, 'fontsize', 20);
set(gcf, 'position', [0,0,600,400]);

text(0.05,1.45,'passband','fontsize',16);
text(0.7,1.45,'transition','fontsize',16);
text(0.7,1.35,'band','fontsize',16);
text(1.4,1.45,'stopband','fontsize',16);
axis([0,4,-0.4,1.5]);

hgexport(gcf, 'example_ls_vs_minimax.png', hgexport('factorystyle'), ...
    'format', 'png');