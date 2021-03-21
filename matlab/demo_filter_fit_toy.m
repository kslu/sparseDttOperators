close all;
clear;

n=12;
M=4;

[Bds,eig_bs]=dttoperators(n,'dct2');
eig_l=2-eig_bs(:,2);
eig_ls=(2-eig_bs(:,2)).^(0:n);
h1=1./(1+exp(4*(eig_l-2)));
h2=exp(-4*(eig_l-1).^2);

hp=zeros(n,M);
hm=zeros(n,M);
gm=zeros(n+1,M);
erm=zeros(1,M);
erp=zeros(1,M);
Hp=cell(1,M);

%% low-pass filter
for l=1:M
    [gm(:,l),hm(:,l),cm{l},xcm{l},erm(l)]=mpgffit(eig_bs,h1,l,[],[],'exhaust');
    [cp{l},hp(:,l),erp(l)]=pgffit(eig_l,h1,l-1);
end

cm3=cm{3};
str1=sprintf('MPGF, p_1(Z^{(%d)},Z^{(%d)})',cm3(2),cm3(3));
cm4=cm{4};
str2=sprintf('MPGF, p_1(Z^{(%d)},Z^{(%d)},Z^{(%d)})',cm4(2),cm4(3),cm4(4));

figure;
hold on; grid on;
plot(eig_l, h1, '-*', 'linewidth', 3, 'markersize', 10);
plot(eig_l, hp(:,3), '-.o', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hp(:,4), '-.s', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hm(:,3), ':+', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hm(:,4), ':x', 'linewidth', 2.5, 'markersize', 10);
legend('Desired', 'PGF, K=2', 'PGF, K=3', str1, str2);

xlabel('Graph frequency');
ylabel('Frequency response');
set(gca, 'fontsize', 20);
set(gcf, 'position', [0,0,600,400]);

% axis([0,4,-0.4,1.5]);
hgexport(gcf, 'example_mpgf_lp.png', hgexport('factorystyle'), ...
    'format', 'png');

%% band-pass filter
for l=1:M
    [gm(:,l),hm(:,l),cm{l},xcm{l},erm(l)]=mpgffit(eig_bs,h2,l,[],[],'exhaust');
    [cp{l},hp(:,l),erp(l)]=pgffit(eig_l,h2,l-1);
end
cm3=cm{3};
str3=sprintf('MPGF, p_1(Z^{(%d)},Z^{(%d)})',cm3(2),cm3(3));
cm4=cm{4};
str4=sprintf('MPGF, p_1(Z^{(%d)},Z^{(%d)},Z^{(%d)})',cm4(2),cm4(3),cm4(4));

figure;
hold on; grid on;
plot(eig_l, h2, '-*', 'linewidth', 3, 'markersize', 10);
plot(eig_l, hp(:,3), '-.o', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hp(:,4), '-.s', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hm(:,3), ':+', 'linewidth', 2.5, 'markersize', 10);
plot(eig_l, hm(:,4), ':x', 'linewidth', 2.5, 'markersize', 10);
legend('Desired', 'PGF, K=2', 'PGF, K=3', str3, str4);

xlabel('Graph frequency');
ylabel('Frequency response');
set(gca, 'fontsize', 20);
set(gcf, 'position', [0,0,600,400]);

hgexport(gcf, 'example_mpgf_bp.png', hgexport('factorystyle'), ...
    'format', 'png');