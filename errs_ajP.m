close all
clear all
fprintf('Loading points \n')
p  = load('points_r1.dat');
ind = 1;
for i = 1 :size(p,1)
        if p(i,1) <= 5
                pp(ind,:) = p(i,:);
                ind = ind + 1;
        end
end
p = pp;
P = 5;
fprintf('Loading solution \n')
ut0 = load('solt_r5-h0.05.dat');

CreateVideo = 0;

ajP = [0.01 0.1 0.5 1.0];

t = linspace(0,10,101);
Markers = {'+','o','*','x','v','d','^','s','>','<'};

figure

errt = zeros(size(ut0,1),length(ajP));
for j = 1 : size(errt,2)
	ut6   = load(['solt_r5-P',num2str(P),'-nL20.0-aP',num2str(ajP(j)),'.dat'])/0.05;
	for tt = 1 : size(ut6,1)
		errt(tt,j) = norm(ut0(tt,:)-ut6(tt,:))/(sqrt(9.7*3));
	end
	semilogy(t(1:size(ut6,1)),errt(1:size(ut6,1),j),strcat('-',Markers{j}))
	hold on
end

legend('$a_P = 0.01$','$a_P = 0.1$','$a_P = 0.5$','$a_P = 1$','$a_P = 3$','$a_P = 5$','$a_P = 9$','interpreter','latex','fontsize',16,'location','southeast')



