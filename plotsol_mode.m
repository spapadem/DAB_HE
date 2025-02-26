clear all
close all

%h = 2.0;
h = 2.0;

Markers = {'+','o','*','x','v','d','^','s','>','<'};

f = 24;

Pmax = 15;
rmax =  5;
err_r = zeros(7,Pmax);
err_i = zeros(7,Pmax);
for r = 1 : rmax
	p = load(['points_r',num2str(r),'-h',num2str(h),'.0.dat']);
	fprintf('r = %2d\n',r)
	ure = load(['exact_f',num2str(f),'-r',num2str(r),'-h',num2str(h),'.0-real.dat']);
	uie = load(['exact_f',num2str(f),'-r',num2str(r),'-h',num2str(h),'.0-imag.dat']);
%	ure = ure/max(ure);
%	uie = uie/max(uie);
	ue = ure + 1i *uie;
	
	for P = 1 : Pmax
		fprintf('Loading P = %2d of %2d\n',P,Pmax)
		ur =  load(['sol_f',num2str(f),'-r',num2str(r),'-P',num2str(P),'-nL20.0-h',num2str(h),'.0-real.dat']);
		ui =  load(['sol_f',num2str(f),'-r',num2str(r),'-P',num2str(P),'-nL20.0-h',num2str(h),'.0-imag.dat']);
%		ur = ur/max(ur);
%		ui = ui/max(ui);
		u = ur + 1i *ui;
		u = u.';
		err_r(r,P) = norm(ur-ure)/norm(ure);
		err_i(r,P) = norm(ui-uie)/norm(uie);
	end
end

figure(1)
for i = 1 : rmax
semilogy(1:Pmax,err_r(i,:),strcat('-',Markers{i}))
hold on
end
grid
 legend('$r = 1$','$r = 2$', '$r = 3$','$r = 4$','$r = 5$','$r = 6$','$r = 7$','$r = 8$','Location','SouthWest','interpreter','latex','FontSize',16)
xlabel('P')
ylabel('E')
title(['Real part'],'interpreter','latex','FontSize',16)

figure(2)
for i = 1 : rmax
semilogy(1:Pmax,err_i(i,:),strcat('-',Markers{i}))
hold on
end
grid
 legend('$r = 1$','$r = 2$', '$r = 3$','$r = 4$','$r = 5$','$r = 6$','r = 7','Location','SouthWest','interpreter','latex','FontSize',16)
xlabel('P')
ylabel('E')
title(['Imaginary part'],'interpreter','latex','FontSize',16)


%return

figure
subplot(2,3,1)
scatter(p(:,1)/20,p(:,2)/20,15,abs(u),'filled'); 
title('Modulus (DAB)','FontSize',16,'interpreter','latex')
xlabel('$x~(\lambda_0)$','interpreter','latex','fontsize',12)
ylabel('$y~(\lambda_0)$','interpreter','latex','fontsize',12)
axis equal; 
axis image;
colormap jet
colorbar
subplot(2,3,2)
scatter(p(:,1)/20,p(:,2)/20,15,(ur),'filled'); 
title('Real part (DAB)','FontSize',16,'interpreter','latex')
xlabel('$x~(\lambda_0)$','interpreter','latex','fontsize',12)
ylabel('$y~(\lambda_0)$','interpreter','latex','fontsize',12)
axis equal; 
axis image;
colormap jet
colorbar
subplot(2,3,3)
scatter(p(:,1)/20,p(:,2)/20,15,(ui),'filled'); 
title('Imaginary part (DAB)','FontSize',16,'interpreter','latex')
xlabel('$x~(\lambda_0)$','interpreter','latex','fontsize',12)
ylabel('$y~(\lambda_0)$','interpreter','latex','fontsize',12)
axis equal; 
axis image;
colormap jet
colorbar

subplot(2,3,4)
scatter(p(:,1)/20,p(:,2)/20,15,abs(ue),'filled'); 
xlabel('$x~(\lambda_0)$','interpreter','latex','fontsize',12)
ylabel('$y~(\lambda_0)$','interpreter','latex','fontsize',12)
axis equal; 
axis image;
colormap jet
colorbar
title('Modulus (Exact)','FontSize',16,'interpreter','latex')

subplot(2,3,5)
scatter(p(:,1)/20,p(:,2)/20,15,ure,'filled');
xlabel('$x~(\lambda_0)$','interpreter','latex','fontsize',12)
ylabel('$y~(\lambda_0)$','interpreter','latex','fontsize',12)
axis equal;
axis image;
colormap jet
colorbar
title('Real Part (Exact)','FontSize',16,'interpreter','latex')

subplot(2,3,6)
scatter(p(:,1)/20,p(:,2)/20,15,uie,'filled');
xlabel('$x~(\lambda_0)$','interpreter','latex','fontsize',12)
ylabel('$y~(\lambda_0)$','interpreter','latex','fontsize',12)
axis equal;
axis image;
colormap jet
colorbar
title('Imaginary Part (Exact)','FontSize',16,'interpreter','latex')



shg


