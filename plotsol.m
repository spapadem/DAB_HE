clear all
%close all

h = 2.0;
%h = 0.5;
NL = 40/h;

f = 24;

Pmax = 20;
rmax =  1;
err_r = zeros(5,Pmax);
err_i = zeros(5,Pmax);
for r =  1 : 2
	p = load(['points_r',num2str(r),'-h',num2str(h),'.0.dat']);
	fprintf('r = %2d\n',r)
	ure = load(['exact_f',num2str(f),'-r',num2str(r),'-h',num2str(h),'.0-real.dat']);
	uie = load(['exact_f',num2str(f),'-r',num2str(r),'-h',num2str(h),'.0-imag.dat']);
	ue = ure + 1i *uie;
	
	for P = 1: 2
		fprintf('Loading P = %2d of %2d\n',P,Pmax)
		ur =  load(['sol_f',num2str(f),'-r',num2str(r),'-PE',num2str(P),'-PN',num2str(P),'-nLE',num2str(NL),'.0-nLN',num2str(NL),'.0-h',num2str(h),'.0-real.dat']);
		ui =  load(['sol_f',num2str(f),'-r',num2str(r),'-PE',num2str(P),'-PN',num2str(P),'-nLE',num2str(NL),'.0-nLN',num2str(NL),'.0-h',num2str(h),'.0-imag.dat']);
%		ur = ur.';
%		ui = ui.';
		u = ur + 1i *ui;
%		u = u.';
		err_r(r,P) = norm(ur-ure)/norm(ure);
		err_i(r,P) = norm(ui-uie)/norm(uie);
	end
end

figure(1)
semilogy(1:Pmax,err_r(1,:),'o-',1:Pmax,err_r(2,:),'r*-',1:Pmax,err_r(3,:),'gd-',1:Pmax,err_r(4,:),'k--',1:Pmax,err_r(5,:),'mp-')
grid
 legend('$r = 1$','$r = 2$', '$r = 3$','$r = 4$','$r = 5$','Location','SouthEast','interpreter','latex','FontSize',16)
xlabel('P')
ylabel('E')
title(['$f = $',num2str(f),' Hz, Real part'],'interpreter','latex','FontSize',16)

figure(2)
semilogy(1:Pmax,err_i(1,:),'o-',1:Pmax,err_i(2,:),'r*-',1:Pmax,err_i(3,:),'gd-',1:Pmax,err_i(4,:),'k--',1:Pmax,err_i(5,:),'mp-')
grid
 legend('$r = 1$','$r = 2$', '$r = 3$','$r = 4$','$r = 5$','Location','SouthEast','interpreter','latex','FontSize',16)
xlabel('P')
ylabel('E')
title(['$f = $',num2str(f),' Hz, Imaginary part'],'interpreter','latex','FontSize',16)

%return


figure(3)
subplot(3,3,1)
scatter(p(:,1),p(:,2),15,abs(u),'filled'); 
title('Modulus (DAB)','FontSize',16,'interpreter','latex')
axis equal; 
axis image;
colormap jet
colorbar
subplot(3,3,2)
scatter(p(:,1),p(:,2),15,(ur),'filled'); 
title('Real part (DAB)','FontSize',16,'interpreter','latex')
axis equal; 
axis image;
colormap jet
colorbar
subplot(3,3,3)
scatter(p(:,1),p(:,2),15,(ui),'filled'); 
title('Imaginary part (DAB)','FontSize',16,'interpreter','latex')
axis equal; 
axis image;
colormap jet
colorbar

subplot(3,3,4)
scatter(p(:,1),p(:,2),15,abs(ue),'filled'); 
axis equal; 
axis image;
colormap jet
colorbar
title('Modulus (Exact)','FontSize',16,'interpreter','latex')
subplot(3,3,5)
scatter(p(:,1),p(:,2),15,ure,'filled');
axis equal;
axis image;
colormap jet
colorbar
title('Real Part (Exact)','FontSize',16,'interpreter','latex')
subplot(3,3,6)
scatter(p(:,1),p(:,2),15,uie,'filled');
axis equal;
axis image;
colormap jet
colorbar
title('Imaginary Part (Exact)','FontSize',16,'interpreter','latex')

subplot(3,3,7)
scatter(p(:,1),p(:,2),15,abs(ue-u),'filled');
axis equal; 
axis image;
colormap jet
colorbar
title('Modulus (Exact)','FontSize',16,'interpreter','latex')
subplot(3,3,8)
scatter(p(:,1),p(:,2),15,ure-ur,'filled');
axis equal;
axis image;
colormap jet
colorbar
title('Real Part (Exact)','FontSize',16,'interpreter','latex')
subplot(3,3,9)
scatter(p(:,1),p(:,2),15,uie-ui,'filled');
axis equal;
axis image;
colormap jet
colorbar
title('Imaginary Part (Exact)','FontSize',16,'interpreter','latex')



shg


