clear all
%close all

h = 2;


f = 73;

Pmax = 5;

rmax = 5;

err_r = zeros(5,Pmax);
err_i = zeros(5,Pmax);
for r = 2
 	p = load(['points_r',num2str(r),'-h',num2str(h),'.0.dat']);
	fprintf('r = %2d\n',r)
	ure = load(['exact_f',num2str(f),'-r',num2str(r),'-h',num2str(h),'.0-real.dat']);
	uie = load(['exact_f',num2str(f),'-r',num2str(r),'-h',num2str(h),'.0-imag.dat']);
%	ure = ure/max(ure);
%	uie = uie/max(uie);
	ue = ure + 1i *uie;
	ue = ue.';
	
	for P =   Pmax
		fprintf('Loading P = %2d of %2d\n',P,Pmax)
		ur =  load(['sol_f',num2str(f),'-r',num2str(r),'-P',num2str(P),'-nL20.0-h',num2str(h),'.0-real.dat']);
		ui = -load(['sol_f',num2str(f),'-r',num2str(r),'-P',num2str(P),'-nL20.0-h',num2str(h),'.0-imag.dat']);
%		u = u.';
%		ur = ur/max(ur);
%		ui = ui/max(ui);
%		fprintf('P = %2d , Real ratio = %6f, Imag. ratio = %6f \n',P,mean(ur./ure),mean(ui./uie))
		u = ur + 1i *ui;
		err_r(r,P) = norm(ur-ure)/norm(ure);
		err_i(r,P) = norm(ui-uie)/norm(uie);
	end
end
%err_r(4,:) = 0;
%err_i(4,:) = 0;
%err_r(5,:) = 0;
%err_i(5,:) = 0;
figure(1)
semilogy(1:Pmax,err_r(1,:),'o-',1:Pmax,err_r(2,:),'r*-',1:Pmax,err_r(3,:),'gd-',1:Pmax,err_r(4,:),'k--',1:Pmax,err_r(5,:),'mp-')
grid
 legend('$r = 1$','$r = 2$', '$r = 3$','$r = 4$','$r=5$','Location','SouthWest','interpreter','latex','FontSize',16)
xlabel('P')
ylabel('E')
title(['$f = $',num2str(f),' Hz, Real part'],'interpreter','latex','FontSize',16)

filename = ['error_f',num2str(f),'_h_',num2str(h),'_nL_',num2str(20),'_real'];
%saveas(gcf,[filename,'.fig'])
%saveas(gcf,[filename,'.jpg'])
%saveas(gcf,[filename,'.eps'],'psc2')



figure(2)
semilogy(1:Pmax,err_i(1,:),'o-',1:Pmax,err_i(2,:),'r*-',1:Pmax,err_i(3,:),'gd-',1:Pmax,err_i(4,:),'k--',1:Pmax,err_i(5,:),'mp-')
grid
 legend('$r = 1$','$r = 2$', '$r = 3$','$r = 4$','$r=5$','Location','SouthWest','interpreter','latex','FontSize',16)
xlabel('P')
ylabel('E')
title(['$f = $',num2str(f),' Hz, Imaginary part'],'interpreter','latex','FontSize',16)

filename = ['error_f',num2str(f),'_h_',num2str(h),'_nL_',num2str(20),'_imag'];
%saveas(gcf,[filename,'.fig'])
%saveas(gcf,[filename,'.jpg'])
%saveas(gcf,[filename,'.eps'],'psc2')

%return

figure(3)
subplot(2,3,1)
scatter(p(:,1),p(:,2),15,abs(u),'filled'); 
title('Modulus (DAB)','FontSize',16,'interpreter','latex')
axis equal; 
axis image;
colormap jet
colorbar
subplot(2,3,2)
scatter(p(:,1),p(:,2),15,(ur),'filled'); 
title('Real part (DAB)','FontSize',16,'interpreter','latex')
axis equal; 
axis image;
colormap jet
colorbar
subplot(2,3,3)
scatter(p(:,1),p(:,2),15,(ui),'filled'); 
title('Imaginary part (DAB)','FontSize',16,'interpreter','latex')
axis equal; 
axis image;
colormap jet
colorbar

subplot(2,3,4)
scatter(p(:,1),p(:,2),15,abs(ue),'filled'); 
axis equal; 
axis image;
colormap jet
colorbar
title('Modulus (Exact)','FontSize',16,'interpreter','latex')

subplot(2,3,5)
scatter(p(:,1),p(:,2),15,ure,'filled');
axis equal;
axis image;
colormap jet
colorbar
title('Real Part (Exact)','FontSize',16,'interpreter','latex')

subplot(2,3,6)
scatter(p(:,1),p(:,2),15,uie,'filled');
axis equal;
axis image;
colormap jet
colorbar
title('Imaginary Part (Exact)','FontSize',16,'interpreter','latex')



shg


