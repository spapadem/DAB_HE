close all
clear all

re = load('errors_wrt_h_r_n.dat');

ndisc = 5;
hxinit = 0.5;
hx = hxinit./2.^((1:ndisc)-1);

loglog(hx,re(1,:),'bo-',hx,re(2,:),'r*-',hx,re(3,:),'gd-',hx,re(4,:),'ks-',hx,re(5,:),'mp-')
legend('$r = 1$','$r = 2$','$r = 3$','$r = 4$','$r = 5$',...
	   'interpreter','LaTeX','location','SouthEast','FontSize',12)
xlabel('$h$','interpreter','LaTeX','Fontsize',14)
ylabel('$\|u-u^h\|_2/\|u\|_2$','interpreter','LaTeX','Fontsize',14)
%title('$u(x) = t^2 \sin(12\pi x) \sin(12\pi y)$','interpreter','LaTeX','Fontsize',16)
grid minor



filename = 'errors_wrt_h_r';
saveas(gcf,[filename,'.fig'])
saveas(gcf,[filename,'.jpg'])
saveas(gcf,[filename,'.eps'],'epsc2')

%exit;
