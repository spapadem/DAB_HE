close all
clear all
xcr = load('solt_r1-P3-nL10.0-cr.dat');

P = 3;

x = zeros(size(xcr,2)/4,1);
y = zeros(size(xcr,2)/4,1);
Us = zeros(size(xcr,2)/4,P+1);
inds = 1 : size(xcr,2)/16;
figure
for t = 1 : size(xcr,1)
	a = 1;
	for i = 1 : 4 : size(xcr,2)-4
		x(a) = xcr(t,i);
		y(a) = xcr(t,i+1);
		Us(a,xcr(t,i+3)+1) = xcr(t,i+2);
		a = a + 1;
	end
%	[x,IX] = sort(x);
%	y = y(IX);
%	Us(:,1) = Us(IX,1);
%	Us(:,2) = Us(IX,2);
%	Us(:,3) = Us(IX,3);
%	Us(:,4) = Us(IX,4);
	title(['$T = $',num2str((t-1)*0.01),' s'],'FontSize',16)
	subplot(4,1,1)
	plot(x(inds    ),Us(inds    ,1),'b.')
	mx = max(abs([max(Us(inds,1)),min(Us(inds,1))]));
	mx = max(mx,eps);
	title('$u$','interpreter','latex','fontsize',14)
	axis([0 4 -mx mx])
	subplot(4,1,2)
	plot(x(inds+41 ),Us(inds+41 ,2),'b.')
	mx = max(abs([max(Us(inds+41,2)),min(Us(inds+41,2))]));
	mx = max(mx,eps);
	title('$\phi_1$','interpreter','latex','fontsize',14)
	axis([0 4 -mx mx])
	subplot(4,1,3)
	plot(x(inds+82 ),Us(inds+82 ,3),'b.')
	mx = max(abs([max(Us(inds+82,3)),min(Us(inds+82,3))]));
	mx = max(mx,eps);
	title('$\phi_2$','interpreter','latex','fontsize',14)
	axis([0 4 -mx mx])
	subplot(4,1,4)
	plot(x(inds+123),Us(inds+123,4),'b.')
	mx = max(abs([max(Us(inds+123,4)),min(Us(inds+123,4))]));
	mx = max(mx,eps);
	title('$\phi_3$','interpreter','latex','fontsize',14)
	axis([0 4 -mx mx])
	drawnow
	MovVec(t) = getframe(gcf);


end


fprintf('Creating video of the solution. \n')
videoname = ['1D_slice.avi'];
myVideo   = VideoWriter(videoname);
myVideo.FrameRate = 60;
myVideo.Quality   = 50;
open(myVideo);
writeVideo(myVideo,MovVec);
close(myVideo);


