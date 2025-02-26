close all
clear all
%fprintf('Loading points \n')
p  = load('points_r1.dat');
ind = 1;
for i = 1 :size(p,1)
        if p(i,1) <= 500
                pp(ind,:) = p(i,:);
               ind = ind + 1;
        end
end
p = pp;

	E   = zeros(5,8);
ind = 1;
for r = [2 3 4 6 8]
	fprintf('r = %d\n',r)
	fprintf('Loading solution \n')
	ut1 = load(['../Bounded_Domain/solt_r',num2str(r),'.dat']);

	errt1 = zeros(size(ut1,1),5);
	for P = 1 : 8
		fprintf('Loading P = %d. \n',P)
		utr = load(['solt_r',num2str(r),'-P',num2str(P),'-nL20.0.dat']);
		E(ind,P) = norm(utr(1:61,:)-ut1(1:61,:))/norm(ut1(1:61,:));
	end
	ind = ind + 1;
end
%t = 0 : 0.1   : 10;
%figure 
%semilogy(t,errt1(:,1),'b-')
%hold on
%semilogy(t,errt1(:,2),'r--')
%semilogy(t,errt1(:,3),'k.-')
%semilogy(t,errt1(:,4),'m.')
%semilogy(t,errt1(:,5),'g-')
%axis([0 5 1e-7 1])
%legend('P=1','P=2','P=3','P=4','P=5','Location','SouthEast')

%E(2,4) = (E(2,3)+E(2,4)+E(2,5))/3;

figure
semilogy(1:8,E(1,:),'bo-')
hold on
semilogy(1:8,E(2,:),'mp-')
semilogy(1:8,E(3,:),'r*-')
semilogy(1:8,E(4,:),'gd-')
semilogy(1:8,E(5,:),'k--')
%semilogy(0:7,E2,'kd-',0:7,E,'bo-',0:7,E10,'r*--')
axis([1 8 1e-8 1e-2])
ylabel('$E$','interpreter','latex','FontSize',16)
xlabel('$P$','interpreter','latex','FontSize',16)
legend('$r = 2$','$r = 3$','$r = 4$','$r = 6$','$r = 8$','interpreter','latex','FontSize',16,'Location','SouthWest')

CreateVideo = 0;
return
errt = zeros(size(ut1,1),4);
ind = 1
for P = 8
	for r = [2 3 4]
    	fprintf('r = %d\n',r)
	    fprintf('Loading solution \n')
	    ut1 = load(['../Bounded_Domain/solt_r',num2str(r),'.dat']);

        fprintf('Loading P = %d. \n',P)
        utr = load(['solt_r',num2str(r),'-P',num2str(P),'-nL20.0.dat']);
		for tt = 1 : size(utr,1)
			errt(tt,ind) = norm(ut1(tt,:)-utr(tt,:))/(sqrt(5*3));
		end
			ind = ind + 1;
	end
	
	figure
	t  = 0 : 0.1 : 10;
	semilogy(t(1:end),errt(:,1),'b-',t(1:end),errt(:,2),'r--', t(1:end),errt(:,3),'g-',t,errt(:,4),'k.-')
	legend('$r = 2$','$r = 4$', '$r = 6$','$r = 8$','Location','SouthEast','interpreter','latex','FontSize',16)
	xlabel('t')
	ylabel('e')
	title(['P = ',num2str(P)])
	grid
	axis([0 10 1e-7 1e-1])
	drawnow
	shg
end


subplot(2,1,1)
scatter(p(:,1),p(:,2),10,squeeze(ut1(end,:)),'filled')
axis equal
axis image
colormap jet
subplot(2,1,2)
scatter(p(:,1),p(:,2),10,squeeze(utr(end,:)),'filled')
axis equal; axis image;
colormap jet


return


for t = 1 :  size(ut1,1)

	subplot(2,1,1)
	scatter(p(:,1),p(:,2),10,squeeze(ut1(t,:)),'filled')
	title(['t: ',num2str(0.01*(t-1))])
    shading interp
    axis equal
    axis image
	cmax = 1*max(ut1(t,:))/1;
	cmin = 1*min(ut1(t,:))/1;
	clim = max(abs([cmax cmin]));
    caxis([-clim clim]);
    colorbar
    colormap jet
    
    subplot(2,1,2)
    scatter(p(:,1),p(:,2),10,squeeze(utr(t,:)),'filled')
    title(['t: ',num2str(0.01*(t-1))])
    shading interp
    axis equal
    axis image
    cmax = 1*max(utr(t,:))/1;
    cmin = 1*min(utr(t,:))/1;
    clim = max(abs([cmax cmin]));
    caxis([-clim clim]);
    colorbar
    colormap jet



	drawnow
    if CreateVideo
        MovVec(t) = getframe(gcf);
    end
end



% Create a video of the solution.
if CreateVideo
    fprintf('Creating video of the solution. \n')
    videoname = ['Order_',num2str(r),'_split.avi'];
    myVideo   = VideoWriter(videoname);
    myVideo.FrameRate = 60;
    myVideo.Quality   = 50;
    open(myVideo);
    writeVideo(myVideo,MovVec);
    close(myVideo);
end

