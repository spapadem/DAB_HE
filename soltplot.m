close all
clear all
fprintf('Loading points \n')
p  = load('points_r1.dat');
ind = 1;
xI = 5;
for i = 1 :size(p,1)
	if p(i,1) <= xI
		pp(ind,:) = p(i,:);
		ind = ind + 1;
	end
end
p = pp;
fprintf('Loading solution \n')
ut = load('solt_r1-P5-nL10.0.dat');
ut1 = load('../Bounded_Domain/solt_r1-h0.1.dat');

%u = ut(1,:);
figure 


CreateVideo = 1;

r = 1;

for t = 1 :  size(ut,1)
	subplot(2,1,1)
	scatter(p(:,1),p(:,2),10,squeeze(ut(t,:)),'filled')
	title(['t: ',num2str(0.01*(t-1))])
	hold on
	plot([xI,xI],[0,3],'k-')
	hold off
    shading interp
    axis equal
    axis image
	cmax = 1*max(ut(t,:))/1;
	cmin = 1*min(ut(t,:))/1;
	clim = max(abs([cmax cmin]));
    caxis([-clim clim]);
    colorbar
    colormap jet
    drawnow

    subplot(2,1,2)
    scatter(p(:,1),p(:,2),10,squeeze(ut1(t,:)),'filled')
    title(['t: ',num2str(0.01*(t-1))])
    hold on
    plot([xI,xI],[0,3],'k-')
    hold off
    shading interp
    axis equal
    axis image
    cmax = 1*max(ut(t,:))/1;
    cmin = 1*min(ut(t,:))/1;
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
    videoname = ['Order_',num2str(r),'.avi'];
    myVideo   = VideoWriter(videoname);
    myVideo.FrameRate = 60;
    myVideo.Quality   = 50;
    open(myVideo);
    writeVideo(myVideo,MovVec);
    close(myVideo);
end

