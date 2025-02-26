close all
%clear all
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
ut0 = load('solt_r1-h0.05.dat');
ut6   = load(['solt_r1-P',num2str(P),'-nL20.0.dat'])/0.05;

CreateVideo = 0;


    errt = zeros(size(ut0,1));
    for tt = 1 : size(ut0,1)
        errt(tt) = norm(ut0(tt,:)-ut6(tt,:))/(sqrt(9.7*3));
    end


for t = [101]
	figure
	subplot(2,1,1)
	scatter(p(:,1),p(:,2),10,squeeze(ut0(t,:)),'filled')
	title(['$t = $ ',num2str(0.1*(t-1)),' s'],'Fontsize',16,'interpreter','latex')
    shading interp
    axis equal
    axis image
	cmax = 1*max(ut0(t,:))/1;
	cmin = 1*min(ut0(t,:))/1;
	clim = max(abs([cmax cmin]));
    caxis([-clim clim]);
    colorbar
    colormap jet
    
    subplot(2,1,2)
    scatter(p(:,1),p(:,2),10,squeeze(ut6(t,:)),'filled')
%	title(['$t:$ ',num2str(0.01*(t-1))],'Fontsize',16,'interpreter','latex')
    shading interp
    axis equal
    axis image
    cmax = 1*max(ut6(t,:))/1;
    cmin = 1*min(ut6(t,:))/1;
    clim = max(abs([cmax cmin]));
    caxis([-clim clim]);
    colorbar
    colormap jet

	filename = ['solt_r4_t',num2str(0.1*(t-1))];
	saveas(gcf,[filename,'.fig'])
	saveas(gcf,[filename,'.eps'],'psc2')



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

