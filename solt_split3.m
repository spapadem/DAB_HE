close all
clear all
fprintf('Loading points \n')
p  = load('points_r3.dat');
ind = 1;
for i = 1 :size(p,1)
        if p(i,1) <= 9.7
                pp(ind,:) = p(i,:);
                ind = ind + 1;
        end
end
p = pp;

fprintf('Loading solution \n')
ut1 = load('../Bounded_Domain/solt_r1.dat');
%ut0  = load('../Bounded_Domain/solt_r3.dat');
%ut6  = load('solt_r3-P6.dat');
%ut61 = load('solt_P6.dat');

%u = ut(1,:);

	errt1 = zeros(size(ut1,1));
	utr   = load(['solt_r1-P',num2str(1),'-nL20.dat']);
	for tt = 1 : size(ut1,1)
	    errt1(tt) = norm(ut1(tt,:)-utr(tt,:))/(sqrt(9.7*3));
	end
%	E2(P+1) = norm(utr2(1:101,:)-ut1(1:101,:))/norm(ut1(1:101,:));

t = 0 : 0.05   : 10;
figure 
semilogy(t,errt1(:,1),'b-')
return
hold on
semilogy(t,errt1(:,2),'r--')
semilogy(t,errt1(:,3),'k.-')
semilogy(t,errt1(:,4),'m.')
semilogy(t,errt1(:,5),'g-')
axis([0 5 1e-7 1])
legend('P=0','P=1','P=2','P=3','P=4','Location','SouthEast')


figure
semilogy(0:7,E,'bo-',0:7,E10,'r*--')
%semilogy(0:7,E2,'kd-',0:7,E,'bo-',0:7,E10,'r*--')
axis([0 7 1e-3 1])
return


CreateVideo = 1;

r = 1;
errt = zeros(size(ut0,1),1);
errt1 = zeros(size(ut01,1),1);
for tt = 1 : size(ut0,1)

	errt(tt) = norm(ut0(tt,:)-ut6(tt,:))/(sqrt(9.7*3));
end

for tt = 1 : size(ut01,1)

    errt1(tt) = norm(ut01(tt,:)-ut61(tt,:))/(sqrt(9.7*3));
end

t  = 0 : 20/300 : 20;
semilogy(t1(1:end-1),errt1,'b-',t(1:end-1),errt,'r--')
legend('r = 1, P = 6','r = 3, P = 6','Location','SouthEast')
xlabel('t')
ylabel('e')
grid
%axis([0 10 1e-6 1e-1])
figure

subplot(2,1,1)
scatter(p(:,1),p(:,2),10,squeeze(ut0(end,:)),'filled')
axis equal
axis image
colormap jet
subplot(2,1,2)
scatter(p(:,1),p(:,2),10,squeeze(ut6(end,:)),'filled')
axis equal; axis image;
colormap jet


return


for t = 1 :  size(ut0,1)

	subplot(2,1,1)
	scatter(p(:,1),p(:,2),10,squeeze(ut0(t,:)),'filled')
	title(['t: ',num2str(0.01*(t-1))])
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
    title(['t: ',num2str(0.01*(t-1))])
    shading interp
    axis equal
    axis image
    cmax = 1*max(ut6(t,:))/1;
    cmin = 1*min(ut6(t,:))/1;
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

