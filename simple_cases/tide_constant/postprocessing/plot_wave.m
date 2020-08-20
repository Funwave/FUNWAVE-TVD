clear all
fdir='../work/output/';

eta=load([fdir 'eta_00000']);
[n,m]=size(eta);
dx=2.0;

x=[0:m-1]*dx;

% wavemaker and sponge
tdw=30*dx;
wd=10.0;
wc=150;
x_wm=[wc-wd wc+wd];
xw_sponge=[0 tdw];
xe_sponge=[x(end)-tdw x(end)];

%ns=1;
%ne=199;

ns=100;
ne=ns;

wid=8;
len=4;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf


% Set up file and options for creating the movie
vidObj = VideoWriter('movie.avi');  % Set filename to write video file
vidObj.FrameRate=10;  % Define the playback framerate [frames/sec]
open(vidObj);

for num=ns:ne

       
fnum=sprintf('%.5d',num);
eta=load([fdir 'eta_' fnum]);

clf
plot(x,eta(1,:))
hold on
plot([x_wm(1) x_wm(1)],[-10 10],'r')
plot([x_wm(2) x_wm(2)],[-10 10],'r')
plot([xw_sponge(1) xw_sponge(1)],[-10 10],'r--')
plot([xw_sponge(2) xw_sponge(2)],[-10 10],'r--')
plot([xe_sponge(1) xe_sponge(1)],[-10 10],'r--')
plot([xe_sponge(2) xe_sponge(2)],[-10 10],'r--')

%axis([0 1024 -1 1])
grid
xlabel('x(m)')
ylabel('eta(m)')
pause(0.1)

 currframe=getframe(gcf);
    writeVideo(vidObj,currframe);  % Get each recorded frame and write it to filename defined above

end

close(vidObj)

