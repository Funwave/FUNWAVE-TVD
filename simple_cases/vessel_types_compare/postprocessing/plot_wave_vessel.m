clear all
fdir='../output/';

eta=load([fdir 'eta_00001']);

[n,m]=size(eta);
dx=1.0;
dy=1.0;
x=[0:m-1]*dx;
y=[0:n-1]*dy;

nfile=[5 10];
min={'5' '10'};

wid=5;
len=8;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf


for num=1:length(nfile)
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);

subplot(length(nfile),1, num)

pcolor(x,y,eta),shading flat
hold on
caxis([-2.5 2.5])
title([' Time = ' min{num} ' sec '])

ylabel(' y (m) ')

xlabel(' x (m) ')


set(gcf,'Renderer','zbuffer')

end
%print -djpeg eta_inlet_shoal_irr.jpg