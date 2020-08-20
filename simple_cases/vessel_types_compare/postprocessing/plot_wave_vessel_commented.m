% Clearing Matlab workspace
clear all

% -----------------------
% ----- User Input ------
% -----------------------

% Directory of output data files
fdir='output/';

% Time series files to plot
nfile=[20 40];

% Time values for series in minutes
times={'20' '40'};



% Determining domain dimensions
eta=load([fdir 'eta_00001']);
[n,m]=size(eta);

% Setting up partition
dx=1.0;
dy=1.0;
x=[0:m-1]*dx;
y=[0:n-1]*dy;

% Dimensions of plot window 
wid=8;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf


for num=1:length(nfile)

% Padding integer values with zeros
% to be 5 letters long e.g. 1 -> 00001
fnum=sprintf('%.5d',nfile(num));

% Loading wave displacement from file
eta=load([fdir 'eta_' fnum]);

% Plotting data on different subplot regions
subplot(length(nfile),1, num)

% Plotting wave displacement
pcolor(x,y,eta),shading flat

% Colorbar range
caxis([-1.5 1.5])

title([' Time = ' times{num} ' sec '])

ylabel(' y (m) ')
xlabel(' x (m) ')

set(gcf,'Renderer','zbuffer')

end
%print -djpeg eta_inlet_shoal_irr.jpg