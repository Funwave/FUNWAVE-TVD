clear all
fdir='output/';


dep=load([fdir 'dep_00000']);

[n,m]=size(dep);
N=2*n-1;
M=m;

dx=1.0;
dy=1.0;
x=[0:M-1]*dx;
y=[0:N-1]*dy;



nfile=[40 80 120];

min={'200' '400' '600'};

wid=8;
len=8;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf

ETA=zeros([N M]);
CH=zeros([N M]);

[ha, pos] = tight_subplot(6,1,[.05 0.5],[.1 .05],[.1 .1]) 
ax=[0 4500 0 120];

for num=1:length(nfile)
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);
mask=load([fdir 'mask_' fnum]);
ch=load([fdir 'C_' fnum]);

eta(mask<1)=NaN;
ch(mask<1)=NaN;


ETA(1:n,:)=eta(:,:);
ETA(n+1:end,:)=eta(n-1:-1:1,:);
CH(1:n,:)=ch(:,:);
CH(n+1:end,:)=ch(n-1:-1:1,:);

axes(ha(2*(num-1)+1));

pcolor(x,y,ETA),shading flat
hold on
caxis([-0.3 1.0])
title([' Time = ' min{num} ' sec '])
axis(ax)

cbar=colorbar;
set(get(cbar,'ylabel'),'String',' \eta (m) ')

ylabel(' y (m) ')

axes(ha(2*(num-1)+2));

pcolor(x,y,CH),shading flat
caxis([0 1.2])

%title([' Time = ' min{num} ' sec '])
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' C (g/L) ')
axis(ax)

if num==length(nfile)
xlabel(' x (m) ')
end

ylabel(' y (m) ')


end
print -djpeg100 wakes_cohesive_sed.jpg