clear all
fdir='../work/output/';


dep=load([fdir 'dep_00001']);

[n,m]=size(dep);
N=2*n-1;
M=m;

dx=1.0;
dy=1.0;
x=[0:M-1]*dx;
y=[0:N-1]*dy;



%nfile=[40 80 120];
nfile=[20];
min={'40'};

wid=8;
len=8;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
clf

ETA=zeros([N M]);
CH=zeros([N M]);
BEDS=zeros([N M]);
BEDB=zeros([N M]);
BB=zeros([N M]);
DEP=zeros([N M]);

DEP(1:n,:)=dep(:,:);
DEP(n+1:end,:)=dep(n-1:-1:1,:);

[ha, pos] = tight_subplot(3,1,[.05 0.5],[.1 .05],[.1 .1]) 
ax=[0 400 0 120];

for num=1:length(nfile)
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);
mask=load([fdir 'mask_' fnum]);
ch=load([fdir 'C_' fnum]);
beds=load([fdir 'DchgS_' fnum]);
bedb=load([fdir 'DchgB_' fnum]);


eta(mask<1)=NaN;
ch(mask<1)=NaN;


ETA(1:n,:)=eta(:,:);
ETA(n+1:end,:)=eta(n-1:-1:1,:);
CH(1:n,:)=ch(:,:);
CH(n+1:end,:)=ch(n-1:-1:1,:);

BEDS(1:n,:)=beds(:,:);
BEDS(n+1:end,:)=beds(n-1:-1:1,:);

BEDB(1:n,:)=bedb(:,:);
BEDB(n+1:end,:)=bedb(n-1:-1:1,:);

BB=BEDS+BEDB;


axes(ha(1));

pcolor(x,y,ETA),shading flat
caxis([-2 2])
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' \eta (m) ')
axis(ax)


axes(ha(2));

pcolor(x,y,log10(CH)),shading flat
caxis([0 0.5])

%title([' Time = ' min{num} ' sec '])
cbar=colorbar;
set(get(cbar,'ylabel'),'String',' log10(C) (g/L) ')
axis(ax)

ylabel(' y (m) ')


axes(ha(3));

contourf(x,y,BB,10,'Edgecolor','none')

cbar=colorbar;
set(get(cbar,'ylabel'),'String',' dZ (m) ')
axis(ax)
caxis([-0.002 0.002])

ylabel(' y (m) ')

xlabel(' x (m) ')

colormap(jet(12))

end


figure
[nn,mm]=size(BB);
for j=1:nn
    B(j)=mean(BB(j,:));
end
subplot(211)
plot(y,B,'LineWidth',2)
grid
xlabel('y (m)')
ylabel('Averaged Bed Change (m)')
axis([0 120 -0.004 0.002])
subplot(212)
plot(y,-DEP(:,1),'LineWidth',2)
hold on
plot([0 120],[0 0],'b--')
grid
xlabel('y (m)')
ylabel('Initial Depth(m)')
axis([0 120 -5.2 2])

