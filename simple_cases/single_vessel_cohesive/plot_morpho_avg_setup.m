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



%nfile=[40 80 120];
nfile=[120];
mint={'80' '160' '240'};

wid=8;
len=4;
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

[ha, pos] = tight_subplot(2,1,[.05 0.5],[.1 .05],[.1 .1]) 
ax=[0 4500 0 120];

for num=1:length(nfile)
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);
mask=load([fdir 'mask_' fnum]);
ch=load([fdir 'C_' fnum]);
beds=load([fdir 'DchgS_' fnum]);

eta(mask<1)=NaN;
ch(mask<1)=NaN;


ETA(1:n,:)=eta(:,:);
ETA(n+1:end,:)=eta(n-1:-1:1,:);
CH(1:n,:)=ch(:,:);
CH(n+1:end,:)=ch(n-1:-1:1,:);

BEDS(1:n,:)=beds(:,:);
BEDS(n+1:end,:)=beds(n-1:-1:1,:);

BB=BEDS;

axes(ha(1));

pcolor(x,y,CH),shading flat
caxis([0 1.2])

cbar=colorbar;
set(get(cbar,'ylabel'),'String',' C (g/L) ')
axis(ax)

ylabel(' y (m) ')

axes(ha(2));

contourf(x,y,BEDS,10)

cbar=colorbar;
set(get(cbar,'ylabel'),'String',' dZ_{sus}  (m) ')
axis(ax)

caxis([-0.0005 0.0005])
ylabel(' y (m) ')

xlabel(' x (m) ')

colormap(jet(12))

end
print -djpeg100 wakes_cohesive_morpho.jpg
figure
[nn,mm]=size(BB);
for j=1:nn
    B(j)=mean(BB(j,:));
end

MaxH=max(ETA,[],2);
MinH=min(ETA,[],2);

subplot(211)
plot(y,B,'LineWidth',2)
grid
xlabel('y (m)')
ylabel('Averaged Bed Change (m)')
axis([0 120 -0.004 0.002])
subplot(212)
plot(y,-DEP(:,1),'LineWidth',2)
hold on
plot([15.7 102.32],[0 0],'b--','LineWidth',1.5)
plot(y,MaxH,'r--','LineWidth',1.5)
plot(y,MinH,'r--','LineWidth',1.5)
grid
xlabel('y (m)')
ylabel('Initial Depth(m)')
axis([0 120 -3.2 2])

print -djpeg100 section_mean_wave_morpho.jpg




