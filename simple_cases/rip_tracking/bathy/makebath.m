% makebath.m
% 
% generate bathymetry for TVD-boussinesq model to study effect of symmetry
% breaking in periodic bathymetry and sediment
%
%--------------------------------------------------------------------------

% number of cells in periodic interval

  ncell=2;
  
% basic call spacing

  ycell=100;
  
% cross-shore scale 

  xcell=60;
  
% basic domain

  dx=1; dy=2; s=1/30;  xc=100;  a=1.5;  sig2=xc*xc/4; bet=2;
  x=(-11:dx:5*xcell); y=(dy/2:dy:ncell*ycell-dy/2);
  [xx,yy]=meshgrid(x,y);

  
  h=s*xx;
  
% initial cell centers

  yc=(0:ncell)*ycell;
  
% break symmetry

%  yc(ncell/2)=yc(ncell/2)-0.25*ycell;
  
  for index=1:ncell+1,
      
      h=h-a*exp( -(  ( xx-xc ).^2 + bet*(yy-yc(index)).^2)/sig2 );
      
  end
  
h(h>8)=8;

  v=[ 0 1 2 3 4 5 6 7 8 9 10 ];
  
  contour(h,v), axis equal
  
save -ASCII depth_a15.txt h
  
figure
wid=6;
len=6;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);
[n m]=size(h);
xx=[0:m-1]*dx;
yy=[0:n-1]*dy;
pcolor(xx,yy,-h),shading flat
colorbar
hold on
v1=[ -10:1:0 ];
contour(xx,yy,-h,v1,'w');
xlabel('x (m) ')
ylabel ('y (m)');



  