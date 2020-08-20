clear all
fdir='../input_files/output/';

nfile=[10:90];
dep=load([fdir 'dep.out']);

for num=1:length(nfile)
clf
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);
etasrn=load([fdir 'etasrn_' fnum]);
mask=load([fdir 'mask_' fnum]);

plot(eta(50,:),'b')
hold on
plot(etasrn(50,:),'c--')
plot(-dep(50,:),'k--')
eta(mask>0)=NaN;
plot(eta(50,:),'ro')
axis([200 250 -0.5 1])
grid
pause
end
