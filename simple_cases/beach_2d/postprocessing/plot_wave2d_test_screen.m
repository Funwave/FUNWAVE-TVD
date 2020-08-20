clear all
fdir='../input_files/output/';

nfile=[20:90];
dep=load([fdir 'dep.out']);

for num=1:length(nfile)
clf
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);
%etasrn=load([fdir 'etasrn_' fnum]);
mask=load([fdir 'mask_' fnum]);

subplot(121)
%pcolor(etasrn-eta),shading flat
caxis([-0.2 0.2])
axis([200 250 1 500])

subplot(122)
eta(mask<1)=NaN;
pcolor(eta),shading flat
caxis([-1 1])
axis([200 250 1 500])
pause(0.1)
end
