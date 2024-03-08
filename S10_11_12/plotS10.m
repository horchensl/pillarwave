
index = 5;

figure(1)
clf
imagesc(time,depth*1e3,raw(:,:,index)./max(raw(:,:,index),[],'all')) % plot 
colormap gray
colorbar
caxis([-0.5 0.5])
ylim([11 26])
xlabel('time [s]')
ylabel('depth [mm]')