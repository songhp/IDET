% Phantom 

%%
clear all
thresh = 1e-4;
imgsize =128;
i=1;
for RadialNum=35:50
[PSNR, Count, t] = Phantom_recovery(imgsize, thresh, RadialNum);
snr(i,:) =PSNR;
iter(i,:) = Count;
time(i,:) =t;
i= i+1;
end
save(['Phantom', num2str(imgsize),'_thresh',num2str(thresh),'.mat']);

figure
% subplot(3,1,1)
plot(35:50, snr(:,1),'pc:', 'LineWidth',2,'MarkerSize',10)
hold on
plot(35:50, snr(:,2),'xb-', 'LineWidth',2,'MarkerSize',10)
plot(35:50, snr(:,3),'or-.', 'LineWidth',2,'MarkerSize',10)
plot(35:50, snr(:,4),'dk--', 'LineWidth',2,'MarkerSize',10)
legend('DORE', 'IDET-k', 'IDET-\beta', 'IDET-\gamma', 4)
ylabel('PSNR (dB)')
xlabel('Nunmbur of Radial Lines')
% xlim([0.1 0.401])
grid on
hold off
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
filename = ['Fig_phantom', num2str(imgsize),'PSNR_thresh',num2str(thresh),'.eps'];
print('-depsc2', filename);
% print -depsc2 Fig_phantom_PSNR.eps

figure
% subplot(3,1,2)
plot(35:50, iter(:,1),'pc:', 'LineWidth',2,'MarkerSize',10)
hold on
plot(35:50, iter(:,2),'xb-', 'LineWidth',2,'MarkerSize',10)
plot(35:50, iter(:,3),'or-.', 'LineWidth',2,'MarkerSize',10)
plot(35:50, iter(:,4),'dk--', 'LineWidth',2,'MarkerSize',10)
legend('DORE', 'IDET-k', 'IDET-\beta', 'IDET-\gamma', 1)
ylabel('Number of Iterations')
xlabel('Nunmbur of Radial Lines')
% xlim([0.1 0.401])
grid on
hold off
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
filename = ['Fig_phantom', num2str(imgsize),'Iter_thresh',num2str(thresh),'.eps'];
print('-depsc2', filename);
% print -depsc2 Fig_phantom_iterations.eps

figure
% subplot(3,1,3)
plot(35:50, time(:,1),'pc:', 'LineWidth',2,'MarkerSize',10)
hold on
plot(35:50, time(:,2),'xb-', 'LineWidth',2,'MarkerSize',10)
plot(35:50, time(:,3),'or-.', 'LineWidth',2,'MarkerSize',10)
plot(35:50, time(:,4),'dk--', 'LineWidth',2,'MarkerSize',10)
legend('DORE', 'IDET-k', 'IDET-\beta', 'IDET-\gamma', 1)
ylabel('CPU Time (s)')
xlabel('Nunmbur of Radial Lines')
% xlim([0.1 0.401])
grid on
hold off
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
filename = ['Fig_phantom', num2str(imgsize),'Time_thresh',num2str(thresh),'.eps'];
print('-depsc2', filename);
% print -depsc2 Fig_phantom_time.eps


%%
clear all
thresh = 1e-3;
imgsize =256;
i=1;
for RadialNum=35:50
[PSNR, Count, t] = Phantom_recovery(imgsize, thresh, RadialNum);
snr(i,:) =PSNR;
iter(i,:) = Count;
time(i,:) =t;
i= i+1;
end
save(['Phantom', num2str(imgsize),'_thresh',num2str(thresh),'.mat']);



figure
% subplot(3,1,1)
plot(35:50, snr(:,1),'pc:', 'LineWidth',2,'MarkerSize',10)
hold on
plot(35:50, snr(:,2),'xb-', 'LineWidth',2,'MarkerSize',10)
plot(35:50, snr(:,3),'or-.', 'LineWidth',2,'MarkerSize',10)
plot(35:50, snr(:,4),'dk--', 'LineWidth',2,'MarkerSize',10)
legend('DORE', 'IDET-k', 'IDET-\beta', 'IDET-\gamma', 4)
ylabel('PSNR (dB)')
xlabel('Nunmbur of Radial Lines')
% xlim([0.1 0.401])
grid on
hold off
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
filename = ['Fig_phantom', num2str(imgsize),'PSNR_thresh',num2str(thresh),'.eps'];
print('-depsc2', filename);
% print -depsc2 Fig_phantom_PSNR.eps

figure
% subplot(3,1,2)
plot(35:50, iter(:,1),'pc:', 'LineWidth',2,'MarkerSize',10)
hold on
plot(35:50, iter(:,2),'xb-', 'LineWidth',2,'MarkerSize',10)
plot(35:50, iter(:,3),'or-.', 'LineWidth',2,'MarkerSize',10)
plot(35:50, iter(:,4),'dk--', 'LineWidth',2,'MarkerSize',10)
legend('DORE', 'IDET-k', 'IDET-\beta', 'IDET-\gamma', 1)
ylabel('Number of Iterations')
xlabel('Nunmbur of Radial Lines')
% xlim([0.1 0.401])
grid on
hold off
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
filename = ['Fig_phantom', num2str(imgsize),'Iter_thresh',num2str(thresh),'.eps'];
print('-depsc2', filename);
% print -depsc2 Fig_phantom_iterations.eps

figure
% subplot(3,1,3)
plot(35:50, time(:,1),'pc:', 'LineWidth',2,'MarkerSize',10)
hold on
plot(35:50, time(:,2),'xb-', 'LineWidth',2,'MarkerSize',10)
plot(35:50, time(:,3),'or-.', 'LineWidth',2,'MarkerSize',10)
plot(35:50, time(:,4),'dk--', 'LineWidth',2,'MarkerSize',10)
legend('DORE', 'IDET-k', 'IDET-\beta', 'IDET-\gamma', 1)
ylabel('CPU Time (s)')
xlabel('Nunmbur of Radial Lines')
% xlim([0.1 0.401])
grid on
hold off
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
filename = ['Fig_phantom', num2str(imgsize),'Time_thresh',num2str(thresh),'.eps'];
print('-depsc2', filename);
% print -depsc2 Fig_phantom_time.eps


