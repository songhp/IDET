load phase.mat

delta = linspace(0.05, 0.5, 16); 
figure;
plot(delta(1:end),phase1(1, 1:end),'bs:', 'LineWidth',3); % alg 
hold on
plot(delta(1:end),phase1(2, 1:end),'k:', 'LineWidth',4); % alg 

plot(delta(1:end),phase3(1, 1:end),'g--', 'LineWidth',3); % alg 
hold on
plot(delta(1:end),phase3(2, 1:end),'b--', 'LineWidth',3); % alg 
hold on
plot(delta(1:end),phase3(3, 1:end),'k--', 'LineWidth',4); % alg 

plot(delta(1:end),phase2(1, 1:end),'g-', 'LineWidth',3); % alg 
hold on
plot(delta(1:end),phase2(2, 1:end),'b-', 'LineWidth',3); % alg 
hold on
plot(delta(1:end),phase2(3, 1:end),'r-', 'LineWidth',4); % alg 


axis([0.05 0.55 0 0.6]);
ylabel('Sparsity k/m');
xlabel('Indeterminacy m/n');
grid on;
axis([0.05 0.5 0.05 0.6]);
legend('DORE', 'IDET-k',  'IDET-\gamma=0.7', 'IDET-\gamma=0.8', 'IDET-\gamma=0.9', ...
'IDET-\beta=0.7', 'IDET-\beta=0.8', 'IDET-\beta=0.9', 4);
% fn = strcat('Fig_PhaseTrans_DORE', '.fig'); 
% saveas(gcf, fn) 
set(gcf,'color','none');
set(gca,'color','none');
set(gcf,'InvertHardCopy','off');
print(gcf,'-depsc2',['Fig_PhaseTrans' '.eps']);



