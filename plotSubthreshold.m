function plotSubthreshold(LP,k,cellID)


figure('Position',[50 50 1500 250]); set(gcf,'color','w');
hold on

plot(LP.V{1,k})

xlabel('time-steps')
ylabel('voltage (mV)')
axis tight
box off


hold on
xlabel('time-steps')
ylabel('voltage (mV)')
title('expon fit to minimum')
axis tight
box off