function plotBadSweep(LP,k,cellID)

figure('Position',[50 50 600 250]); set(gcf,'color','w');
hold on
plot(LP.V{1,k},'k','LineWidth',0.25)
xlabel('time-steps')
ylabel('voltage (mV)')
axis tight
export_fig(['D:\genpath\',cellID,' bad sweep ',int2str(k)],'-pdf','-r100');
close