function plotQCfailed(LP,k,cellID,qc)

figure('Position',[50 50 600 250]); set(gcf,'color','w');
hold on
plot(LP.V{1,k},'k','LineWidth',0.25)
xlabel('time-steps')
ylabel('voltage (mV)')
title(num2str(qc.logicVec))
axis tight
export_fig(['D:\genpath\',cellID,' ',int2str(k),' qc fail ',num2str(qc.logicVec)],'-pdf','-r100');
close