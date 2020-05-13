function [LP] = betweenSweepVmQC(LP,cellID)
%{
betweenSweepVmQC
%}
for k = 1:length(LP.V)                                                    % for each sweep
    rmp(1,k) = LP.stats{k,1}.qc.restVPre;
    rmp(2,k) = LP.stats{k,1}.qc.restVPost;
end

LP.rmp = rmp;

figure('Position',[50 50 300 250]); set(gcf,'color','w');
plot(LP.rmp')
xlabel('sweep #')
ylabel('resting membrane potential (mV)')
legend({'pre-stim','post-stim'})
box off
axis tight
ylim([-80 -40])

% save figure
export_fig(['D:\genpath\',cellID,' rmp'],'-pdf','-r100');
close