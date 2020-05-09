function plotSuprathreshold(LP,sp,k,cellID)

%{
plotSuprathreshold
%}

figure('Position',[50 50 1500 250]); set(gcf,'color','w');
subplot(1,4,1:3)
hold on
plot(LP.V{1,k},'k','LineWidth',0.25)
plot(sp.peakTime,LP.V{1,k}(sp.peakTime),'.r','markersize',16)
plot(sp.maxdVdtTime,LP.V{1,k}(sp.maxdVdtTime),'.c','markersize',16)        % plot max dV/dt
plot(sp.thresholdTime,LP.V{1,k}(sp.thresholdTime),'.g','markersize',16)                % threshold
plot(sp.thresholdRefTime,LP.V{1,k}(sp.thresholdRefTime),'.g','markersize',10)          % refined threshold
plot(sp.troughTime,LP.V{1,k}(sp.troughTime),'.b','markersize',16)                         % trough
plot(sp.halfHeightTimeUpPT,LP.V{1,k}(sp.halfHeightTimeUpPT),'.k','markersize',16)      % half height time up
plot(sp.halfHeightTimeDownPT,LP.V{1,k}(sp.halfHeightTimeDownPT),'.k','markersize',16)  % half height time down
plot(sp.halfHeightTimeUpTP,LP.V{1,k}(sp.halfHeightTimeUpTP),'.y','markersize',16)      % half height time up
plot(sp.halfHeightTimeDownTP,LP.V{1,k}(sp.halfHeightTimeDownTP),'.y','markersize',16)  % half height time down
xlabel('time-steps')
ylabel('voltage (mV)')
axis tight

subplot(1,4,4)
hold on
plot(LP.V{1,k},'k','LineWidth',0.25)
plot(sp.peakTime,LP.V{1,k}(sp.peakTime),'.r','markersize',16)
plot(sp.maxdVdtTime,LP.V{1,k}(sp.maxdVdtTime),'.c','markersize',16)                       % plot max dV/dt
plot(sp.thresholdTime,LP.V{1,k}(sp.thresholdTime),'.g','markersize',16)                   % threshold
plot(sp.thresholdRefTime,LP.V{1,k}(sp.thresholdRefTime),'.g','markersize',10)             % refined threshold
plot(sp.troughTime,LP.V{1,k}(sp.troughTime),'.b','markersize',16)                         % trough
plot(sp.halfHeightTimeUpPT,LP.V{1,k}(sp.halfHeightTimeUpPT),'.k','markersize',16)         % half height time up
plot(sp.halfHeightTimeDownPT,LP.V{1,k}(sp.halfHeightTimeDownPT),'.k','markersize',16)     % half height time down
plot(sp.halfHeightTimeUpTP,LP.V{1,k}(sp.halfHeightTimeUpTP),'.y','markersize',16)         % half height time up
plot(sp.halfHeightTimeDownTP,LP.V{1,k}(sp.halfHeightTimeDownTP),'.y','markersize',16)     % half height time down
xlabel('time-steps')
ylabel('voltage (mV)')
axis tight
xlim([sp.peakTime(1)-(2/LP.acquireRes) sp.troughTime(1)+(2/LP.acquireRes)])

% save figure
export_fig(['D:\genpath\',cellID,' ',int2str(k),' suprathreshold parameters'],'-pdf','-r100');
close