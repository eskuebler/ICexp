function supraStats = processSuprathresholdLongPulse(LP,cellID,params,k)

%{
processSuprathresholdLongPulse
%}

% get spike times based on voltage trace
LP.putSpTimes = LP.stimOn(1,k)+...
    find(LP.V{1,k}(LP.stimOn(1,k):LP.stimOff(1,k))>=params.thresholdV)-1;   % 0 mV threshold
% 
% % spikes based on dV/dt (smooth V with Bessel filter ?)
% dVdt = diff(LP.V{1,k}(1:end))/LP.acquireRes;
% tempSP = find(dVdt > (20/LP.acquireRes));
% c = 1;
% for j = 1:length(tempSP)-1
%     if sum(dVdt(1,tempSP(j):putSP.dVdttempSP(j+1)) < 0) > 0
%         putSP.dVdt(c) = putSP.dVdt(j);
%         c = c + 1;
%     end
% end
% clear c j
% >20mV/ms and returns below 0mV/ms between events

if ~isempty(LP.putSpTimes)
    
    [int4Peak,LP.putSpTimes2] = int4APs(LP.putSpTimes);
    
    if ~isempty(LP.putSpTimes2)
        [peak,peakTime] = estimatePeak(LP,int4Peak,k);                       % estimate of peak
        [dVdtNthreshold,peak,peakTime,LP] = estimateMaxdVdtNthreshold ...
            (LP,peak,peakTime,k,params);
        [peak,peakTime,dVdtNthreshold,LP] = QCpeakNthreshold(LP,...
            peak,peakTime,dVdtNthreshold,params);
        [dVdtNthresholdR] = refineThreshold(LP,dVdtNthreshold,peakTime,k);
        [peak,peakTime,dVdtNthreshold,dVdtNthresholdR] = QCpeakNthreshold2(LP.acquireRes,...
            peak,peakTime,dVdtNthreshold,dVdtNthresholdR,params);
        [trough,troughTime] = estimateTrough(LP,peakTime,k);
        [peak,peakTime,trough,troughTime,dVdtNthreshold,dVdtNthresholdR] = ...
            QCpeakNtrough(peak,peakTime,trough,troughTime,params,dVdtNthreshold,dVdtNthresholdR);
        [spParams,peak,peakTime,trough,troughTime,dVdtNthreshold,dVdtNthresholdR,LP] = ...
            estimateSpParams(LP,peak,peakTime,trough,troughTime,dVdtNthreshold,dVdtNthresholdR,k,params);
        [spTrainParams] = estimateAPTrainParams(LP,dVdtNthresholdR,spParams,k);
        % estimate plateau potential
        % assessPersistence

        % waveforms
        params.windowBeg = round(3/LP.acquireRes);                                 % time prior to action potential
        params.windowEnd = round(10/LP.acquireRes);                                % time after action potential
        params.windowBegS = round(2/params.sampleRTdt);                              % time prior to action potential
        params.windowEndS = round(6/params.sampleRTdt);                              % time after action potential
        wf = getWaveforms(LP,params,peakTime,k);
        
        supraStats.spTimes = LP.putSpTimes2;
        supraStats.spWaveforms = wf;
        supraStats.peak = peak; 
        supraStats.peakTime = peakTime;
        supraStats.maxdVdt = dVdtNthreshold.maxdVdt;
        supraStats.maxdVdtTime = dVdtNthreshold.maxdVdtTime;
        supraStats.threshold = dVdtNthreshold.threshold;
        supraStats.thresholdTime = dVdtNthreshold.thresholdTime;
        supraStats.thresholdRef = dVdtNthresholdR.thresholdRef;
        supraStats.thresholdRefTime = dVdtNthresholdR.thresholdRefTime;
        supraStats.trough = trough;
        supraStats.troughTime = troughTime;
        supraStats.peak2trough = (troughTime-peakTime).*LP.acquireRes;
        supraStats.heightPT = spParams.heightPT;
        supraStats.halfHeightTimeUpPT = spParams.halfHeightTimeUpPT;
        supraStats.halfHeightTimeDownPT = spParams.halfHeightTimeDownPT;
        supraStats.fullWidthPT = spParams.fullWidthPT;
        supraStats.heightTP = spParams.heightTP;
        supraStats.halfHeightTimeUpTP = spParams.halfHeightTimeUpTP;
        supraStats.halfHeightTimeDownTP = spParams.halfHeightTimeDownTP;
        supraStats.fullWidthTP = spParams.fullWidthTP;
        supraStats.peakUpStroke = spParams.peakUpStroke;
        supraStats.peakDownStroke = spParams.peakDownStroke;
        supraStats.peakStrokeRatio = spParams.peakStrokeRatio;
        supraStats.latency = spTrainParams.latency;
        supraStats.meanFR50 = spTrainParams.meanFR50;
        supraStats.meanFR100 = spTrainParams.meanFR100;
        supraStats.meanFR250 = spTrainParams.meanFR250;
        supraStats.meanFR500 = spTrainParams.meanFR500;
        supraStats.meanFR750 = spTrainParams.meanFR750;
        supraStats.meanFR1000 = spTrainParams.meanFR1000;
        supraStats.peakAdapt = spTrainParams.peakAdapt;
        supraStats.ISI = spTrainParams.ISI;
        supraStats.instaRate = spTrainParams.instaRate;
        supraStats.meanISI = spTrainParams.meanISI;
        supraStats.cvISI = spTrainParams.cvISI;
        supraStats.adaptIndex = spTrainParams.adaptIndex;
        supraStats.adaptIndex2 = spTrainParams.adaptIndex2;
        supraStats.peakAdapt2 = spTrainParams.peakAdapt2;
        supraStats.delay = spTrainParams.delay;
        supraStats.burst = spTrainParams.burst;
        if exist('fast_trough')
            supraStats.fastTrough = spParams.fast_trough;
            supraStats.fastTroughDur = spParams.fast_trough_dur;
        end
        if exist('slow_trough')
            supraStats.slowTrough = spParams.slow_trough;
            supraStats.slowTroughDur = spParams.slow_trough_dur;
        end
        supraStats.sweepID = k;
        
        figure('Position',[50 50 1500 250]); set(gcf,'color','w');
        subplot(1,4,1:3)
        hold on
        plot(LP.V{1,k},'k','LineWidth',0.25)
        plot(peakTime,LP.V{1,k}(peakTime),'.r','markersize',16)
        plot(dVdtNthreshold.maxdVdtTime,...
            LP.V{1,k}(dVdtNthreshold.maxdVdtTime),'.c','markersize',16)                       % plot max dV/dt
        plot(dVdtNthreshold.thresholdTime,...
            LP.V{1,k}(dVdtNthreshold.thresholdTime),'.g','markersize',16)                   % threshold
        plot(dVdtNthresholdR.thresholdRefTime,...
            LP.V{1,k}(dVdtNthresholdR.thresholdRefTime),'.g','markersize',10)             % refined threshold
        plot(troughTime,LP.V{1,k}(troughTime),'.b','markersize',16)                         % trough
        plot(spParams.halfHeightTimeUpPT,...
            LP.V{1,k}(spParams.halfHeightTimeUpPT),'.k','markersize',16)         % half height time up
        plot(spParams.halfHeightTimeDownPT,...
            LP.V{1,k}(spParams.halfHeightTimeDownPT),'.k','markersize',16)     % half height time down
        plot(spParams.halfHeightTimeUpTP,...
            LP.V{1,k}(spParams.halfHeightTimeUpTP),'.y','markersize',16)         % half height time up
        plot(spParams.halfHeightTimeDownTP,...
            LP.V{1,k}(spParams.halfHeightTimeDownTP),'.y','markersize',16)     % half height time down
        xlabel('time-steps')
        ylabel('voltage (mV)')
        axis tight
        subplot(1,4,4)
        hold on
        plot(LP.V{1,k},'k','LineWidth',0.25)
        plot(peakTime,LP.V{1,k}(peakTime),'.r','markersize',16)
        plot(dVdtNthreshold.maxdVdtTime,...
            LP.V{1,k}(dVdtNthreshold.maxdVdtTime),'.c','markersize',16)                       % plot max dV/dt
        plot(dVdtNthreshold.thresholdTime,...
            LP.V{1,k}(dVdtNthreshold.thresholdTime),'.g','markersize',16)                   % threshold
        plot(dVdtNthresholdR.thresholdRefTime,...
            LP.V{1,k}(dVdtNthresholdR.thresholdRefTime),'.g','markersize',10)             % refined threshold
        plot(troughTime,LP.V{1,k}(troughTime),'.b','markersize',16)                         % trough
        plot(spParams.halfHeightTimeUpPT,...
            LP.V{1,k}(spParams.halfHeightTimeUpPT),'.k','markersize',16)         % half height time up
        plot(spParams.halfHeightTimeDownPT,...
            LP.V{1,k}(spParams.halfHeightTimeDownPT),'.k','markersize',16)     % half height time down
        plot(spParams.halfHeightTimeUpTP,...
            LP.V{1,k}(spParams.halfHeightTimeUpTP),'.y','markersize',16)         % half height time up
        plot(spParams.halfHeightTimeDownTP,...
            LP.V{1,k}(spParams.halfHeightTimeDownTP),'.y','markersize',16)     % half height time down
        xlabel('time-steps')
        ylabel('voltage (mV)')
        axis tight
        xlim([peakTime(1)-(2/LP.acquireRes) troughTime(1)+(2/LP.acquireRes)])
        
        % save figure
        export_fig(['D:\genpath\',cellID,' suprathreshold parameters ',int2str(k)],'-pdf','-r100');
        close
                
        clear int4Peak PutAPTimes peak peakTime maxdVdt maxdVdtTime threshold thresholdTime thresholdRef trough troughTime ...
            heightPT halfHeightTimeUpPT halfHeightTimeDownPT fullWidthPT heightTP halfHeightTimeUpTP halfHeightTimeDownTP fullWidthTP ...
            peakUpStroke peakDownStroke peakStrokeRatio fast_trough slow_trough int4Peak2 rheobaseHeight bestRest addIndOut ...
            fast_trough_dur slow_trough_dur
    end
%     % response processing
%     spike_times = s.APtimes{k,N.goodSweeps};
%     raster = zeros(1,stimOff-stimOn);
%     raster(1,spike_times) = 1;
%     rasterDs = zeros(1,2200);
%     downsampleBinary
%     rasterConv = conv(rasterDs,params.gaussF,'same');
%     j = length(rasterConv);
%     rasterConv = [rasterConv,zeros(1,params.rasterRes)];
%     rasterPSTH = zeros(1,2200+100);
%     for t = params.rasterRes+1:j
%         rasterPSTH(t) = sum(rasterConv(t-params.rasterRes:t+params.rasterRes));
%     end
%     clear t j
%     % y = [rasterPSTH,zeros(1,100)];
%     %                                 y = smoothdata(rasterPSTH,'gaussian',100);
%     % Y = [Y; y(900:2100)];
% 
%     % plot result
%     if params.plotOpt == 1
%         figure('Position', [20,20,400,800]);
%         subplot(5,1,1)
%         plot(raster,'k')
%         title('spike raster')
%         axis tight;
%         box off;
%         set(gca,'LineWidth',2);
%     %                         xlim([800 2200])
%         subplot(5,1,2)
%         plot(rasterDs,'k')
%         title('down-sampled raster')
%         axis tight;
%     %                         xlim([800 2200])
%         box off;
%         set(gca,'LineWidth',2);
%         subplot(5,1,3)
%         plot(rasterConv,'k')
%         title('20 ms Gaussian')
%         axis tight;
%     %                         xlim([800 2200])
%         box off;
%         set(gca,'LineWidth',2);
%         subplot(5,1,4)
%         plot(rasterPSTH,'k')
%         title('75 ms window')
%         axis tight;
%     %                         xlim([800 2200])
%         box off;
%         set(gca,'LineWidth',2);
%         subplot(5,1,5)
%         plot(smoothdata(rasterPSTH,'gaussian',100),'k')
%         title('smoothdata function 100 ms')
%         xlabel('time (ms)')
%         axis tight;
%     %                         xlim([800 2200])
%         box off;
%         set(gca,'LineWidth',2);
%         drawnow
%         pause(1)
%         close
else
    supraStats.spTimes = NaN;
    supraStats.spWaveforms = NaN;
    supraStats.peak = NaN; 
    supraStats.peakTime = NaN;
    supraStats.maxdVdt = NaN;
    supraStats.maxdVdtTime = NaN;
    supraStats.threshold = NaN;
    supraStats.thresholdTime = NaN;
    supraStats.thresholdRef = NaN;
    supraStats.thresholdRefTime = NaN;
    supraStats.trough = NaN;
    supraStats.troughTime = NaN;
    supraStats.peak2trough = NaN;
    supraStats.heightPT = NaN;
    supraStats.halfHeightTimeUpPT = NaN;
    supraStats.halfHeightTimeDownPT = NaN;
    supraStats.fullWidthPT = NaN;
    supraStats.heightTP = NaN;
    supraStats.halfHeightTimeUpTP = NaN;
    supraStats.halfHeightTimeDownTP = NaN;
    supraStats.fullWidthTP = NaN;
    supraStats.peakUpStroke = NaN;
    supraStats.peakDownStroke = NaN;
    supraStats.peakStrokeRatio = NaN;
    supraStats.fastTrough = NaN;
    supraStats.fastTroughDur = NaN;
    supraStats.slowTrough = NaN;
    supraStats.slowTroughDur = NaN;
    supraStats.latency = NaN;
    supraStats.meanFR50 = NaN;
    supraStats.meanFR100 = NaN;
    supraStats.meanFR250 = NaN;
    supraStats.meanFR500 = NaN;
    supraStats.meanFR750 = NaN;
    supraStats.meanFR1000 = NaN;
    supraStats.peakAdapt = NaN;
    supraStats.ISI = NaN;
    supraStats.instaRate = NaN;
    supraStats.meanISI = NaN;
    supraStats.cvISI = NaN;
    supraStats.adaptIndex = NaN;
    supraStats.adaptIndex2 = NaN;
    supraStats.peakAdapt2 = NaN;
    supraStats.delay = NaN;
    supraStats.burst = NaN;
    supraStats.sweepID = k;
    
    figure('Position',[50 50 1500 250]); set(gcf,'color','w');
    hold on
    plot(LP.V{1,k},'k','LineWidth',0.25)
    xlabel('time-steps')
    ylabel('voltage (mV)')
    axis tight
    export_fig(['D:\genpath\',cellID,' suprathreshold parameters ',int2str(k)],'-pdf','-r100');
    close
end