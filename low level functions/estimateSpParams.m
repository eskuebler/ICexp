function [sp,LP] = estimateSpParams(LP,sp,k,params,cellID,folder)

for i = 1:length(sp.peakTime)
    % AP height based on peak to trough (i.e., potassium dynamics - possibly noisy)
    heightPT(i) = sp.peak(i) - sp.trough(i);
    peakMinusHeight = sp.peak(i)-(heightPT(i)/2);
    temp2 = find(LP.V{1,k}(sp.thresholdRefTime(i):sp.peakTime(i))<peakMinusHeight, 1, 'last');          
    temp2 = sp.thresholdRefTime(i) + temp2;
    if ~isempty(temp2)
        halfHeightTimeUpPT(i) = temp2; 
        temp2 = find(LP.V{1,k}(sp.peakTime(i):sp.troughTime(i))<peakMinusHeight, 1, 'first');          
        temp2 = sp.peakTime(i) + temp2;
        if ~isempty(temp2)
            halfHeightTimeDownPT(i) = temp2(1);
            fullWidthPT(i) = (halfHeightTimeDownPT(i) - halfHeightTimeUpPT(i))*LP.acquireRes;
        else
            halfHeightTimeUpPT(i) = NaN;
            halfHeightTimeDownPT(i) = NaN;
            fullWidthPT(i) = NaN;
        end
    else
        halfHeightTimeUpPT(i) = NaN;
        halfHeightTimeDownPT(i) = NaN;
        fullWidthPT(i) = NaN;
    end

    % AP height based on threshold to peak (i.e., sodium dynamics - possibly reliable)
    heightTP(i) = sp.peak(i) - sp.thresholdRef(i);
    peakMinusHeight = sp.peak(i)-(heightTP(i)/2);
    temp2 = find(LP.V{1,k}(sp.thresholdRefTime(i):sp.peakTime(i))<peakMinusHeight, 1, 'last');          
    temp2 = sp.thresholdRefTime(i) + temp2;
    if ~isempty(temp2)
        halfHeightTimeUpTP(i) = temp2;
        temp2 = find(LP.V{1,k}(sp.peakTime(i):sp.troughTime(i))<peakMinusHeight, 1, 'first');          
        temp2 = sp.peakTime(i) + temp2;
        if ~isempty(temp2)
            halfHeightTimeDownTP(i) = temp2(1);
            fullWidthTP(i) = (halfHeightTimeDownTP(i) - halfHeightTimeUpTP(i))*LP.acquireRes;
        else
            halfHeightTimeUpTP(i) = NaN;
            halfHeightTimeDownTP(i) = NaN;
            fullWidthTP(i) = NaN;
        end
    else
        halfHeightTimeUpTP(i) = NaN;
        halfHeightTimeDownTP(i) = NaN;
        fullWidthTP(i) = NaN;
    end

    % compute peak stroke ratio
    maxTemp = max(sp.dVdt(sp.thresholdRefTime(i):sp.peakTime(i)));
    minTemp = min(sp.dVdt(sp.peakTime(i):sp.troughTime(i)-1));
    if ~isempty(maxTemp) && ~isempty(minTemp)
        peakUpStroke(i) = maxTemp;
        peakDownStroke(i) = minTemp;
        peakStrokeRatio(i) = peakUpStroke(i) / peakDownStroke(i);
    end
    
%     return2zero = find(LP.dVdt{1,k}(sp.peakTime(i)+1:sp.peakTime(i+1)-(2/LP.acquireRes))>=0,1,'first');
    return2zero(i) = find(LP.dVdt{1,k}(sp.peakTime(i)+1:sp.peakTime(i)+1+(100/LP.acquireRes))>=0,1,'first');
    
    % Short (5ms) and long (between events) troughs
%     restingPot = mean(LP.V{1,k}(LP.stimOn(1,k)-(550/LP.acquireRes):LP.stimOn(1,k)-(50/LP.acquireRes))); 
    if i == 1 % if first spike of sweep it is a possible rheobase spike
        if i < length(sp.peakTime) % if there are more spikes to follow we can use the interval
            % first we need to check for large deflections in V w dV/dt to 
            % ensure no spikelets have occurred or spike(s) may have been 
            % omitted due to QC in previous analysis steps (i.e.,
            % threshold)
            findVec = find(LP.dVdt{1,k}(sp.peakTime(i)+1:sp.peakTime(i+1)-(1.5/LP.acquireRes))>20);
            if isempty(findVec) % if there are no large deflections
                if (sp.peakTime(i+1)-sp.peakTime(i))*LP.acquireRes > params.shortTwin
                    [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(params.shortTwin/LP.acquireRes)));
                else
                    [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i+1)));
                end   
                [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i+1)));
                st_t(i) = sp.peakTime(i)+slow_trough_dur(i); % plotting parameter
            else % if there are large deflections
                if findVec(1)*LP.acquireRes > 5 % if the duration from 
                    % the peak to the large deflection < 5ms
                    [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(params.shortTwin/LP.acquireRes)));
                else % otherwise look for minimum prior to large deflection
                    [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+findVec(1)));
                end
                [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+findVec(1)));
                st_t(i) = sp.peakTime(i)+slow_trough_dur(i); % plotting parameter
            end
            clear findVec % clear for next iteration

%             clc
%             disp([slow_trough(i),slow_trough_dur(i)*LP.acquireRes])
%             disp([fast_trough(i),fast_trough_dur(i)*LP.acquireRes])

            figure('Position',[50 50 600 400]); set(gcf,'color','w');
            subplot(2,1,1)
            plot(LP.V{1,k}(sp.peakTime(i)-(5/LP.acquireRes):sp.peakTime(i+1)),'k')
            hold on
            scatter(1+(5/LP.acquireRes),sp.peak(i),'g')
            scatter(fast_trough_dur(i)+(5/LP.acquireRes),fast_trough(i),'r')
            scatter(slow_trough_dur(i)+(5/LP.acquireRes),slow_trough(i),10,'b')
            annotation('textbox',[.3 .7 0 0],'String',fast_trough_dur(i)*LP.acquireRes,'FitBoxToText','on','LineStyle','none','color','r');
            annotation('textbox',[.3 .9 0 0],'String',slow_trough_dur(i)*LP.acquireRes,'FitBoxToText','on','LineStyle','none');
            annotation('textbox',[.52 .7 0 0],'String',fast_trough(i),'FitBoxToText','on','LineStyle','none','color','r');
            annotation('textbox',[.52 .9 0 0],'String',slow_trough(i),'FitBoxToText','on','LineStyle','none');
            xlabel('time')
            ylabel('voltage (mV)')
            axis tight
            subplot(2,1,2)
            hold on
            plot(LP.dVdt{1,k}(sp.peakTime(i)-(5/LP.acquireRes):sp.peakTime(i+1)),'k')
            scatter((5/LP.acquireRes)+1,LP.dVdt{1,k}(sp.peakTime(i)+1),'g')
            scatter(return2zero(i)+(5/LP.acquireRes),0,'r')
            annotation('textbox',[.3 .3 0 0],'String',return2zero(i)*LP.acquireRes,'FitBoxToText','on','LineStyle','none');
            xlabel('time')
            ylabel('dVdt')
            axis tight
            export_fig(['D:\figs2\p2t ',cellID,' ',int2str(k)],'-pdf','-r600');
            close
            
        else % there is only one spike there is no interval to use like above
            [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(params.shortTwin/LP.acquireRes)));
            % again check for large deflections and ensure we look in a
            % window prior to that occurring
            findVec = find(LP.dVdt{1,k}(sp.peakTime(i)+1:LP.stimOff(1,k))>20);
            if isempty(findVec)
                if LP.stimOff(1,k)-sp.peakTime(i) < 300/LP.acquireRes % window end of stim or 300 ms max?
                    [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):LP.stimOff(1,k)));
                else
                    [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(300/LP.acquireRes)));
                end
                st_t(i) = sp.peakTime(i)+slow_trough_dur(i);
            else
                [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+findVec(1)));
                st_t(i) = sp.peakTime(i)+slow_trough_dur(i);
            end
            clear findVec
            
%             clc
%             disp([slow_trough(i),slow_trough_dur(i)*LP.acquireRes])
%             disp([fast_trough(i),fast_trough_dur(i)*LP.acquireRes])

            figure('Position',[50 50 600 400]); set(gcf,'color','w');
            subplot(2,1,1)
            plot(LP.V{1,k}(sp.peakTime(i)-(5/LP.acquireRes):sp.peakTime(i)+(100/LP.acquireRes)),'k')
            hold on
            scatter(1+(5/LP.acquireRes),sp.peak(i),'g')
            scatter(fast_trough_dur(i)+(5/LP.acquireRes),fast_trough(i),'r')
            scatter(slow_trough_dur(i)+(5/LP.acquireRes),slow_trough(i),10,'b')
            annotation('textbox',[.3 .7 0 0],'String',fast_trough_dur(i)*LP.acquireRes,'FitBoxToText','on','LineStyle','none','color','r');
            annotation('textbox',[.3 .9 0 0],'String',slow_trough_dur(i)*LP.acquireRes,'FitBoxToText','on','LineStyle','none');
            annotation('textbox',[.52 .7 0 0],'String',fast_trough(i),'FitBoxToText','on','LineStyle','none','color','r');
            annotation('textbox',[.52 .9 0 0],'String',slow_trough(i),'FitBoxToText','on','LineStyle','none');
            xlabel('time')
            ylabel('voltage (mV)')
            axis tight
            subplot(2,1,2)
            hold on
            plot(LP.dVdt{1,k}(sp.peakTime(i)-(5/LP.acquireRes):sp.peakTime(i)+(100/LP.acquireRes)),'k')
            scatter((5/LP.acquireRes)+1,LP.dVdt{1,k}(sp.peakTime(i)+1),'g')
            scatter(return2zero(i)+(5/LP.acquireRes),0,'r')
            annotation('textbox',[.3 .3 0 0],'String',return2zero(i)*LP.acquireRes,'FitBoxToText','on','LineStyle','none');
            xlabel('time')
            ylabel('dVdt')
            axis tight
            export_fig(['D:\figs2\p2t ',cellID,' ',int2str(k)],'-pdf','-r600');
            close
            
        end
    elseif i < length(sp.peakTime)
        if (sp.peakTime(i+1)-sp.peakTime(i))*LP.acquireRes > params.shortTwin
            [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(params.shortTwin/LP.acquireRes)));
        else
            [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i+1)));
        end   
        [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i+1)));
        st_t(i) = sp.peakTime(i)+slow_trough_dur(i);
        
%         clc
%         disp([slow_trough(i),slow_trough_dur(i)*LP.acquireRes])
%         disp([fast_trough(i),fast_trough_dur(i)*LP.acquireRes])
%         
%         figure('Position',[50 50 600 400]); set(gcf,'color','w');
%         plot(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i+1)),'k')
%         hold on
%         scatter(fast_trough_dur(i),fast_trough(i),'r')
%         scatter(slow_trough_dur(i),slow_trough(i),10,'b')
%         xlabel('time')
%         ylabel('voltage (mV)')
%         axis tight
%         close
    else % final last spike (less care is taken for this spike as much of 
        % what is going on could be confused with stimulus offset and/or a 
        % number of other conductances
        
        [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(params.shortTwin/LP.acquireRes)));
        if LP.stimOff(1,k)-sp.peakTime(i) > 0
            if LP.stimOff(1,k)-sp.peakTime(i) < 300/LP.acquireRes % window end stim or 300 ms?
                [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):LP.stimOff(1,k)));
            else
                [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(300/LP.acquireRes)));
            end
            st_t(i) = sp.peakTime(i)+slow_trough_dur(i);
        else
            slow_trough(i) = NaN;
            slow_trough_dur(i) = NaN;
            st_t(i) = NaN;
        end
%         clc
%         disp([slow_trough(i),slow_trough_dur(i)*LP.acquireRes])
%         disp([fast_trough(i),fast_trough_dur(i)*LP.acquireRes])
%         
%         figure('Position',[50 50 600 400]); set(gcf,'color','w');
%         plot(LP.V{1,k}(sp.peakTime(i):LP.stimOff(1,k)),'k')
%         hold on
%         scatter(fast_trough_dur(i),fast_trough(i),'r')
%         scatter(slow_trough_dur(i),slow_trough(i),10,'b')
%         xlabel('time')
%         ylabel('voltage (mV)')
%         axis tight
%         close
    end
end 

% remove NaNs
% halfHeightTimeUpPT(isnan(halfHeightTimeUpPT)) = [];
% halfHeightTimeDownPT(isnan(halfHeightTimeDownPT)) = [];
% halfHeightTimeUpTP(isnan(halfHeightTimeUpTP)) = [];
% halfHeightTimeDownTP(isnan(halfHeightTimeDownTP)) = [];

if length(heightTP)>1
    idx = heightTP<params.percentRheobaseHeight*heightTP(1);
    LP.qcRemovals.percentRheobaseHeight = LP.putSpTimes2(idx);
    LP.qcRemovals.QCmatpercentRheobaseHeight = heightTP < ...
        params.percentRheobaseHeight * heightTP(1);
    LP.putSpTimes2(idx) = [];
    sp.peak(idx) = []; sp.peakTime(idx) = []; 
    sp.trough(idx) = []; sp.troughTime(idx) = [];
    sp.thresholdRef(idx) = [];
    sp.thresholdRefTime(idx) = [];
    sp.threshold(idx) = [];
    sp.thresholdTime(idx) = [];
    sp.maxdVdt(idx) = [];
    sp.maxdVdtTime(idx) = [];
    heightPT(idx) = [];
    fullWidthPT(idx) = [];
    heightTP(idx) = [];
    fullWidthTP(idx) = [];
    halfHeightTimeUpPT(idx) = [];
    halfHeightTimeDownPT(idx) = [];
    halfHeightTimeUpTP(idx) = [];
    halfHeightTimeDownTP(idx) = [];
    peakUpStroke(idx) = [];
    peakDownStroke(idx) = [];
    peakStrokeRatio(idx) = [];
    fast_trough(idx) = [];
    fast_trough_dur(idx) = [];
    slow_trough(idx) = [];
    slow_trough_dur(idx) = [];
    st_t(idx) = [];
else
    LP.qcRemovals.QCmatpercentRheobaseHeight = zeros(length(LP.putSpTimes2),1);
    LP.qcRemovals.percentRheobaseHeight = NaN;
end

if fullWidthTP(1) <= 0.7        % narrow spiking
    idx = abs(sp.peak-sp.thresholdRef)<params.minDiffThreshold2PeakN;
    LP.qcRemovals.minDiffThreshold2PeakN = LP.putSpTimes2(idx);
    LP.qcRemovals.QCmatminDiffThreshold2PeakN = abs(sp.peak-sp.thresholdRef)<params.minDiffThreshold2PeakN;
    LP.putSpTimes2(idx) = [];
    sp.peak(idx) = []; sp.peakTime(idx) = []; 
    sp.trough(idx) = []; sp.troughTime(idx) = [];
    sp.thresholdRef(idx) = [];
    sp.thresholdRefTime(idx) = [];
    sp.threshold(idx) = [];
    sp.thresholdTime(idx) = [];
    sp.maxdVdt(idx) = [];
    sp.maxdVdtTime(idx) = [];
    heightPT(idx) = [];
    fullWidthPT(idx) = [];
    heightTP(idx) = [];
    fullWidthTP(idx) = [];
    halfHeightTimeUpPT(idx) = [];
    halfHeightTimeDownPT(idx) = [];
    halfHeightTimeUpTP(idx) = [];
    halfHeightTimeDownTP(idx) = [];
    peakUpStroke(idx) = [];
    peakDownStroke(idx) = [];
    peakStrokeRatio(idx) = [];
    fast_trough(idx) = [];
    fast_trough_dur(idx) = [];
    slow_trough(idx) = [];
    slow_trough_dur(idx) = [];
    st_t(idx) = [];
else                                                        % broad spiking
    idx = abs(sp.peak-sp.thresholdRef)<params.minDiffThreshold2PeakB;
    LP.qcRemovals.minDiffThreshold2PeakB = LP.putSpTimes2(idx);
    LP.qcRemovals.QCmatminDiffThreshold2PeakB = abs(sp.peak-sp.thresholdRef)<params.minDiffThreshold2PeakB;
    LP.putSpTimes2(idx) = [];
    sp.peak(idx) = []; sp.peakTime(idx) = []; 
    sp.trough(idx) = []; sp.troughTime(idx) = [];
    sp.thresholdRef(idx) = [];
    sp.thresholdRefTime(idx) = [];
    sp.threshold(idx) = [];
    sp.thresholdTime(idx) = [];
    sp.maxdVdt(idx) = [];
    sp.maxdVdtTime(idx) = [];
    heightPT(idx) = [];
    fullWidthPT(idx) = [];
    heightTP(idx) = [];
    fullWidthTP(idx) = [];
    halfHeightTimeUpPT(idx) = [];
    halfHeightTimeDownPT(idx) = [];
    halfHeightTimeUpTP(idx) = [];
    halfHeightTimeDownTP(idx) = [];
    peakUpStroke(idx) = [];
    peakDownStroke(idx) = [];
    peakStrokeRatio(idx) = [];
    fast_trough(idx) = [];
    fast_trough_dur(idx) = [];
    slow_trough(idx) = [];
    slow_trough_dur(idx) = [];
    st_t(idx) = [];
end

if params.plot_all == 1
    figure('Position',[50 50 600 400]); set(gcf,'color','w');
    plot(LP.V{1,k},'k')
    hold on
    scatter(st_t,slow_trough)
    xlabel('time')
    ylabel('voltage (mV)')
    axis tight
    box off
    export_fig([folder(1:length(folder)-8),cellID,' ',int2str(k),' trough'],'-pdf','-r100');
    close
end

sp.heightPT = heightPT;
sp.halfHeightTimeUpPT = halfHeightTimeUpPT;
sp.halfHeightTimeDownPT = halfHeightTimeDownPT;
sp.fullWidthPT = fullWidthPT;
sp.heightTP = heightTP;
sp.halfHeightTimeUpTP = halfHeightTimeUpTP;
sp.halfHeightTimeDownTP = halfHeightTimeDownTP;
sp.fullWidthTP = fullWidthTP;
sp.peakUpStroke = peakUpStroke;
sp.peakDownStroke = peakDownStroke;
sp.peakStrokeRatio = peakStrokeRatio;
sp.fast_trough = fast_trough;
sp.fast_trough_dur = fast_trough_dur*LP.acquireRes;
sp.slow_trough = slow_trough;
sp.slow_trough_dur = slow_trough_dur*LP.acquireRes;