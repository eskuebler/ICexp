function [sp,LP] = estimateSpParams(LP,sp,k,params)

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
    
    % Short (5ms) and long (between events) troughs
    restingPot = mean(LP.V{1,k}(LP.stimOn(1,k)-(550/LP.acquireRes):LP.stimOn(1,k)-(50/LP.acquireRes)));
    [fast_trough(i),fast_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(5/LP.acquireRes)));
    if i < length(sp.peakTime)
        [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i+1)));
    else
        [slow_trough(i),slow_trough_dur(i)] = min(LP.V{1,k}(sp.peakTime(i):sp.peakTime(i)+(100/LP.acquireRes)));
    end
end 

% remove NaNs
halfHeightTimeUpPT(isnan(halfHeightTimeUpPT)) = [];
halfHeightTimeDownPT(isnan(halfHeightTimeDownPT)) = [];
halfHeightTimeUpTP(isnan(halfHeightTimeUpTP)) = [];
halfHeightTimeDownTP(isnan(halfHeightTimeDownTP)) = [];

if length(heightPT)>1
    idx = heightPT<params.percentRheobaseHeight*heightPT(1);
    LP.qcRemovals.percentRheobaseHeight = LP.putSpTimes2(idx);
    LP.qcRemovals.QCmatpercentRheobaseHeight = heightPT < ...
        params.percentRheobaseHeight * heightPT(1);
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
    peakUpStroke(idx) = [];
    peakDownStroke(idx) = [];
    peakStrokeRatio(idx) = [];
    fast_trough(idx) = [];
    fast_trough_dur(idx) = [];
    slow_trough(idx) = [];
    slow_trough_dur(idx) = [];
else
    LP.qcRemovals.QCmatpercentRheobaseHeight = zeros(length(LP.putSpTimes2),1);
    LP.qcRemovals.percentRheobaseHeight = NaN;
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