function [sp,LP] = estimateMaxdVdtNthreshold(LP,sp,k,params)

% estimateMaxdVdtNthreshold

% [1] of maximum change in voltage (used to find threshold)
% [2] initial threshold estimate (see below for more accurate estimate)
dVdt = diff(LP.V{1,k})/LP.acquireRes;                          % derivative

% filter with Gaussian filter

threshold = zeros(1,length(LP.putSpTimes2)); 
thresholdTime = zeros(1,length(LP.putSpTimes2));
countRemoved = 1;
for i = 1:length(LP.putSpTimes2)
    if i > length(LP.putSpTimes2)-countRemoved
        break
    end
    % find max dV/dt
    if length(LP.putSpTimes2) == 1
    	[maxdVdt(i), maxdVdtTime(i)] = max(dVdt(1:sp.peakTime(i))); % max change in voltage
        thresholdTime(i) = find(dVdt(1:maxdVdtTime(i))<(0.05*maxdVdt(i)), 1, 'last');
        threshold(i) = LP.V{1,k}(thresholdTime(i));
    elseif length(LP.putSpTimes2) > 1
        if i == 1
            [maxdVdt(i), maxdVdtTime(i)] = max(dVdt(1:sp.peakTime(i)-1)); % max change in voltage
            thresholdTime(i) = find(dVdt(1:maxdVdtTime(i))<(0.05*maxdVdt(i)), 1, 'last');
            threshold(i) = LP.V{1,k}(thresholdTime(i));
        elseif i > 1
            [maxdVdt(i), maxdVdtTime(i)] = max(dVdt(sp.peakTime(i-1):sp.peakTime(i)-1)); % max change in voltage
            maxdVdtTime(i) = maxdVdtTime(i) + sp.peakTime(i-1);
            if ~isempty(find(dVdt(sp.peakTime(i-1):maxdVdtTime(i))<(0.05*maxdVdt(i)), 1, 'last'))
                thresholdTime(i) = find(dVdt(sp.peakTime(i-1):maxdVdtTime(i))<(0.05*maxdVdt(i)), 1, 'last');
                thresholdTime(i) = thresholdTime(i)+sp.peakTime(i-1);
                threshold(i) = LP.V{1,k}(thresholdTime(i));
                
                % if it gets to this point take an average
                
            else
                maxdVdt(i) = [];
                maxdVdtTime(i) = [];
                LP.putSpTimes2(i) = [];
                sp.peak(i) = [];
                sp.peakTime(i) = [];
                countRemoved = countRemoved + 1;
            end
        end
    end
end

%{
QCpeakNthreshold
%}
 
idx = threshold>params.maxThreshold;
LP.qcRemovals.maxThreshold = LP.putSpTimes2(idx);
LP.putSpTimes2(idx) = [];
sp.peak(idx) = []; sp.peakTime(idx) = [];
threshold(idx) = [];
thresholdTime(idx) = [];
maxdVdt(idx) = [];
maxdVdtTime(idx) = [];

diffthreshold2peak = abs(threshold-sp.peak);
idx2 = diffthreshold2peak < params.minDiffThreshold2Peak;
LP.qcRemovals.diffthreshold2peak = LP.putSpTimes2(idx2);
LP.putSpTimes2(idx2) = [];
sp.peak(idx2) = []; sp.peakTime(idx2) = []; 
threshold(idx2) = [];
thresholdTime(idx2) = [];
maxdVdt(idx2) = [];
maxdVdtTime(idx2) = [];

diffthreshold2peakT = (thresholdTime-sp.peakTime)/LP.acquireRes;
idx3 = diffthreshold2peakT > params.maxDiffThreshold2PeakT;
LP.qcRemovals.diffthreshold2peakT = LP.putSpTimes2(idx3);
LP.putSpTimes2(idx3) = [];
sp.peak(idx3) = []; sp.peakTime(idx3) = [];
threshold(idx3) = [];
thresholdTime(idx3) = [];
maxdVdt(idx3) = [];
maxdVdtTime(idx3) = [];

sp.dVdt = dVdt;
sp.threshold = threshold;
sp.thresholdTime = thresholdTime;
sp.maxdVdt = maxdVdt;
sp.maxdVdtTime = maxdVdtTime;