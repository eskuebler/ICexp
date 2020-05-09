function [sp] = refineThreshold(LP,sp,k)

% refined threshold (based on average max dV/dt)

thresholdRef = zeros(1,length(LP.putSpTimes2)); thresholdRefTime = zeros(1,length(LP.putSpTimes2));
AVGmaxdVdt = mean(sp.maxdVdt); 
for i = 1:length(LP.putSpTimes2)
    if length(LP.putSpTimes2) == 1
        thresholdRefTime(i) = find(sp.dVdt(1:sp.maxdVdtTime(i))<(0.05*AVGmaxdVdt), 1, 'last');
        thresholdRef(i) = LP.V{1,k}(thresholdRefTime(i));
    elseif length(LP.putSpTimes2) > 1
        if i == 1
            [maxdVdt(i), maxdVdtTime(i)] = max(sp.dVdt(1:sp.peakTime(i)-1)); % max change in voltage
            thresholdRefTime(i) = find(sp.dVdt(1:maxdVdtTime(i))<(0.05*AVGmaxdVdt), 1, 'last');
            thresholdRef(i) = LP.V{1,k}(thresholdRefTime(i));
        elseif i > 1
            [maxdVdt(i), maxdVdtTime(i)] = max(sp.dVdt(sp.peakTime(i-1):sp.peakTime(i)-1)); % max change in voltage
            maxdVdtTime(i) = maxdVdtTime(i) + sp.peakTime(i-1);
            thresholdRefTime(i) = find(sp.dVdt(sp.peakTime(i-1):maxdVdtTime(i))<(0.05*AVGmaxdVdt), 1, 'last');
            thresholdRefTime(i) = thresholdRefTime(i)+sp.peakTime(i-1);
            thresholdRef(i) = LP.V{1,k}(thresholdRefTime(i));
        end
    end
end

sp.thresholdRef = thresholdRef;
sp.thresholdRefTime = thresholdRefTime;
% b.maxdVdt = maxdVdt;
% b.maxdVdtTime = maxdVdtTime;