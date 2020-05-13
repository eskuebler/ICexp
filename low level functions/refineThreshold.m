function [sp] = refineThreshold(LP,sp,k)

% refined threshold (based on average max dV/dt)

thresholdRef = zeros(1,length(LP.putSpTimes2)); 
thresholdRefTime = zeros(1,length(LP.putSpTimes2));
AVGmaxdVdt = mean(sp.maxdVdt);                                              % max average

for i = 1:length(LP.putSpTimes2)                                            % for each spike time
    if length(LP.putSpTimes2) == 1
        thresholdRefTime(i) = find(sp.dVdt(1:sp.maxdVdtTime(i)) < ...
            (0.05*AVGmaxdVdt), 1, 'last');                                  % 5% of max average dV/dt
        thresholdRef(i) = LP.V{1,k}(thresholdRefTime(i));                   % voltage at threshold
        hold on
        plot(LP.V{1,k})
        plot(sp.dVdt)
        scatter([1 10],[AVGmaxdVdt AVGmaxdVdt],'b')
        scatter(thresholdRefTime(i),thresholdRef(i))
        scatter(thresholdTime(i),threshold(i))
        xlim([thresholdTime(i)-10 sp.peakTime(i)+10])
        pause(1)
        close
    else
        if i == 1
            [~, maxdVdtTime(i)] = max(sp.dVdt(1:sp.peakTime(i)-1));         % max change in voltage
            thresholdRefTime(i) = find(sp.dVdt(1:maxdVdtTime(i)) < ...
                (0.05*AVGmaxdVdt), 1, 'last');
            thresholdRef(i) = LP.V{1,k}(thresholdRefTime(i));
            
%             hold on
%             plot(LP.V{1,k})
%             plot(dVdt)
%             scatter(maxdVdtTime(i),maxdVdt(i))
%             scatter(thresholdTime(i),threshold(i))
%             scatter(thresholdTime(i),dVdt(thresholdTime(i)))
%             xlim([thresholdTime(i)-10 sp.peakTime(i)+10])
%             pause(1)
%             close
        elseif i > 1
            
            % this needs the same if statements
            
            [maxdVdt(i), maxdVdtTime(i)] = max(sp.dVdt(sp.peakTime(i-1):sp.peakTime(i)-1)); % max change in voltage
            maxdVdtTime(i) = maxdVdtTime(i) + sp.peakTime(i-1);
            if ~isempty(find(sp.dVdt(sp.peakTime(i-1):maxdVdtTime(i))<(0.05*AVGmaxdVdt), 1, 'last'))
                thresholdRefTime(i) = find(sp.dVdt(sp.peakTime(i-1):maxdVdtTime(i))<(0.05*AVGmaxdVdt), 1, 'last');
                thresholdRefTime(i) = thresholdRefTime(i)+sp.peakTime(i-1);
                thresholdRef(i) = LP.V{1,k}(thresholdRefTime(i));
            else
                'yo'
                thresholdRefTime(i) = find(sp.dVdt(sp.peakTime(i-1):maxdVdtTime(i))<(0.05*AVGmaxdVdt), 1, 'last');
                thresholdRefTime(i) = thresholdRefTime(i)+sp.peakTime(i-1);
                thresholdRef(i) = LP.V{1,k}(thresholdRefTime(i));
            end
        end
    end
end

sp.thresholdRef = thresholdRef;
sp.thresholdRefTime = thresholdRefTime;
% b.maxdVdt = maxdVdt;
% b.maxdVdtTime = maxdVdtTime;