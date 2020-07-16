% estimateAPparameters
if voltage_data(end)==0 || s.stimType(k,N.goodSweeps)==3
    for i = 1:length(Put2APTimes)
        % AP height based on peak to trough (i.e., potassium dynamics - possibly noisy)
        heightPT(i) = peak(i) - trough(i);
        peakMinusHeight = peak(i)-(heightPT(i)/2);
        temp2 = find(voltage_data(thresholdRefTime(i):peakTime(i))<peakMinusHeight, 1, 'last');          
        temp2 = thresholdRefTime(i) + temp2;
        halfHeightTimeUpPT(i) = temp2;
        temp2 = find(voltage_data(peakTime(i):troughTime(i))<peakMinusHeight, 1, 'first');          
        temp2 = peakTime(i) + temp2;
        halfHeightTimeDownPT(i) = temp2(1);
        fullWidthPT(i) = (halfHeightTimeDownPT(i) - halfHeightTimeUpPT(i))*acquireRes;

        % AP height based on threshold to peak (i.e., sodium dynamics - possibly reliable)
        heightTP(i) = peak(i) - thresholdRef(i);
        peakMinusHeight = peak(i)-(heightTP(i)/2);
        temp2 = find(voltage_data(thresholdRefTime(i):peakTime(i))<peakMinusHeight, 1, 'last');          
        temp2 = thresholdRefTime(i) + temp2;
        halfHeightTimeUpTP(i) = temp2;
        temp2 = find(voltage_data(peakTime(i):troughTime(i))<peakMinusHeight, 1, 'first');          
        temp2 = peakTime(i) + temp2;
        halfHeightTimeDownTP(i) = temp2(1);
        fullWidthTP(i) = (halfHeightTimeDownTP(i) - halfHeightTimeUpTP(i))*acquireRes;

        % compute peak stroke ratio
        peakUpStroke(i) = max(dVdt(thresholdRefTime(i):peakTime(i)));
        peakDownStroke(i) = min(dVdt(peakTime(i):troughTime(i)-1));
        peakStrokeRatio(i) = peakUpStroke(i) / peakDownStroke(i);
    end
    clear i temp temp2 dVdt
else
    if s.stimType(k,N.goodSweeps)==1
        for i = 1:length(Put2APTimes)
            % AP height based on peak to trough (i.e., potassium dynamics - possibly noisy)
            heightPT(i) = peak(i) - trough(i);
            peakMinusHeight = peak(i)-(heightPT(i)/2);
            temp2 = find(voltage_data(thresholdRefTime(i):peakTime(i))<peakMinusHeight, 1, 'last');          
            temp2 = thresholdRefTime(i) + temp2;
            if ~isempty(temp2)
                halfHeightTimeUpPT(i) = temp2; 
                temp2 = find(voltage_data(peakTime(i):troughTime(i))<peakMinusHeight, 1, 'first');          
                temp2 = peakTime(i) + temp2;
                if ~isempty(temp2)
                    halfHeightTimeDownPT(i) = temp2(1);
                    fullWidthPT(i) = (halfHeightTimeDownPT(i) - halfHeightTimeUpPT(i))*acquireRes;
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
            heightTP(i) = peak(i) - thresholdRef(i);
            peakMinusHeight = peak(i)-(heightTP(i)/2);
            temp2 = find(voltage_data(thresholdRefTime(i):peakTime(i))<peakMinusHeight, 1, 'last');          
            temp2 = thresholdRefTime(i) + temp2;
            if ~isempty(temp2)
                halfHeightTimeUpTP(i) = temp2;
                temp2 = find(voltage_data(peakTime(i):troughTime(i))<peakMinusHeight, 1, 'first');          
                temp2 = peakTime(i) + temp2;
                if ~isempty(temp2)
                    halfHeightTimeDownTP(i) = temp2(1);
                    fullWidthTP(i) = (halfHeightTimeDownTP(i) - halfHeightTimeUpTP(i))*acquireRes;
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
            maxTemp = max(dVdt(thresholdRefTime(i):peakTime(i)));
            minTemp = min(dVdt(peakTime(i):troughTime(i)-1));
            if ~isempty(maxTemp) && ~isempty(minTemp)
                peakUpStroke(i) = maxTemp;
                peakDownStroke(i) = minTemp;
                peakStrokeRatio(i) = peakUpStroke(i)  / peakDownStroke(i);
            end
            % After hyper-polarization potentials & repolarization time
            restingPot = mean(voltage_data(stimOn-(600/acquireRes):stimOn-(100/acquireRes)));
            [fast_trough(i),fast_trough_dur(i)] = min(voltage_data(peakTime(i):peakTime(i)+(5/acquireRes)));
            if i < length(Put2APTimes)
                [slow_trough(i),slow_trough_dur(i)] = min(voltage_data(peakTime(i):peakTime(i+1)));
            else
                [slow_trough(i),slow_trough_dur(i)] = min(voltage_data(peakTime(i):peakTime(i)+(150/acquireRes)));
            end
        end 
        clear i temp temp2 dVdt
        
    elseif s.stimType(k,N.goodSweeps)==2    
        for i = 1:length(Put2APTimes)
            % AP height based on peak to trough (i.e., potassium dynamics - possibly noisy)
            heightPT(i) = peak(i) - trough(i);
            peakMinusHeight = peak(i)-(heightPT(i)/2);
            temp2 = find(voltage_data(thresholdRefTime(i):peakTime(i))<peakMinusHeight, 1, 'last');          
            temp2 = thresholdRefTime(i) + temp2;
            if ~isempty(temp2)
                halfHeightTimeUpPT(i) = temp2; 
                temp2 = find(voltage_data(peakTime(i):troughTime(i))<peakMinusHeight, 1, 'first');          
                temp2 = peakTime(i) + temp2;
                if ~isempty(temp2)
                    halfHeightTimeDownPT(i) = temp2(1);
                    fullWidthPT(i) = (halfHeightTimeDownPT(i) - halfHeightTimeUpPT(i))*acquireRes;
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
            heightTP(i) = peak(i) - thresholdRef(i);
            peakMinusHeight = peak(i)-(heightTP(i)/2);
            temp2 = find(voltage_data(thresholdRefTime(i):peakTime(i))<peakMinusHeight, 1, 'last');          
            temp2 = thresholdRefTime(i) + temp2;
            if ~isempty(temp2)
                halfHeightTimeUpTP(i) = temp2;
                temp2 = find(voltage_data(peakTime(i):troughTime(i))<peakMinusHeight, 1, 'first');          
                temp2 = peakTime(i) + temp2;
                if ~isempty(temp2)
                    halfHeightTimeDownTP(i) = temp2(1);
                    fullWidthTP(i) = (halfHeightTimeDownTP(i) - halfHeightTimeUpTP(i))*acquireRes;
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
            maxTemp = max(dVdt(thresholdRefTime(i):peakTime(i)));
            minTemp = min(dVdt(peakTime(i):troughTime(i)-1));
            if ~isempty(maxTemp) && ~isempty(minTemp)
                peakUpStroke(i) = maxTemp;
                peakDownStroke(i) = minTemp;
                peakStrokeRatio(i) = peakUpStroke(i)  / peakDownStroke(i);
            end
            
            % After hyper-polarization potentials & repolarization time
            restingPot = mean(voltage_data(stimOn-(600/acquireRes):stimOn-(100/acquireRes)));
            [fast_trough(i),fast_trough_dur(i)] = min(voltage_data(peakTime(i):peakTime(i)+(5/acquireRes)));
            if i < length(Put2APTimes)
                [slow_trough(i),slow_trough_dur(i)] = min(voltage_data(peakTime(i):peakTime(i+1)));
            end
            %{
            hold on
            plot(voltage_data(peakTime(i)-(30/acquireRes):end),'k')
            scatter((30/acquireRes),peak,'k','MarkerFaceColor','w')
            scatter(fast_trough_dur(i)+(30/acquireRes),fast_trough(i),'k','MarkerFaceColor','w')
            if exist('slow_trough')==1
                scatter(slow_trough_dur(i)+(30/acquireRes),slow_trough(i),'k','MarkerFaceColor','w')
            end
            
            if exist('slow_trough')==1 && (slow_trough(i)-restingPot)<1
                repolarizationT(i) = (slow_trough_dur(i)-fast_trough_dur(i));
                repolarizationTmV(i) = voltage_data(peakTime(i)+slow_trough_dur(i));
            else
                if i==1 && length(Put2APTimes) > 1
                    targetT = peakTime(i)+fast_trough_dur(i);
                    vec = find(abs(voltage_data(targetT:peakTime(i+1))-restingPot)<1); 
                    repolarizationT = vec(1);
    %                 scatter(loc,voltage_data(targetT-(30/acquireRes)+loc),'k','MarkerFaceColor','w')
                elseif i==1 && length(Put2APTimes) == 1
                    targetT = peakTime(i)+fast_trough_dur(i);
                    vec = find(abs(voltage_data(targetT:end)-restingPot)<1);
                    repolarizationT = vec(1);
    %                 scatter(loc,voltage_data(targetT-(30/acquireRes)+loc),'k','MarkerFaceColor','w')
                end
            end
            pause(1)
            close
            %}
        end 
        clear i temp temp2 dVdt
    end
    fast_trough_dur = fast_trough_dur * acquireRes;
    if exist('slow_trough_dur')==1
        slow_trough_dur = slow_trough_dur * acquireRes;
    end
end
    
% remove NaNs
halfHeightTimeUpPT(isnan(halfHeightTimeUpPT)) = [];
halfHeightTimeDownPT(isnan(halfHeightTimeDownPT)) = [];
halfHeightTimeUpTP(isnan(halfHeightTimeUpTP)) = [];
halfHeightTimeDownTP(isnan(halfHeightTimeDownTP)) = [];

if params.plotAcc == 1

    xlabel('time')
    ylabel('voltage (mV)')
    axis tight;
	xlim([stimOn-(50/acquireRes) stimOff+(50/acquireRes)])									% zoom in on APs
    box off;
end
