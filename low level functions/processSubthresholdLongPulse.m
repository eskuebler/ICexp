%{
processSubthresholdLongPulse
%}

figure('Position',[50 50 1500 250]); set(gcf,'color','w');
subplot(1,5,1)
plot(LP.V{1,k})
xlabel('time-steps')
ylabel('voltage (mV)')
title('raw subthreshold trace')
axis tight
box off

LP.subSweepAmps(LP.qcSubSweeps) = LP.sweepAmps(k,1);

% estimate minimum voltage
[LP.minV(LP.qcSubSweeps),LP.minVt(LP.qcSubSweeps)] = min(LP.V{1,k}(1,LP.stimOn(1,k):LP.stimOff(1,k)));
LP.maxSubDeflection(LP.qcSubSweeps) = LP.minV(LP.qcSubSweeps)-LP.restVPre(1,k);
LP.minVt(LP.qcSubSweeps) = LP.minVt(LP.qcSubSweeps)+LP.stimOn(1,k);

% time constant (rest to minimum V)
y = double(LP.V{1,k}(LP.stimOn(1,k):LP.minVt(LP.qcSubSweeps))');
x = linspace(1,LP.minVt(LP.qcSubSweeps)-LP.stimOn(1,k),length(y))';
subplot(1,5,2)
hold on
xlabel('time-steps')
ylabel('voltage (mV)')
title('expon fit to minimum')
axis tight
box off
if length(y)>=4
    [f,gof] = fit(x,y,'exp2');
    if gof.rsquare > 0.75          % Label NaN if rsquared < 0
        plot(f,x,y,'k-')
        temp = .63*(abs(f(1)-f(length(x))));
        vecin = find(f(1:length(x))<(f(1)-temp), 1, 'first');
        scatter(vecin+1,LP.V{1,k}(LP.stimOn(1,k))-temp)
        legend({'data','model','63%'},'location','northeast')
        if ~isempty(vecin)
            LP.tauMin(LP.qcSubSweeps) = vecin*LP.acquireRes;
            LP.tauMinamp(LP.qcSubSweeps) = LP.sweepAmps(k,1);
        end
    else
        plot(x,y)
    end
end
clear x y f temp gof vecin

% sag & sag ratio
LP.subSteadyState(LP.qcSubSweeps) = mean(LP.V{1,k}(LP.stimOff(1,k)-round(50/LP.acquireRes):LP.stimOff(1,k)-1));
LP.sag(LP.qcSubSweeps) =  abs(LP.subSteadyState(LP.qcSubSweeps)-LP.minV(LP.qcSubSweeps));
LP.sagRatio(LP.qcSubSweeps) = LP.minV(LP.qcSubSweeps)/LP.subSteadyState(LP.qcSubSweeps);

subplot(1,5,3)
hold on
plot((LP.minVt(LP.qcSubSweeps):LP.stimOff(1,k)-1),...
    LP.V{1,k}(LP.minVt(LP.qcSubSweeps):LP.stimOff(1,k)-1))
plot((LP.stimOff(1,k)-round(50/LP.acquireRes):LP.stimOff(1,k)-1),...
    LP.V{1,k}(LP.stimOff(1,k)-round(50/LP.acquireRes):LP.stimOff(1,k)-1))
xlabel('time-steps')
ylabel('voltage (mV)')
title('steady state & sag')
axis tight
box off

% time constant based on fit between minimum and steady state
y = double(LP.V{1,k}(LP.minVt(LP.qcSubSweeps):LP.stimOff(1,k)-1)');
x = linspace(1,LP.stimOff(1,k)-LP.minVt(LP.qcSubSweeps),length(y))';
subplot(1,5,4)
hold on
xlabel('time-steps')
ylabel('voltage (mV)')
title('expon fit to steady state')
axis tight
box off
if length(y)>=4
    [f,gof] = fit(x,y,'exp2');
    if gof.rsquare > 0.75          % Label NaN if rsquared < 0
        plot(f,x,y,'k-')
        temp = .63*(abs(f(1)-f(length(x))));
        vecin = find(f(1:length(x))>(f(1)+temp), 1, 'first');
        scatter(vecin+1,LP.V{1,k}(LP.minVt(LP.qcSubSweeps))+temp)
        legend({'data','model','63%'},'location','southeast')
        if ~isempty(vecin)
            LP.tauSS(LP.qcSubSweeps) = vecin*LP.acquireRes;
            LP.tauSSamp(LP.qcSubSweeps) = LP.sweepAmps(k,1);
        end
    else
        plot(x,y)
    end
end
clear x y f temp gof vecin

% rebound slope
[val,loc] = max(LP.V{1,k}(LP.stimOff(1,k):LP.stimOff(1,k)+round(params.reboundWindow/LP.acquireRes)));
x = (loc:loc+round(params.reboundFitWindow/LP.acquireRes))-loc;
[f,gof] = polyfit(x,LP.V{1,k}(LP.stimOff+loc:LP.stimOff+loc+round(params.reboundFitWindow/LP.acquireRes))',1);
LP.reboundSlope(LP.qcSubSweeps) = f(1);
LP.reboundDepolarization(LP.qcSubSweeps) = abs(LP.V{1,k}(LP.stimOff(1,k)+loc)-...
    LP.V{1,k}(LP.stimOff(1,k)+loc+round(params.reboundFitWindow/LP.acquireRes)));

% plot rebound slope
subplot(1,5,5)
hold on
plot(LP.V{1,k}(LP.stimOff(1,k):LP.stimOff(1,k)+loc+round(params.reboundFitWindow/LP.acquireRes)))
scatter(loc,val)
scatter(round(params.reboundFitWindow/LP.acquireRes)+loc,mean(LP.V{1,k}(end-(3/LP.acquireRes):end)))
plot(x+loc,(f(1)*x+f(2))','k','LineWidth',1)
title('rebound slope')
ylabel('mV')
xlabel('time')
axis tight
box off
clear val loc x f

% save figure
export_fig(['D:\genpath\',cellID,' subthreshold parameters ',int2str(k)],'-pdf','-r100');
close

% rebound spikes
reboundAPTimes = find(LP.V{1,k}(LP.stimOff(1,k):LP.stimOff(1,k)+(params.reboundSpWindow/LP.acquireRes))>=params.thresholdV); 
if ~isempty(reboundAPTimes)                 % if no spikes
    diffPutAPTime = diff(reboundAPTimes);
    rebound2APTimes = [];
    tag = 1;
    dCount = 1;
    for i = 1:length(reboundAPTimes)-1
        if diffPutAPTime(i) ~= 1
            int4Peak{dCount} = reboundAPTimes(tag):reboundAPTimes(i);
            rebound2APTimes(dCount) = reboundAPTimes(tag);
            tag = i+1;
            dCount = dCount + 1;
        end
    end
    int4Peak{dCount} = reboundAPTimes(tag):reboundAPTimes(end);
    rebound2APTimes(dCount) = reboundAPTimes(tag);
    LP.reboundAPs(LP.qcSubSweeps) = rebound2APTimes;
    clear diffPutAPTime tag dCount i reboundAPTimes rebound2APTimes int4Peak
end
clear reboundAPTimes

LP.qcSubSweepIdx(LP.qcSubSweeps,1) = k;
LP.qcSubSweeps = LP.qcSubSweeps + 1;
