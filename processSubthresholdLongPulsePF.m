function subStats = processSubthresholdLongPulsePF(LP,params,qc,k,cellID)

%{
processSubthresholdLongPulse
- analysis of subthreshold sweeps
- QC parameters set in processICsweepsParFor
- parameters computed:
    minimum V
    current input
    maximum deflection
    time constant (two different ways)
        - expon fit from rest to minimum
        - expon fit from minimum to steady state
    steady state
    sag
    sag ratio
- finally analysis of rebound spikes
%}

figure('Position',[50 50 1100 250]); set(gcf,'color','w');
hold on

plot(LP.V{1,k},'k')
xlabel('time-steps')
ylabel('voltage (mV)')
axis tight
box off

subStats.subSweepID = k;
subStats.subSweepAmps = LP.sweepAmps(k,1);

% estimate minimum voltage
[subStats.minV,subStats.minVt] = ...
    min(LP.V{1,k}(1,LP.stimOn(1,k):LP.stimOff(1,k)));
subStats.maxSubDeflection = subStats.minV-qc.restVPre;
subStats.minVt = subStats.minVt+LP.stimOn(1,k);

% time constant (rest to minimum V)
y = double(LP.V{1,k}(LP.stimOn(1,k):subStats.minVt)');
x = double(linspace(1,subStats.minVt-LP.stimOn(1,k),length(y))');
if length(y)>=4
    [f,gof] = fit(x,y,'exp2');
    if gof.rsquare > 0.75          % Label NaN if rsquared < 0
        plot(x+LP.stimOn(1,k),f(x),'r-.','LineWidth',2)
        temp = .63*(abs(f(1)-f(length(x))));
        vecin = find(f(1:length(x))<(f(1)-temp), 1, 'first');
        if ~isempty(vecin)
            scatter(vecin(1)+1+LP.stimOn(1,k),LP.V{1,k}(LP.stimOn(1,k))-temp,'r','filled')
            subStats.tauMin = vecin(1)*LP.acquireRes;
            subStats.tauMinamp = LP.sweepAmps(k,1);
            subStats.tauMinGF = 1;
        else
            subStats.tauMinGF = 0;
        end
    else
        subStats.tauMinGF = 0;
    end
end

% sag & sag ratio
subStats.subSteadyState = mean(LP.V{1,k}(LP.stimOff(1,k)-round(50/LP.acquireRes):LP.stimOff(1,k)-1));
subStats.sag = abs(subStats.subSteadyState-subStats.minV);
subStats.sagRatio = subStats.minV/subStats.subSteadyState;

% time constant based on fit between minimum and steady state
y = double(LP.V{1,k}(subStats.minVt:LP.stimOff(1,k)-1)');
x = double(linspace(1,LP.stimOff(1,k)-subStats.minVt,length(y))');
if length(y)>=4
    [f,gof] = fit(x,y,'exp2');
    if gof.rsquare > 0.75          % Label NaN if rsquared < 0
        plot(x+subStats.minVt,f(x),'b-.','LineWidth',2)
        temp = .63*(abs(f(1)-f(length(x))));
        vecin = find(f(1:length(x))>(f(1)+temp), 1, 'first');
        if ~isempty(vecin)
            scatter(vecin(1)+1+subStats.minVt,LP.V{1,k}(subStats.minVt)+temp,'b','filled')
            subStats.tauSS = vecin(1)*LP.acquireRes;
            subStats.tauSSamp = LP.sweepAmps(k,1);
            subStats.tauSSGF = 1;
        else
            subStats.tauSSGF = 0;
        end
    else
        subStats.tauSSGF = 0;
    end
end

plot((LP.stimOff(1,k)-round(50/LP.acquireRes):LP.stimOff(1,k)-1),...
    LP.V{1,k}(LP.stimOff(1,k)-round(50/LP.acquireRes):LP.stimOff(1,k)-1),'g-.','LineWidth',2)

% rebound slope
[val,loc] = max(LP.V{1,k}(LP.stimOff(1,k):...
  LP.stimOff(1,k)+round(params.reboundWindow/LP.acquireRes)));
x = (loc:loc+round(params.reboundFitWindow/LP.acquireRes))-loc;
[f,~] = polyfit(x,LP.V{1,k}(LP.stimOff(1,k)+loc:...
	LP.stimOff(1,k)+loc+round(params.reboundFitWindow/LP.acquireRes))',1);
subStats.reboundSlope = f(1);
subStats.reboundDepolarization = abs(LP.V{1,k}(LP.stimOff(1,k)+loc)-...
    LP.V{1,k}(LP.stimOff(1,k)+loc+round(params.reboundFitWindow/LP.acquireRes)));
plot(x+loc+LP.stimOff(1,k),(f(1)*x+f(2))','c-.','LineWidth',2)
scatter(loc+LP.stimOff(1,k),val,'g','filled')
scatter(round(params.reboundFitWindow/LP.acquireRes)+loc+LP.stimOff(1,k),mean(LP.V{1,k}(end-(3/LP.acquireRes):end)),'g','filled')

% save figure
export_fig(['D:\genpath\',cellID,' ',int2str(k),' subthreshold parameters'],'-pdf','-r100');
close

% rebound spikes
reboundAPTimes = find(LP.V{1,k}(LP.stimOff(1,k):...
    LP.stimOff(1,k)+(params.reboundSpWindow/LP.acquireRes))>=params.thresholdV); 
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
    subStats.reboundAPs = rebound2APTimes;
else
    subStats.reboundAPs = NaN;
end
