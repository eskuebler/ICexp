%{
summarySubthresholdLP
%}

if LP.qcSubSweeps > 2
    


    % resistance
    subSweeps = find(LP.subSweepAmps>-100);     % find subthreshold sweeps > -100 pA
    y = LP.minV(subSweeps)';                    % minimum voltage
    x = LP.subSweepAmps(subSweeps)';            % current input
    f = polyfit(x,y,1);
    LP.resistance = f(1) * (10^3);
    subplot(1,3,1)
    hold on
    plot(x,(f(1)*x+f(2))','k','LineWidth',1)
    scatter(x,y,'r')
    legend('off')
    xlabel('input current (pA)')
    ylabel('membrane potential (mV)')
    box off
    axis tight
    clear subSweeps x y f

    % time constant (rest to minimum)
    if isfield(LP,'tauMinamp') && length(LP.tauMinamp) > 2
        subSweeps = find(LP.tauMinamp>-100);        % find subthreshold sweeps > -100 pA
        y = LP.tauMin(subSweeps)';                  % time constant estimate
        x = LP.tauMinamp(subSweeps)';               % current input
        LP.tauMinAvg = mean(y);
        f = polyfit(x,y,1);
        subplot(1,3,2)
        hold on
        plot(x,(f(1)*x+f(2))','k','LineWidth',1)
        scatter(x,y,'r')
        legend('off')
        xlabel('input current (pA)')
        ylabel('membrane potential (mV)')
        box off
        axis tight
        ylim([0 100])
        clear subSweeps x y f
    end

    % time constant (minimum to steady state)
    if isfield(LP,'tauSSamp') && length(LP.tauSSamp) > 2
        subSweeps = find(LP.tauSSamp>-100);         % find subthreshold sweeps > -100 pA
        y = LP.tauSS(subSweeps)';                   % time constant estimate
        x = LP.tauSSamp(subSweeps)';                % current input
        LP.tauMinAvg = mean(y);
        f = polyfit(x,y,1);
        subplot(1,3,3)
        hold on
        plot(x,(f(1)*x+f(2))','k','LineWidth',1)
        scatter(x,y,'r')
        legend('off')
        xlabel('input current (pA)')
        ylabel('membrane potential (mV)')
        box off
        axis tight
        ylim([0 100])
        clear subSweeps x y f
    end

    % save figure
    export_fig(['D:\genpath\',cellID,' subthreshold summary'],'-pdf','-r100');
    close

end