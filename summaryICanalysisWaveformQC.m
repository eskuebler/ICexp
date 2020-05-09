%{
summaryICanalysis
%}

clear; close all; clc;

load('pyramidal cell IDs.mat')
pyrID = ID; clear ID
load('cell_types_specimen_details.mat','donor__species','specimen__id',...
    'line_name','tag__dendrite_type','ef__ri','ef__tau',...
    'ef__threshold_i_long_square','cell_reporter_status',...
    'structure__acronym');

cellList = dir('D:\genpath\*.mat');

input_current_s = zeros(length(cellList),25);
firing_rate_s = zeros(length(cellList),25);
for n = 1:length(cellList)
    clc; disp(n)
    if length(cellList(n).name) == 17
        cellID = cellList(n).name(1:length(cellList(n).name)-4);
        if sum(cellID(5:6)=='MJ')==2
            ID(n,:) = cellID;
            continue
        end
    else
        cellID = cellList(n).name(10:length(cellList(n).name)-4);
    end
    if length(cellID) == 9
        % match species ID # to get Allen parameters
        k  = find(specimen__id==str2num(cellID));

        % get Allen Institute parameters
        ID(n,1) = specimen__id(k);
        specimen(n,1) = donor__species(k);
        struct(n,1) = structure__acronym(k);
        transline(n,1) = line_name(k);
        reporterStatus(n,1) = cell_reporter_status(k);
        dendrite_type(n,1) = tag__dendrite_type(k);
        AIr(n,1) = ef__ri(k);
        AItau(n,1) = ef__tau(k);
        AIrheobase(n,1) = round(ef__threshold_i_long_square(k));
        if sum(ismember(pyrID,ID(n,1)))
            pyramidalID(n,1) = 1;
        else
            pyramidalID(n,1) = 0;
        end
    else
        ID(n,:) = cellID;
        specimen(n,:) = 'NHP';
        struct(n,:) = 'PFC';
        transline(n,:) = 'N/A';
        reporterStatus(n,:) = 'N/A';
        dendrite_type(n,:) = 'N/A';
        AIr(n,1) = NaN;
        AItau(n,1) = NaN;
        AIrheobase(n,1) = NaN;
        pyramidalID(n,1) = NaN;
    end

    % load our analysis
    load(['D:\genpath\',cellList(n).name]);
    
    % subthreshold summary parameters
    resistance(n,1) = a.LP.subSummary.resistance;
    time_constant(n,1) = a.LP.subSummary.tauMin;
    time_constant2(n,1) = a.LP.subSummary.tauSS;
    if a.LP.fullStruct == 1
        k = find(round(a.LP.sweepAmps) == -90);                % find -90 pA input
        if length(k)>1
            k = k(1);
        end
        if ~isempty(k) && isfield(a.LP,'stats') && length(a.LP.stats)>k && ~isempty(a.LP.stats{k,1})
            rebound_slope(n,1) = 0;
            rebound_depolarization(n,1) = 0;
            nb_rebound_sp(n,1) = length(a.LP.stats{k,1}.reboundAPs);
            sag(n,1) = a.LP.stats{k,1}.sag;
            sag_ratio(n,1) = a.LP.stats{k,1}.sagRatio;
            Vrest(n,1) = a.LP.stats{k,1}.qc.restVPre;
        else
            k = find(round(a.LP.sweepAmps) == -70);
            if length(k)>1
                k = k(1);
            end
            if ~isempty(k) && isfield(a.LP,'stats') && length(a.LP.stats)>k && ~isempty(a.LP.stats{k,1})
                rebound_slope(n,1) = 0;
                rebound_depolarization(n,1) = 0;
                nb_rebound_sp(n,1) = length(a.LP.stats{k,1}.reboundAPs);
                sag(n,1) = a.LP.stats{k,1}.sag;
                sag_ratio(n,1) = a.LP.stats{k,1}.sagRatio;
                Vrest(n,1) = a.LP.stats{k,1}.qc.restVPre;
            else
                rebound_slope(n,1) = NaN;
                rebound_depolarization(n,1) = NaN;
                nb_rebound_sp(n,1) = NaN;
                sag(n,1) = NaN;
                sag_ratio(n,1) = NaN;
                Vrest(n,1) = NaN;
            end
        end

        % find rheobase and parameters of first spike
        if isfield(a.LP,'stats')
            [B,I] = sort(round(a.LP.sweepAmps(1:length(a.LP.stats))));
            temp = find(B>0);
            idx = I(temp);
            amp = B(temp);
            input_current_s(n,1:length(amp)) = round(a.LP.sweepAmps(idx));
            for k = 1:length(idx)
                if ~isempty(a.LP.stats{idx(k),1}) && sum(~isnan(a.LP.stats{idx(k),1}.spTimes))>0
                    firing_rate_s(n,k) = a.LP.stats{idx(k),1}.meanFR1000;
                end
            end
            for k = 1:length(idx)
                if ~isempty(a.LP.stats{idx(k),1}) && sum(~isnan(a.LP.stats{idx(k),1}.spTimes))>0
                    % QC vector
                    broad = a.LP.stats{idx(k),1}.fullWidthPT(1) > 1.5;
                    lowPeak = a.LP.stats{idx(k),1}.peak(1) < 10;
                    lowHeight = a.LP.stats{idx(k),1}.heightPT(1) < 40;
                    qcVec = broad&&lowPeak||lowHeight;
                    if qcVec == 0
                        rheobaseLP(n,1) = amp(k);
                        delay(n,1) = a.LP.stats{idx(k),1}.delay;
                        burst(n,1) = a.LP.stats{idx(k),1}.burst;
                        latency(n,1) = a.LP.stats{idx(k),1}.latency;
                        cv_ISI(n,1) = a.LP.stats{idx(k),1}.cvISI;
                        adaptation1(n,1) = a.LP.stats{idx(k),1}.adaptIndex;
                        adaptation2(n,1) = a.LP.stats{idx(k),1}.adaptIndex2;
                        peak_adaptation1(n,1) = a.LP.stats{idx(k),1}.peakAdapt;
                        peak_adaptation2(n,1) = a.LP.stats{idx(k),1}.peakAdapt2;
                        peakLP(n,1) = a.LP.stats{idx(k),1}.peak(1);
                        thresholdLP(n,1) = a.LP.stats{idx(k),1}.thresholdRef(1);
                        half_width_threshold_peak(n,1) = a.LP.stats{idx(k),1}.fullWidthTP(1);
                        half_width_peak_trough(n,1) = a.LP.stats{idx(k),1}.fullWidthPT(1);
                        height_threshold_peak(n,1) = a.LP.stats{idx(k),1}.heightTP(1);
                        height_peak_trough(n,1) = a.LP.stats{idx(k),1}.heightPT(1);
                        peak_up_stroke(n,1) = a.LP.stats{idx(k),1}.peakUpStroke(1);
                        peak_down_stroke(n,1) = a.LP.stats{idx(k),1}.peakDownStroke(1);
                        peak_stroke_ratio(n,1) = a.LP.stats{idx(k),1}.peakStrokeRatio(1);
                        trough(n,1) = a.LP.stats{idx(k),1}.trough(1);
                        fastTrough(n,1) = a.LP.stats{idx(k),1}.fastTroughDur(1);
                        slowTrough(n,1) = a.LP.stats{idx(k),1}.slowTroughDur(1);
                        if size(a.LP.stats{idx(k),1}.waves,1) > 1
                            wf(n,:) = mean(a.LP.stats{idx(k),1}.waves);
                        else
                            wf(n,:) = a.LP.stats{idx(k),1}.waves;
                        end
                        break
                    else
                        rheobaseLP(n,1) = NaN;
                        break
                    end
                else
                    rheobaseLP(n,1) = NaN;
                end
            end
            maxFiringRate(n,1) = max(firing_rate_s(n,:));

            clear B I idx amp temp k
        end
    end
end

T = table(ID,specimen,struct,transline,reporterStatus,dendrite_type,pyramidalID,...
    resistance,time_constant,time_constant2,rebound_slope,nb_rebound_sp,...
    sag,sag_ratio,input_current_s,firing_rate_s,rheobaseLP,delay,burst,latency,cv_ISI,...
    adaptation1,adaptation2,peak_adaptation1,peak_adaptation2,peakLP,thresholdLP,...
    half_width_threshold_peak,half_width_peak_trough,height_threshold_peak,height_peak_trough,...
    peak_up_stroke,peak_down_stroke,peak_stroke_ratio,trough,fastTrough,slowTrough,maxFiringRate);
writetable(T,'cell type details (22 02 2020).xlsx','Sheet','Sheet1','WriteRowNames',true)
T = table(wf);
writetable(T,'cell type details (22 02 2020).xlsx','Sheet','Sheet2','WriteRowNames',true)
system('taskkill /F /IM EXCEL.EXE');

save('cell type details (22 02 2020).mat','ID','specimen','struct','transline','reporterStatus','dendrite_type','pyramidalID',...
    'resistance','time_constant','time_constant2','rebound_slope','nb_rebound_sp',...
    'sag','sag_ratio','input_current_s','firing_rate_s','rheobaseLP','delay','burst','latency','cv_ISI',...
    'adaptation1','adaptation2','peak_adaptation1','peak_adaptation2','peakLP','thresholdLP',...
    'half_width_threshold_peak','half_width_peak_trough','height_threshold_peak','height_peak_trough',...
    'peak_up_stroke','peak_down_stroke','peak_stroke_ratio','trough','fastTrough','slowTrough','maxFiringRate');

listnan = isnan(resistance);
resistance(listnan) = [];
AIr(listnan) = [];
listnan = isnan(time_constant);
time_constant(listnan) = [];
AItau(listnan) = [];
listnan = isnan(rheobaseLP);
rheobaseLP(listnan) = [];
AIrheobase(listnan) = [];

figure('Position',[50 50 800 250]); set(gcf,'color','w');
subplot(1,3,1)
scatter(AIr,resistance)
xlabel('resistance (AI)')
ylabel('resistance (JMT)')
[r_RI,p_Ri] = corr(AIr,resistance);
subplot(1,3,2)
scatter(AItau,time_constant)
xlabel('time constant (AI)')
ylabel('time constant (JMT)')
[r_tau,p_tau] = corr(AItau,time_constant);
subplot(1,3,3)
scatter(AIrheobase,rheobaseLP)
xlabel('rheobase (AI)')
ylabel('rheobase (JMT)')
[r_rheobase,p_rheobase] = corr(AIrheobase,rheobaseLP);
