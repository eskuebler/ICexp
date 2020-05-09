%{
summaryICanalysis
%}

clear; close all; clc;

% QC waveform (count bad spikes, register by QC)
    % first normal then morphed
    % if bad spikes come after good spike re want to
    % remove 

%     ADD Acuire res
    
% free parameters (across sweep QC)
minNmaxThres = 10;
origStdThresMax = 3;
origStdThresMin = 0.7;
savefilename = 'cell type details (2020 05 04)';
savefilenameInd = 'cell type details (2020 05 04) sp qc';
savefilenameIndBin = 'cell type details (2020 05 04) sp qc binary';

% load data
load('dendrite_layer.mat')
layerID = ID; dendrite_typeMJ = dendrite_type; clear ID dendrite_type
load('pyramidal cell IDs.mat')
pyrID = ID; clear ID
load('cell_types_specimen_details.mat','donor__species','specimen__id',...
    'line_name','tag__dendrite_type','ef__ri','ef__tau',...
    'ef__threshold_i_long_square','cell_reporter_status',...
    'structure__acronym','structure__layer');

% list of cells
cellList = dir('D:\genpath\*.mat');

% initialize
input_current_s = zeros(length(cellList),25);
firing_rate_s = zeros(length(cellList),25);
removedListStd = []; rmvdStdCount = 1;
removedListMinMax = []; rmvdMMCount = 1;
spqcmat = zeros(length(cellList),60);

for n = 1:281%length(cellList)
    clc; disp(n)
    sweepIDcount = 1;
    
    if length(cellList(n).name) == 17
        cellID = cellList(n).name(1:length(cellList(n).name)-4);
    else
        cellID = cellList(n).name(10:length(cellList(n).name)-4);
    end

    
    % pull data for AIBS cell
    if length(cellID) == 9
        % match species ID # to get Allen parameters
        k  = find(specimen__id==str2num(cellID));

        % get Allen Institute parameters
        ID(n,:) = ([num2str(specimen__id(k)),'__AI']);
        specimen(n,1) = donor__species(k);
        struct(n,1) = structure__acronym(k);
        cortical_layer(n,1) = categorical(structure__layer(k));
        transline(n,1) = line_name(k);
        reporterStatus(n,1) = cell_reporter_status(k);
        dendrite_type(n,1) = tag__dendrite_type(k);
        AIr(n,1) = ef__ri(k);
        AItau(n,1) = ef__tau(k);
        AIrheobase(n,1) = round(ef__threshold_i_long_square(k));
        if sum(ismember(pyrID,str2num(ID(n,1:9))))
            pyramidalID(n,1) = 1;
        else
            pyramidalID(n,1) = 0;
        end
        access_resistance(n,1) = NaN;
        temperature(n,1) = 0;
        
    % pull data for PCTD cell
    else
        ID(n,:) = cellID;
        specimen(n,:) = categorical(cellstr('NHP'));
        struct(n,:) = categorical(cellstr('PFC'));
        transline(n,:) = categorical(cellstr('N/A'));
        reporterStatus(n,:) = categorical(cellstr('N/A'));
        AIr(n,1) = NaN;
        AItau(n,1) = NaN;
        AIrheobase(n,1) = NaN;
        pyramidalID(n,1) = NaN;
        if sum(ismember(layerID,cellID))
            k = find(ismember(layerID,cellID)==1);
            if AccesResistance(k,1) ~= 'N/A'
                access_resistance(n,1) = ...
                    str2num(char(AccesResistance(k,1)));
            end
            dendrite_type(n,1) = dendrite_typeMJ(k,1);
            cortical_layer(n,1) = Layer(k,1);
            if Temperature(k,1) ~= 'N/A'
                temperature(n,1) = str2num(char(Temperature(k,1)));
            end
        else
            access_resistance(n,1) = NaN;
            dendrite_type(n,1) = categorical(cellstr('N/A'));
            cortical_layer(n,1) = categorical(cellstr('N/A'));
            temperature(n,1) = 0;
        end
    end

    % load our analysis
    load(['D:\genpath\',cellList(n).name]);
    
    % qc stuff
    qc_logic = zeros(1,6);

    if a.LP.fullStruct == 1
        qcID{n,1} = cellID;
        qc_V_vec(n,1:length(a.LP.rmp(1,:))) = round(a.LP.rmp(1,:),2);
        qc_V_vecDelta(n,1:length(a.LP.rmp(1,:))) = ...
            round(a.LP.rmp(1,1) - a.LP.rmp(1,:),2);
        spqcmatn = zeros(length(a.LP.sweepAmps),6);
        spqcmatnbinary = nan(20,300); spqcmatnbinaryid = nan(20,300);
        input_current_spqc = zeros(20,1);
        for k = 1:length(a.LP.sweepAmps)
            % quality control parameters
            qc_restVpre(n,k) = round(a.LP.stats{k,1}.qc.restVPre,2);
            qc_restVpost(n,k) = round(a.LP.stats{k,1}.qc.restVPost,2);
            qc_restVdiffpreNpost(n,k) = round(a.LP.stats{k,1}.qc.diffV_b_e,2);
            qc_rmse_pre_lt(n,k) = round(a.LP.stats{k,1}.qc.rmse_pre,2);
            qc_rmse_post_lt(n,k) = round(a.LP.stats{k,1}.qc.rmse_post,2);
            qc_rmse_pre_st(n,k) = round(a.LP.stats{k,1}.qc.rmse_pre_st,2);
            qc_rmse_post_st(n,k) = round(a.LP.stats{k,1}.qc.rmse_post_st,2);
            qc_logic = qc_logic+a.LP.stats{k,1}.qc.logicVec;
            vec = find(a.LP.stats{k,1}.qc.logicVec==1);
            if length(vec) == 1
                qc_class_mat(n,k) = vec;
            elseif length(vec) == 2
                result = strcat(num2str(vec(1)),num2str(vec(2)));
                result = str2num(result);
                qc_class_mat(n,k) = result;
            elseif length(vec) == 3
                result = strcat(num2str(vec(1)),num2str(vec(2)),...
                    num2str(vec(3)));
                result = str2num(result);
                qc_class_mat(n,k) = result; 
            elseif length(vec) == 4
                result = strcat(num2str(vec(1)),num2str(vec(2)),...
                    num2str(vec(3)),num2str(vec(4)));
                result = str2num(result);
                qc_class_mat(n,k) = result; 
            elseif length(vec) == 5
                result = strcat(num2str(vec(1)),num2str(vec(2)),...
                    num2str(vec(3)),num2str(vec(4)),num2str(vec(5)));
                result = str2num(result);
                qc_class_mat(n,k) = result; 
            elseif length(vec) == 6
                result = strcat(num2str(vec(1)),num2str(vec(2)),...
                    num2str(vec(3)),num2str(vec(4)),num2str(vec(5)),...
                    num2str(vec(6)));
                result = str2num(result);
                qc_class_mat(n,k) = result; 
            end
            if isfield(a.LP,'stats') && k<=length(a.LP.stats) && ...
                    ~isempty(a.LP.stats{k,1}) && ...
                    sum(a.LP.stats{k,1}.qc.logicVec) == 0
                sweepID(n,k) = k;
                sweepBinary(n,k) = 1;
                sweepBinaryOrig(1,k) = 1;
            else
                sweepID(n,k) = 0;
                sweepBinary(n,k) = 0;
                sweepBinaryOrig(1,k) = 0;
            end
            processSpQC
        end
        qc_logic_mat(n,1:6) = qc_logic;
        processBwSweepsQC                                                   % across sweep QC
    end
    if size(qc_class_mat,1)~=n
        qc_class_mat(n,:) = 0;
    end
    
    if exist('input_current_spqc','var')
        
        % save QC for cell
        % good spikes vs bad spikes
        % ratio 0.3
        
        T = table(input_current_spqc,spqcmatnbinary);
        writetable(T,[savefilenameInd,' ',cellID,'.xlsx'],'Sheet',...
            'Sheet1','WriteRowNames',true)
        T = table(input_current_spqc,spqcmatnbinaryid);
        writetable(T,[savefilenameInd,' ',cellID,'.xlsx'],'Sheet',...
            'Sheet2','WriteRowNames',true)
        clear input_current_spqc spqcmatnbinary spqcmatnbinaryid
    end
    
    if size(sweepBinary,1)==n && sum(sweepBinary(n,:))>0
                
%         if n == 294
%            'YO' 
%         end
        % subthreshold summary parameters
        acquireRes(n,1) = double(a.LP.acquireRes);
        resistance(n,1) = round(double(a.LP.subSummary.resistance),2);
        time_constant(n,1) = round(double(a.LP.subSummary.tauMin),2);
        time_constant2(n,1) = round(double(a.LP.subSummary.tauSS),2);
        if a.LP.fullStruct == 1
            k = find(ismember(sweepID(n,:),...
                find(round(double(a.LP.sweepAmps)) == -90))==1);                         % find -90 pA input
            if length(k)>1
                k = k(1);
            end
            if ~isempty(k)
                getSubthresholdStats
            else
                k = find(ismember(sweepID(n,:),...
                    find(round(double(a.LP.sweepAmps)) == -70))==1);
                if length(k)>1
                    k = k(1);
                end
                if ~isempty(k)
                    getSubthresholdStats
                else
                    k = find(ismember(sweepID(n,:),...
                        find(round(double(a.LP.sweepAmps)) == -50))==1);
                    if length(k)>1
                        k = k(1);
                    end
                    if ~isempty(k)
                        getSubthresholdStats
                    else
                        subamp(n,1) = NaN;
                        submin(n,1) = NaN;
                        rebound_slope(n,1) = NaN;
                        rebound_depolarization(n,1) = NaN;
                        nb_rebound_sp(n,1) = 0;
                        sag(n,1) = NaN;
                        steadystate(n,1) = NaN;
                        sag_ratio(n,1) = NaN;
                    end
                end
            end

            % find rheobase and parameters of first spike
            [B,I] = sort(round(double(a.LP.sweepAmps(1:length(a.LP.stats)))));
            input_current_all(n,1:length(B)) = B;
            int_vec = find(B>0);
            temp = int_vec(find(ismember(int_vec,sweepID(n,:))==1));
            idx = I(temp);
            amp = B(temp);
            input_current_s(n,1:length(amp)) = round(double(a.LP.sweepAmps(idx)));

            for k = 1:length(idx)
                if isfield(a.LP.stats{idx(k),1},'spTimes') && ...
                        sum(~isnan(a.LP.stats{idx(k),1}.spTimes))>0
                    % spike train parameters
                    firing_rate_s(n,k) = a.LP.stats{idx(k),1}.meanFR1000;
                    train_delay(n,k) = round(double(a.LP.stats{idx(k),1}.delay),2);
                    train_burst(n,k) = round(double(a.LP.stats{idx(k),1}.burst),2);
                    train_latency(n,k) = round(double(a.LP.stats{idx(k),1}.latency),2);
                    train_cv_ISI(n,k) = round(double(a.LP.stats{idx(k),1}.cvISI),2);
                    train_adaptation1(n,k) = round(double(a.LP.stats{idx(k),1}.adaptIndex),2);
                    train_adaptation2(n,k) = round(double(a.LP.stats{idx(k),1}.adaptIndex2),2);
                    train_peak_adaptation1(n,k) = round(double(a.LP.stats{idx(k),1}.peakAdapt),2);
                    train_peak_adaptation2(n,k) = round(double(a.LP.stats{idx(k),1}.peakAdapt2),2);
                end
            end

            for k = 1:length(idx)
                if isfield(a.LP.stats{idx(k),1},'spTimes') && ...
                        sum(~isnan(a.LP.stats{idx(k),1}.spTimes))>0
                    rheobaseLP(n,1) = amp(k);
                    delay(n,1) = round(double(a.LP.stats{idx(k),1}.delay),2);
                    burst(n,1) = round(double(a.LP.stats{idx(k),1}.burst),2);
                    latency(n,1) = round(double(a.LP.stats{idx(k),1}.latency),2);
                    cv_ISI(n,1) = round(double(a.LP.stats{idx(k),1}.cvISI),2);
                    adaptation1(n,1) = round(double(a.LP.stats{idx(k),1}.adaptIndex),2);
                    adaptation2(n,1) = round(double(a.LP.stats{idx(k),1}.adaptIndex2),2);
                    peak_adaptation1(n,1) = round(double(a.LP.stats{idx(k),1}.peakAdapt),2);
                    peak_adaptation2(n,1) = round(double(a.LP.stats{idx(k),1}.peakAdapt2),2);
                    peakLP(n,1) = round(double(a.LP.stats{idx(k),1}.peak(1)),2);
                    thresholdLP(n,1) = round(double(a.LP.stats{idx(k),1}.thresholdRef(1)),2);
                    half_width_threshold_peak(n,1) = round(double(a.LP.stats{idx(k),1}.fullWidthTP(1)),2);
                    half_width_peak_trough(n,1) = round(double(a.LP.stats{idx(k),1}.fullWidthPT(1)),2);
                    height_threshold_peak(n,1) = round(double(a.LP.stats{idx(k),1}.heightTP(1)),2);
                    height_peak_trough(n,1) = round(double(a.LP.stats{idx(k),1}.heightPT(1)),2);
                    peak_up_stroke(n,1) = round(double(a.LP.stats{idx(k),1}.peakUpStroke(1)),2);
                    peak_down_stroke(n,1) = round(double(a.LP.stats{idx(k),1}.peakDownStroke(1)),2);
                    peak_stroke_ratio(n,1) = round(double(a.LP.stats{idx(k),1}.peakStrokeRatio(1)),2);
                    trough(n,1) = round(double(a.LP.stats{idx(k),1}.trough(1)),2);
                    fastTrough(n,1) = a.LP.stats{idx(k),1}.fastTroughDur(1);
                    slowTrough(n,1) = a.LP.stats{idx(k),1}.slowTroughDur(1);
                    if size(a.LP.stats{idx(k),1}.waves,1) > 1
                        wf(n,:) = round(mean(a.LP.stats{idx(k),1}.waves),2);
                    else
                        wf(n,:) = round(a.LP.stats{idx(k),1}.waves,2);
                    end
                    break
                else
                    rheobaseLP(n,1) = NaN;
                end
            end
            maxFiringRate(n,1) = max(firing_rate_s(n,:));

            clear B I idx amp temp k
        end
    end
end

% ind = find(diffMinMaxV~=0);
% figure('Position',[50 50 300 250]); set(gcf,'color','w');
% histogram(diffMinMaxV(ind),40,'FaceColor','k');
% line([minNmaxThres,minNmaxThres],[0,30], ...
%             'color','r','linewidth',1,'linestyle','--');
% xlabel('diff b/w min and max V')
% ylabel('probability')
% axis tight
% box off
% % export_fig(['qc mean diff min max resting V'],'-pdf','-r100');
% % close
% 
% ind = find(meanOrigV~=0);
% figure('Position',[50 50 300 250]); set(gcf,'color','w');
% histogram(meanOrigV(ind),40,'FaceColor','k');
% xlabel('mean in V across sweep')
% ylabel('probability')
% axis tight
% box off
% % export_fig(['qc mean resting V'],'-pdf','-r100');
% % close
% 
% ind = find(stdOrigV~=0);
% figure('Position',[50 50 300 250]); set(gcf,'color','w');
% hold on
% histogram(stdOrigV(ind),40,'FaceColor','k');
% line([origStdThresMax,origStdThresMax],[0,20], ...
%             'color','r','linewidth',1,'linestyle','--');
% line([origStdThresMin,origStdThresMin],[0,20], ...
%             'color','g','linewidth',1,'linestyle','--');
% xlabel('std in V across sweep')
% ylabel('probability')
% axis tight
% box off
% % export_fig(['qc std resting V'],'-pdf','-r100');
% close

% QCanalysis

% all parameters
T = table(ID,specimen,struct,cortical_layer,transline,reporterStatus,...
    dendrite_type,pyramidalID,temperature,resistance,access_resistance,...
    time_constant,time_constant2,subamp,submin,rebound_slope,nb_rebound_sp,sag,...
    steadystate,sag_ratio,rheobaseLP,delay,burst,latency,cv_ISI,adaptation1,...
    adaptation2,peak_adaptation1,peak_adaptation2,peakLP,thresholdLP,...
    half_width_threshold_peak,half_width_peak_trough,...
    height_threshold_peak,height_peak_trough,...
    peak_up_stroke,peak_down_stroke,peak_stroke_ratio,trough,...
    fastTrough,slowTrough,maxFiringRate);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet1',...
        'WriteRowNames',true)
% sweeps that pass qc
T = table(ID,qc_class_mat);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet2',...
        'WriteRowNames',true)
% spike trains all sweeps
T = table(ID,input_current_s,firing_rate_s,train_delay,train_burst,...
    train_latency,train_cv_ISI,train_adaptation1,train_adaptation2,...
    train_peak_adaptation1,train_peak_adaptation2);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet3',...
        'WriteRowNames',true)
% qc parameters for all sweeps
T = table(ID,meanOrigV,diffMinMaxV,stdOrigV,qc_V_vec,qc_restVpre,qc_restVpost,qc_restVdiffpreNpost,...
    qc_rmse_pre_lt,qc_rmse_post_lt,qc_rmse_pre_st,qc_rmse_post_st);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet4',...
        'WriteRowNames',true)
% waveforms
T = table(ID,wf);
    writetable(T,[savefilename,'.xlsx'],'Sheet','Sheet5',...
        'WriteRowNames',true)
% close excel
system('taskkill /F /IM EXCEL.EXE');

% plot our analysis versus Allen Institute
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
