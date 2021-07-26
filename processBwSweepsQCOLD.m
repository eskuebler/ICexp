%{
processBwSweepsQC
%}

ind = find(sweepBinary(n,:)==1);                                            % index of sweeps that passed other QC
vec = qc_V_vec(n,ind);                                                      % resting voltages for QC'd sweeps
if length(vec) > 1                                                          % if one sweep don't analyze
    diffMinMaxV(n,1) = round(abs(min(vec)-max(vec)),2);                     % diff b/w min & max
    if diffMinMaxV(n,1) < minNmaxThres                                      % if diff b/w min & max > thres
        stdOrigV(n,1) = round(std(vec),2);                                  % s.d. original voltages
        meanOrigV(n,1) = round(mean(vec),2);                                % mean original voltages
        if stdOrigV(n,1) > origStdThresMax                                  % if s.d. > max threshold
            removedListStd{rmvdStdCount,1} = cellID;                        % add ID to s.d. remove cell list
            rmvdStdCount = rmvdStdCount + 1;
            sweepID(n,:) = 0;                                               % clear count QC vector
            sweepBinary(n,:) = 0;                                           % clear binary QC vector
            qc_class_mat(n,:) = 8;
        elseif stdOrigV(n,1) >= origStdThresMin && ...
                stdOrigV(n,1) <= origStdThresMax                            % if s.d. b/w both min and max thres
            outlierVec = find(abs(vec-meanOrigV(n,1)) > ...
                (stdOrigV(n,1)*1.75));                                      % find sweeps > 1.75 s.d. +/- mean
            sweepBinary(n,ind(outlierVec)) = 0;                             % remove failed sweeps from binary QC vector
            qc_class_mat(n,ind(outlierVec)) = 9; 
%             figure('Position',[50 50 300 250]); set(gcf,'color','w');       % generate figure
%             hold on
%             scatter(ind,vec,'k')
%             scatter(ind(outlierVec),vec(outlierVec),'r')
%             line([1,length(sweepBinaryOrig)],[meanOrigV(n,1),meanOrigV(n,1)], ...
%                     'color','b','linewidth',1,'linestyle','--');
%             xlabel('current input')
%             xticks(1:length(sweepBinaryOrig))
%             xticklabels({a.LP.sweepAmps(1:length(sweepBinaryOrig),1)'})
%             xtickangle(90)
%             ylabel('resting V (mV)')
%             axis tight
%             ylim([-90 -40])
%             export_fig(['D:\test\', ...
%                 cellList(n).name(1:length(cellList(n).name)-4), ...
%                 ' rmp w outliers'],'-pdf','-r100');                         % save figure
%             close                                                           % close figure
        elseif stdOrigV(n,1) < origStdThresMin                              % if s.d. < min thres
%             figure('Position',[50 50 300 250]); set(gcf,'color','w');       % generate figure
%             hold on
%             scatter(ind,vec,'k')
%             line([1,length(sweepBinaryOrig)],[meanOrigV(n,1),meanOrigV(n,1)], ...
%                     'color','b','linewidth',1,'linestyle','--');
%             xlabel('current input')
%             xticks(1:length(sweepBinaryOrig))
%             xticklabels({a.LP.sweepAmps(1:length(sweepBinaryOrig),1)'})
%             xtickangle(90)
%             ylabel('resting V (mV)')
%             axis tight
%             ylim([-90 -40])
%             export_fig(['D:\test\', ...
%                 cellList(n).name(1:length(cellList(n).name)-4), ...
%                 ' rmp w outliers'],'-pdf','-r100');                         % save figure
%             close                                                           % close figure
        end
    else
        removedListMinMax{rmvdMMCount,1} = cellID;                          % add to min/max remove cell list
        rmvdMMCount = rmvdMMCount + 1;
        sweepID(n,:) = 0;                                                   % clear count QC vector
        sweepBinary(n,:) = 0;                                               % clear binary QC vector
        qc_class_mat(n,:) = 7;
    end
    diffBwSwpQC = sum(sweepBinaryOrig)-sum(sweepBinary(n,:));               % diff b/w number of passed sweeps
    qc_logic_mat(n,7) = diffBwSwpQC;                                        % add diff to QC count matrix
%     qc_class_mat(n,:) = 7;
end
clear diffBwSwpQC ind vec outlierVec sweepBinaryOrig
