%{
QCanalysis
%}

p = loadParams;

figure('Position',[50 50 400 600]); set(gcf,'color','w');
subplot(3,4,1)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_restVpre(n,1) qc_restVpost(n,1)],'k','linewidth',0.25)
end
plot([1 2],[-50 -50],'r-.','linewidth',1)
xlabel('condition')
ylabel('voltage (mV)')
xticks(1:2)
xticklabels({'pre','post'})
subplot(3,4,2)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_restVpre(n,6) qc_restVpost(n,6)],'k','linewidth',0.25)
end
plot([1 2],[-50 -50],'r-.','linewidth',1)
xlabel('condition')
ylabel('voltage (mV)')
xticks(1:2)
xticklabels({'pre','post'})
subplot(3,4,3)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_restVpre(n,12) qc_restVpost(n,12)],'k','linewidth',0.25)
end
plot([1 2],[-50 -50],'r-.','linewidth',1)
xlabel('condition')
ylabel('voltage (mV)')
xticks(1:2)
xticklabels({'pre','post'})
subplot(3,4,4)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_restVpre(n,20) qc_restVpost(n,20)],'k','linewidth',0.25)
end
plot([1 2],[-50 -50],'r-.','linewidth',1)
xlabel('condition')
ylabel('voltage (mV)')
xticks(1:2)
xticklabels({'pre','post'})
subplot(3,4,5)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_rmse_pre_st(n,1) qc_rmse_post_st(n,1)],'k','linewidth',0.25)
end
plot([1 2],[p.RMSEst p.RMSEst],'r-.','linewidth',1)
xlabel('condition')
ylabel('RMS (short term)')
xticks(1:2)
ylim([0 0.7])
xticklabels({'pre','post'})
subplot(3,4,6)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_rmse_pre_st(n,6) qc_rmse_post_st(n,6)],'k','linewidth',0.25)
end
plot([1 2],[p.RMSEst p.RMSEst],'r-.','linewidth',1)
xlabel('condition')
ylabel('RMS (short term)')
xticks(1:2)
ylim([0 0.7])
xticklabels({'pre','post'})
subplot(3,4,7)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_rmse_pre_st(n,12) qc_rmse_post_st(n,12)],'k','linewidth',0.25)
end
plot([1 2],[p.RMSEst p.RMSEst],'r-.','linewidth',1)
xlabel('condition')
ylabel('RMS (short term)')
xticks(1:2)
ylim([0 0.7])
xticklabels({'pre','post'})
subplot(3,4,8)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_rmse_pre_st(n,20) qc_rmse_post_st(n,20)],'k','linewidth',0.25)
end
plot([1 2],[p.RMSEst p.RMSEst],'r-.','linewidth',1)
xlabel('condition')
ylabel('RMS (short term)')
xticks(1:2)
ylim([0 0.7])
xticklabels({'pre','post'})
subplot(3,4,9)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_rmse_pre_lt(n,1) qc_rmse_post_lt(n,1)],'k','linewidth',0.25)
end
plot([1 2],[p.RMSElt p.RMSElt],'r-.','linewidth',1)
xlabel('condition')
ylabel('RMS (long term)')
xticks(1:2)
ylim([0 14])
xticklabels({'pre','post'})
subplot(3,4,10)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_rmse_pre_lt(n,6) qc_rmse_post_lt(n,6)],'k','linewidth',0.25)
end
plot([1 2],[p.RMSElt p.RMSElt],'r-.','linewidth',1)
xlabel('condition')
ylabel('RMS (long term)')
xticks(1:2)
ylim([0 14])
xticklabels({'pre','post'})
subplot(3,4,11)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_rmse_pre_lt(n,12) qc_rmse_post_lt(n,12)],'k','linewidth',0.25)
end
plot([1 2],[p.RMSElt p.RMSElt],'r-.','linewidth',1)
xlabel('condition')
ylabel('RMS (long term)')
xticks(1:2)
ylim([0 14])
xticklabels({'pre','post'})
subplot(3,4,12)
hold on
for n = 1:281%length(cellList)
    plot([1 2],[qc_rmse_pre_lt(n,20) qc_rmse_post_lt(n,20)],'k','linewidth',0.25)
end
plot([1 2],[p.RMSElt p.RMSElt],'r-.','linewidth',1)
xlabel('condition')
ylabel('RMS (long term)')
xticks(1:2)
ylim([0 14])
xticklabels({'pre','post'})
% close

figure('Position',[50 50 900 900]); set(gcf,'color','w');
ind = qc_rmse_pre_lt(:,1)==0;
subplot(5,4,1)
histogram(qc_rmse_pre_lt(~ind,1),40,'FaceColor','k','Normalization','probability');
line([p.RMSElt,p.RMSElt],[0,0.4], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('lt RMS pre-stim')
    axis tight
    box off
ind = qc_rmse_pre_lt(:,6)==0;
subplot(5,4,2)
histogram(qc_rmse_pre_lt(~ind,6),40,'FaceColor','k','Normalization','probability');
line([p.RMSElt,p.RMSElt],[0,0.4], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('lt RMS pre-stim')
    axis tight
    box off
ind = qc_rmse_pre_lt(:,12)==0;
subplot(5,4,3)
histogram(qc_rmse_pre_lt(~ind,12),40,'FaceColor','k','Normalization','probability');
line([p.RMSElt,p.RMSElt],[0,0.4], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('lt RMS pre-stim')
    axis tight
    box off
ind = qc_rmse_pre_lt(:,20)==0;
subplot(5,4,4)
histogram(qc_rmse_pre_lt(~ind,20),40,'FaceColor','k','Normalization','probability');
line([p.RMSElt,p.RMSElt],[0,0.4], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('lt RMS pre-stim')
    axis tight
    box off
ind = qc_rmse_post_lt(:,1)==0;
subplot(5,4,5)
histogram(qc_rmse_post_lt(~ind,1),40,'FaceColor','k','Normalization','probability');
line([p.RMSElt,p.RMSElt],[0,0.4], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('lt RMS post-stim')
    axis tight
    box off
ind = qc_rmse_post_lt(:,6)==0;
subplot(5,4,6)
histogram(qc_rmse_post_lt(~ind,6),40,'FaceColor','k','Normalization','probability');
line([p.RMSElt,p.RMSElt],[0,0.4], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('lt RMS post-stim')
    axis tight
    box off
ind = qc_rmse_post_lt(:,12)==0;
subplot(5,4,7)
histogram(qc_rmse_post_lt(~ind,12),40,'FaceColor','k','Normalization','probability');
line([p.RMSElt,p.RMSElt],[0,0.4], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('lt RMS post-stim')
    axis tight
    box off
ind = qc_rmse_post_lt(:,20)==0;
subplot(5,4,8)
histogram(qc_rmse_post_lt(~ind,20),40,'FaceColor','k','Normalization','probability');
line([p.RMSElt,p.RMSElt],[0,0.4], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('lt RMS post-stim')
    axis tight
    box off
ind = qc_rmse_pre_st(:,1)==0;
subplot(5,4,9)
histogram(qc_rmse_pre_st(~ind,1),40,'FaceColor','k','Normalization','probability');
line([p.RMSEst,p.RMSEst],[0,0.3], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('short-term RMS')
    axis tight
    box off
ind = qc_rmse_pre_st(:,6)==0;
subplot(5,4,10)
histogram(qc_rmse_pre_st(~ind,6),40,'FaceColor','k','Normalization','probability');
line([p.RMSEst,p.RMSEst],[0,0.3], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('short-term RMS')
    axis tight
    box off
ind = qc_rmse_pre_st(:,12)==0;
subplot(5,4,11)
histogram(qc_rmse_pre_st(~ind,12),40,'FaceColor','k','Normalization','probability');
line([p.RMSEst,p.RMSEst],[0,0.3], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('short-term RMS')
    axis tight
    box off
ind = qc_rmse_pre_st(:,20)==0;
subplot(5,4,12)
histogram(qc_rmse_pre_st(~ind,20),40,'FaceColor','k','Normalization','probability');
line([p.RMSEst,p.RMSEst],[0,0.3], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('short-term RMS')
    axis tight
    box off
ind = qc_restVpre(:,1)==0;
subplot(5,4,13)
histogram(qc_restVpre(~ind,1),40,'FaceColor','k','Normalization','probability');
line([p.minimumRestingPot,p.minimumRestingPot],[0,0.05], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('V rest pre')
    axis tight
    box off
ind = qc_restVpre(:,6)==0;
subplot(5,4,14)
histogram(qc_restVpre(~ind,6),40,'FaceColor','k','Normalization','probability');
line([p.minimumRestingPot,p.minimumRestingPot],[0,0.05], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('V rest pre')
    axis tight
    box off
ind = qc_restVpre(:,12)==0;
subplot(5,4,15)
histogram(qc_restVpre(~ind,12),40,'FaceColor','k','Normalization','probability');
line([p.minimumRestingPot,p.minimumRestingPot],[0,0.05], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('V rest pre')
    axis tight
    box off
ind = qc_restVpre(:,20)==0;
subplot(5,4,16)
histogram(qc_restVpre(~ind,20),40,'FaceColor','k','Normalization','probability');
line([p.minimumRestingPot,p.minimumRestingPot],[0,0.05], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('V rest pre')
    axis tight
    box off
ind = qc_restVdiffpreNpost(:,1)==0;
subplot(5,4,17)
histogram(qc_restVdiffpreNpost(~ind,1),40,'FaceColor','k','Normalization','probability');
line([p.maxDiffBwBeginEnd,p.maxDiffBwBeginEnd],[0,0.2], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('delta V (abs(pre-post))')
    axis tight
    box off
ind = qc_restVdiffpreNpost(:,6)==0;
subplot(5,4,18)
histogram(qc_restVdiffpreNpost(~ind,6),40,'FaceColor','k','Normalization','probability');
line([p.maxDiffBwBeginEnd,p.maxDiffBwBeginEnd],[0,0.2], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('delta V (abs(pre-post))')
    axis tight
    box off
ind = qc_restVdiffpreNpost(:,12)==0;
subplot(5,4,19)
histogram(qc_restVdiffpreNpost(~ind,12),40,'FaceColor','k','Normalization','probability');
line([p.maxDiffBwBeginEnd,p.maxDiffBwBeginEnd],[0,0.2], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('delta V (abs(pre-post))')
    axis tight
    box off
ind = qc_restVdiffpreNpost(:,20)==0;
subplot(5,4,20)
histogram(qc_restVdiffpreNpost(~ind,20),40,'FaceColor','k','Normalization','probability');
line([p.maxDiffBwBeginEnd,p.maxDiffBwBeginEnd],[0,0.2], ...
            'color','r','linewidth',1,'linestyle','--');
    ylabel('probability')
    xlabel('delta V (abs(pre-post))')
    axis tight
    box off
% close

figure('Position',[50 50 900 200]); set(gcf,'color','w');
subplot(1,4,1)
scatter(qc_restVpre(1:281,1),qc_restVpost(1:281,1),10,'filled','k')
xlabel('V_r_e_s_t pre')
ylabel('V_r_e_s_t post')
xlim([-100 -20])
ylim([-100 -20])
box off
subplot(1,4,2)
scatter(qc_restVpre(1:281,6),qc_restVpost(1:281,6),10,'filled','k')
xlabel('V_r_e_s_t pre')
ylabel('V_r_e_s_t post')
xlim([-100 -20])
ylim([-100 -20])
box off
subplot(1,4,3)
scatter(qc_restVpre(1:281,12),qc_restVpost(1:281,12),10,'filled','k')
xlabel('V_r_e_s_t pre')
ylabel('V_r_e_s_t post')
xlim([-100 -20])
ylim([-100 -20])
box off
subplot(1,4,4)
scatter(qc_restVpre(1:281,20),qc_restVpost(1:281,20),10,'filled','k')
xlabel('V_r_e_s_t pre')
ylabel('V_r_e_s_t post')
xlim([-100 -20])
ylim([-100 -20])
box off

figure('Position',[0 0 300 1000]); set(gcf,'color','w');
imagesc(qc_logic_mat(1:281,:))
xlabel('QC criteria')
xticks(1:7)
xticklabels({'stRMS_p_r_e','stRMS_p_o_s_t',...
    'ltRMS_p_r_e','ltRMS_p_o_s_t',...
    'V_p_r_e-V_p_o_s_t<10','V_p_r_e>-50','b/w sweeps'})
xtickangle(45)
yticks(1:4:281)
ylabel('cell')
colorbar
colormap('gray');
box off

figure('Position',[50 50 750 200]); set(gcf,'color','w');
subplot(1,3,1)
hold on
for n = 1:281
    scatter(qc_rmse_pre_st(n,1),qc_rmse_post_st(n,1),3,'k')
end
plot([0.2 0.2],[0 0.2],'r','linewidth',0.25)
plot([0 0.2],[0.2 0.2],'r','linewidth',0.25)
xlabel('short term RMS pre')
ylabel('short term RMS post')
xlim([0 0.4])
ylim([0 0.4])
subplot(1,3,2)
hold on
for n = 1:281
    scatter(qc_rmse_pre_lt(n,1),qc_rmse_post_lt(n,1),3,'k')
end
plot([0.75 0.75],[0 0.75],'r','linewidth',0.25)
plot([0 0.75],[0.75 0.75],'r','linewidth',0.25)
xlabel('long term RMS pre')
ylabel('long term RMS post')
xlim([0 2])
ylim([0 2])
subplot(1,3,3)
hold on
for n = 1:281
    scatter(qc_rmse_pre_st(n,1),qc_rmse_pre_lt(n,1),3,'k')
end
plot([0.2 0.2],[0 0.5],'r','linewidth',0.25)
plot([0 0.2],[0.5 0.5],'r','linewidth',0.25)
xlabel('short term RMS')
ylabel('long term RMS')
xlim([0 0.25])
ylim([0 6.5])

figure('Position',[50 50 300 250]); set(gcf,'color','w');
plot(qc_V_vecDelta(1:281,1:20)')
xlabel('sweep #')
ylabel('diff(sweep(1),sweep(n))')
axis tight
box off

close all


classes = unique(qc_class_mat(1:281,:));
classes = classes(2:end);
qc_mat_classes = zeros(281,length(classes));
for n = 1:size(qc_class_mat(1:281,1))
    for k = 1:length(classes)
         qc_mat_classes(n,k) = qc_mat_classes(n,k)+length(find(qc_class_mat(n,:)==classes(k)));
    end
%     figure('Position',[50 50 600 250]); set(gcf,'color','w');
%     scatter(1:length(classes),qc_mat_classes(n,:),'k','filled')
%     xlabel('QC criteria combination')
%     ylabel('removal count')
%     xticks(1:length(classes))
%     xticklabels({classes})
%     xtickangle(45)
%     axis tight
%     ylim([0 25])
%     export_fig(['D:\test\', ...
%         cellList(n).name(1:length(cellList(n).name)-4), ...
%         ' qc removal counts'],'-pdf','-r100');
%     close
end

figure('Position',[50 50 600 800]); set(gcf,'color','w');
imagesc(qc_mat_classes)
xlabel('QC criteria combination')
ylabel('neuron')
colormap('gray')
xticks(1:length(classes))
xticklabels({classes})
xtickangle(45)
colorbar
box off

figure('Position',[50 50 600 250]); set(gcf,'color','w');
bar(sum(qc_mat_classes))
xlabel('QC criteria combination')
ylabel('count')
xticks(1:length(classes))
xticklabels({classes})
xtickangle(45)
axis tight
box off
