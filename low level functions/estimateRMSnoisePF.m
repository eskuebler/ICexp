%{
estimateRMSnoise
%}

% vectors to estimate noise
vec_pre = double(a.LP.V{1,k}(a.LP.stimOn-(500/a.LP.acquireRes)-1:a.LP.stimOn-1));
% t_pre = LP.stimOn-(500/LP.acquireRes)-1:LP.stimOn-1;
vec_post = double(a.LP.V{1,k}(a.LP.stimOff+(250/a.LP.acquireRes):a.LP.stimOff+(750/a.LP.acquireRes)-1));
% t_post = LP.stimOff+(250/LP.acquireRes):LP.stimOff+(750/LP.acquireRes)-1;

% figure('Position',[50 50 250 250]); set(gcf,'color','w');
% hold on
% plot(vec_pre)
% plot(vec_post)
% xlabel('time-steps')
% ylabel('voltage (mV)')
% axis tight
% ylim([-100 -30])
% legend({'pre-stim','post-stim'})
% export_fig(['D:\genpath\',cellID,' RMS noise vectors ',int2str(k)],'-pdf','-r100');
% close

% long-term noise
a.LP.restVPre(1,k) = mean(vec_pre);
rmse_pre = sqrt(mean((vec_pre - a.LP.restVPre(1,k)).^2));
a.LP.restVPost(1,k) = mean(vec_post);
rmse_post = sqrt(mean((vec_post - a.LP.restVPost(1,k)).^2));

% short-term noise
stWin = 1.5/a.LP.acquireRes;
winCount = 1;
for i = 1:stWin/2:length(vec_post)-stWin
    yhat = mean(vec_pre(1,i:i+stWin));
    rmse_pre_st(winCount) = sqrt(mean((vec_pre(1,i:i+stWin) - yhat).^2));
%     clear yhet
    yhat = mean(vec_post(1,i:i+stWin));
    rmse_post_st(winCount) = sqrt(mean((vec_post(1,i:i+stWin) - yhat).^2));
%     clear yhat
    winCount = winCount+1;
end

% figure('Position',[50 50 250 250]); set(gcf,'color','w');
% hold on
% plot(rmse_pre_st)
% plot(rmse_post_st)
rmse_pre_st = mean(rmse_pre_st);
rmse_post_st = mean(rmse_post_st);

% Vm difference pre-stimulus and post-stimulus
diffV_b_e = abs(a.LP.restVPre(1,k)-a.LP.restVPost(1,k)); % differnce between end and stim onset
