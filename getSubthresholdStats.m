% collect subthreshold data

subamp(n,1) = round(double(a.LP.stats{k,1}.subSweepAmps),2);
submin(n,1) = round(double(a.LP.stats{k,1}.minV),2);
rebound_slope(n,1) = round(double(a.LP.stats{k,1}.reboundSlope),2);
rebound_depolarization(n,1) = round(double(a.LP.stats{k,1}.reboundDepolarization),2);
sag(n,1) = round(double(a.LP.stats{k,1}.sag),2);
steadystate(n,1) = round(double(a.LP.stats{k,1}.subSteadyState),2);
sag_ratio(n,1) = round(double(a.LP.stats{k,1}.sagRatio),2);
if sum(isnan(a.LP.stats{k,1}.reboundAPs))==0
    nb_rebound_sp(n,1) = length(a.LP.stats{k,1}.reboundAPs);
else
    nb_rebound_sp(n,1) = 0;
end