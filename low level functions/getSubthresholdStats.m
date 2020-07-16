% collect subthreshold data

IC.subamp(n,1) = round(double(a.LP.stats{k,1}.subSweepAmps),2);
IC.submin(n,1) = round(double(a.LP.stats{k,1}.minV),2);
IC.rebound_slope(n,1) = round(double(a.LP.stats{k,1}.reboundSlope),2);
IC.rebound_depolarization(n,1) = round(double(a.LP.stats{k,1}.reboundDepolarization),2);
IC.sag(n,1) = round(double(a.LP.stats{k,1}.sag),2);
IC.steadystate(n,1) = round(double(a.LP.stats{k,1}.subSteadyState),2);
IC.sag_ratio(n,1) = round(double(a.LP.stats{k,1}.sagRatio),2);
IC.Vrest_sag_sweep(n,1) = double(a.LP.stats{k, 1}.qc.restVPre);  % Michelle change get Vm prestim of sweep for sag ratio
if sum(isnan(a.LP.stats{k,1}.reboundAPs))==0
    IC.nb_rebound_sp(n,1) = length(a.LP.stats{k,1}.reboundAPs);
else
    IC.nb_rebound_sp(n,1) = 0;
end