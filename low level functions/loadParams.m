function p = loadParams

% spike detection
p.thresholdV = -20;                              % detection threshold for V-based spikes
p.thresholdDVDT = 20;                          % detection threshold for dV/dt-based spikes

% swep-wise root mean square quality control parameters
p.RMSEst = .2;                               % maximum RMSE measure short term
p.RMSElt = 1.05;                                % maximum RMSE measure long term
p.maxDiffBwBeginEnd = 3;                      % maximum difference between beginning and end of sweep
p.maxDiffBwSweeps = 10;                         % maximum difference b/w sweeps
p.minimumRestingPot = -50;                     % minimum resting potential

% rebound slope and spike parameters
p.reboundWindow = 100;                         % window to find maximum rebound peak
p.reboundFitWindow = 150;                      % window from max rebound peak to fit / acquireRes
p.reboundSpWindow = 50;                        % window to look for rebound spikes (ms)

% target sampling parameters
p.sampleRT = 5e4;                              % sample rate we want
p.sampleRTdt = 1000/p.sampleRT;                % sample rate we want

% spike-wise quality control parameters
p.minDiffThreshold2Peak = 20;                  % max diff in V bw threshold and peak
p.maxDiffThreshold2PeakT = 2;                  % max diff in t bw threshold and peak
p.minDiffPeak2Trough = 30;                     % max diff in V bw peak and trough
p.maxDiffPeak2TroughT = 10;                    % max diff in t bw peak and trough
p.percentRheobaseHeight = .3;                 % APs must be X percent of Rheobase height
p.maxThreshold = -20;                          % above this value APs are eliminated (mV)
p.minTrough = -30;                             % above this value APs are eliminated (mV)