
shift = 0.1;

% 1 second window
tsf = 1;
c = 1;
for ts = 0:shift:tEnd-tsf
    fR(c,6) = sum(sp>ts & sp<(ts+tsf));
    c = c + 1;
end

% 750 ms window
tsf = 0.75;
c = 1;
for ts = 0:shift:tEnd-tsf
    fR(c,5) = sum(sp>ts & sp<(ts+tsf))/tsf;
    c = c + 1;
end

% 500 ms window
tsf = 0.5;
c = 1;
for ts = 0:shift:tEnd-tsf
    fR(c,4) = sum(sp>ts & sp<(ts+tsf))/tsf;
    c = c + 1;
end

% 250 ms window
tsf = 0.25;
c = 1;
for ts = 0:shift:tEnd-tsf
    fR(c,3) = sum(sp>ts & sp<(ts+tsf))/tsf;
    c = c + 1;
end

% 100 ms window
tsf = 0.1;
c = 1;
for ts = 0:shift:tEnd-tsf
    fR(c,2) = sum(sp>ts & sp<(ts+tsf))/tsf;
    c = c + 1;
end

% 50 ms window
tsf = 0.05;
c = 1;
for ts = 0:shift:tEnd-tsf
    fR(c,1) = sum(sp>ts & sp<(ts+tsf))/tsf;
    c = c + 1;
end