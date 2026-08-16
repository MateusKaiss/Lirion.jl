% 3-bus T=2 case exercising the pmin / qmax / qmin columns of
% mpc.gen_timeseries. Gen 1 has a must-run floor of 300 MW in t=1
% (pmin_ts[1,1] = 300) and is free in t=2. Gen 2 has a reactive
% ceiling tightened to 50 MVAr in t=2 (qmax_ts[2,2] = 50).
% Used to verify that time-varying bounds are applied selectively
% only when strictly tighter than the static gen table.

function mpc = case3_gen_ts_bounds
mpc.version = '2';
mpc.baseMVA = 100.0;

mpc.bus = [
    1  3   110.0  40.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
    2  2   110.0  40.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
    3  2    95.0  50.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
];

mpc.gen = [
    1  100.0  0.0  1000.0  -1000.0  1.0  100.0  1  2000.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'T'  0.0;
    2  100.0  0.0  1000.0  -1000.0  1.0  100.0  1  2000.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'T'  0.0;
    3    0.0  0.0  1000.0  -1000.0  1.0  100.0  1     0.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'T'  0.0;
];

mpc.gencost = [
    2  0.0  0.0  3  0.110  5.0  0.0;
    2  0.0  0.0  3  0.085  1.2  0.0;
    2  0.0  0.0  3  0.000  0.0  0.0;
];

mpc.branch = [
    1  3  0.065  0.62  0.45  9000.0  9000.0  9000.0  0.0  0.0  1  -60.0  60.0;
    3  2  0.025  0.75  0.70  9000.0  9000.0  9000.0  0.0  0.0  1  -60.0  60.0;
    1  2  0.042  0.90  0.30  9000.0  9000.0  9000.0  0.0  0.0  1  -60.0  60.0;
];

mpc.bus_timeseries = [
    1  1  110.0  40.0;
    1  2  110.0  40.0;
    1  3   95.0  50.0;
    2  1  110.0  40.0;
    2  2  110.0  40.0;
    2  3   95.0  50.0;
];

%% gen_timeseries — columns: time, gen_id, pmax, pmin, qmax, qmin (NaN = use static).
%   Gen 1 must-run 300 MW at t=1 (pmin_ts overrides the static pmin=0).
%   Gen 2 reactive ceiling tightened to 50 MVAr at t=2 (qmax_ts overrides static 1000).
mpc.gen_timeseries = [
    1  1  NaN  300.0  NaN  NaN;
    2  2  NaN    NaN   50.0  NaN;
];
