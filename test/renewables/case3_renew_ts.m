% 3-bus multi-period case (T=3) with gen 2 tagged as renewable and a
% time-varying availability profile. Gen 2's pmax_ts is 50/150/200 MW
% across periods 1/2/3 (e.g. ramp-up solar profile).
% Used to verify the constraint  p_gen[g,t] + curt[g,t] == pmax_ts[g,t].

function mpc = case3_renew_ts
mpc.version = '2';
mpc.baseMVA = 100.0;

mpc.bus = [
    1  3   110.0  40.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
    2  2   110.0  40.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
    3  2    95.0  50.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
];

mpc.gen = [
    1  1000.0  0.0  1000.0  -1000.0  1.0  100.0  1  2000.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'T'  0.0;
    2  1000.0  0.0  1000.0  -1000.0  1.0  100.0  1   200.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'R'  0.0;
    3     0.0  0.0  1000.0  -1000.0  1.0  100.0  1     0.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'T'  0.0;
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

%% bus timeseries (same demand in every period for test clarity)
mpc.bus_timeseries = [
    1  1  110.0  40.0;
    1  2  110.0  40.0;
    1  3   95.0  50.0;
    2  1  110.0  40.0;
    2  2  110.0  40.0;
    2  3   95.0  50.0;
    3  1  110.0  40.0;
    3  2  110.0  40.0;
    3  3   95.0  50.0;
];

%% gen availability profile: gen 2 (solar), ramp-up across periods
%    time   gen_id   pmax_available
mpc.gen_timeseries = [
    1  2   50.0;
    2  2  150.0;
    3  2  200.0;
];
