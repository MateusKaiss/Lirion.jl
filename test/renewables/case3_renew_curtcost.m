% 3-bus single-period case with gen 2 (renewable) massively over-sized
% (pmax = 1000 MW >> demand) and a high curtailment_cost (50 /MW).
% Verifies that the curtailment_cost term in the objective incentivizes
% the optimizer to dispatch renewable energy rather than curtail it.

function mpc = case3_renew_curtcost
mpc.version = '2';
mpc.baseMVA = 100.0;

mpc.bus = [
    1  3   110.0  40.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
    2  2   110.0  40.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
    3  2    95.0  50.0  0.0  0.0  1  1.0  0.0  240.0  1  1.1  0.9;
];

mpc.gen = [
    1  1000.0  0.0  1000.0  -1000.0  1.0  100.0  1  2000.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'T'   0.0;
    2  1000.0  0.0  1000.0  -1000.0  1.0  100.0  1  1000.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'R'  50.0;
    3     0.0  0.0  1000.0  -1000.0  1.0  100.0  1     0.0  0.0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  'T'   0.0;
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
