function mpc = case4_3zone_star
%CASE4_3ZONE_STAR  Tiny 4-bus star to exercise N-way SplitBus.
%  Bus 1 (zone 1) is the hub: branch 1-2 crosses into zone 2 and branch 1-3
%  crosses into zone 3, so under SplitBus bus 1 is a boundary bus shared by
%  THREE zones (its own + zone 2 + zone 3). Cheap generation lives only at the
%  hub, so the neighbor zones must import power across the junction -- this
%  forces a non-trivial net-injection (Sigma p_z = 0) consensus over 3 copies.

mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1.0	0	230	1	1.1	0.9;
	4	1	30	0	0	0	1	1.0	0	230	1	1.1	0.9;
	2	2	50	0	0	0	2	1.0	0	230	2	1.1	0.9;
	3	2	40	0	0	0	3	1.0	0	230	3	1.1	0.9;
];

%% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
mpc.gen = [
	1	0	0	200	-200	1.0	100	1	300	0;
	2	0	0	100	-100	1.0	100	1	100	0;
	3	0	0	100	-100	1.0	100	1	100	0;
];

%% fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax
mpc.branch = [
	1	4	0.005	0.1	0	300	300	300	0	0	1	-360	360;
	1	2	0.005	0.1	0	300	300	300	0	0	1	-360	360;
	1	3	0.005	0.1	0	300	300	300	0	0	1	-360	360;
];

%% model startup shutdown ncost c2 c1 c0
mpc.gencost = [
	2	0	0	3	0	10	0;
	2	0	0	3	0	100	0;
	2	0	0	3	0	100	0;
];
