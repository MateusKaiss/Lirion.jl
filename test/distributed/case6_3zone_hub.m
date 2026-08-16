function mpc = case6_3zone_hub
%CASE6_3ZONE_HUB  6-bus, 3-zone case to exercise N-way SplitBus on SOCWR/AC/DC.
%  Bus 1 (zone 1) is the hub: branch 1-3 crosses into zone 2 and branch 1-5
%  crosses into zone 3, so under SplitBus bus 1 is a boundary bus shared by
%  THREE zones. Each zone owns two buses and at least one generator; the hub
%  carries the cheap unit so the neighbour zones import across the junction,
%  forcing a non-trivial net-injection (Sigma p_z = 0) consensus over 3 copies.
%  Impedances and reactive headroom are kept moderate so the SOCWR relaxation
%  is well-conditioned (Clarabel-friendly).

mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
mpc.bus = [
	1	3	20	10	0	0	1	1.0	0	230	1	1.1	0.9;
	2	1	40	20	0	0	1	1.0	0	230	1	1.1	0.9;
	3	2	50	25	0	0	2	1.0	0	230	2	1.1	0.9;
	4	1	30	15	0	0	2	1.0	0	230	2	1.1	0.9;
	5	2	45	22	0	0	3	1.0	0	230	3	1.1	0.9;
	6	1	35	18	0	0	3	1.0	0	230	3	1.1	0.9;
];

%% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin
mpc.gen = [
	1	0	0	200	-200	1.0	100	1	300	0;
	3	0	0	100	-100	1.0	100	1	120	0;
	5	0	0	100	-100	1.0	100	1	120	0;
];

%% fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax
mpc.branch = [
	1	2	0.04	0.50	0.30	300	300	300	0	0	1	-30	30;
	1	3	0.05	0.60	0.30	300	300	300	0	0	1	-30	30;
	3	4	0.04	0.50	0.20	300	300	300	0	0	1	-30	30;
	1	5	0.05	0.60	0.30	300	300	300	0	0	1	-30	30;
	5	6	0.04	0.50	0.20	300	300	300	0	0	1	-30	30;
];

%% model startup shutdown ncost c2 c1 c0
mpc.gencost = [
	2	0	0	3	0	10	0;
	2	0	0	3	0	40	0;
	2	0	0	3	0	40	0;
];
