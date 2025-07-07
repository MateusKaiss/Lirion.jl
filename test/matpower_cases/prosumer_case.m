function mpc = case3_timeseries_example
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 0.0	 0.0	 0.0	 0.0	 1	1.00000	0.00000	240.0	1	1.10000	0.90000;
	2	 2	 0.0	 0.0	 0.0	 0.0	 1	1.00000	0.00000	240.0	1	1.10000	0.90000;
	3	 2	 0.0	 0.0	 0.0	 0.0	 1	1.00000	0.00000	240.0	1	1.10000	0.90000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
mpc.gen = [
	1	1000.0	0.0	1000.0	-1000.0	1.0	100.0	1	2000.0	0.0;
	2	1000.0	0.0	1000.0	-1000.0	1.0	100.0	1	2000.0	0.0;
	3	   0.0	0.0	1000.0	-1000.0	1.0	100.0	1	   0.0	0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0.0	0.0	3	0.110000	5.000000	0.000000;
	2	0.0	0.0	3	0.085000	1.200000	0.000000;
	2	0.0	0.0	3	0.000000	0.000000	0.000000;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	3	0.065	0.62	0.45	9000.0	9000.0	9000.0	0.0	0.0	1	-30.0	30.0;
	3	2	0.025	0.75	0.7	50.0	50.0	50.0	0.0	0.0	1	-30.0	30.0;
	1	2	0.042	0.9	0.3	9000.0	9000.0	9000.0	0.0	0.0	1	-30.0	30.0;
];

%% bus timeseries data
%	time	bus	pd	qd
mpc.bus_timeseries = [
	1	1	110.0	40.0;
	1	2	90.0	30.0;
	1	3	95.0	50.0;
	2	1	120.0	41.0;
	2	2	95.0	31.0;
	2	3	100.0	51.0;
];

%% consumer timeseries data
%	time	bus	demand	q_demand
mpc.consumer_timeseries = [
	1	2	50.0	10.0;
	1	3	40.0	5.0;
	2	2	55.0	12.0;
	2	3	42.0	6.0;
];

%% prosumer static data
%	bus	generation	demand	q_generation	q_demand	cost_injection	inverter_capacity
mpc.prosumer = [
	2	10.0	5.0	2.0	1.0	50.0	100.0;
	3	20.0	10.0	4.0	2.0	55.0	100.0;
];

%% prosumer timeseries data
%	time	bus	generation	demand	q_generation	q_demand
mpc.prosumer_timeseries = [
	1	2	15.0	6.0	2.5	1.2;
	1	3	22.0	11.0	4.2	2.1;
	2	2	16.0	6.5	2.8	1.3;
	2	3	24.0	12.0	4.5	2.2;
];

%% generator cost timeseries data
%   time   gen_id   model   startup   shutdown   ncost   [c(n-1) ... c0]   tou   feed_in_tariff
mpc.gencost_timeseries = [
    1   1   2   0.0   0.0   3   0.110000   5.000000   0.000000   1.0   0.0;
    1   2   2   0.0   0.0   3   0.085000   1.200000   0.000000   1.0   0.0;
    2   1   2   0.0   0.0   3   0.120000   5.200000   0.000000   1.2   0.0;
    2   2   2   0.0   0.0   3   0.088000   1.250000   0.000000   1.2   0.0;
];