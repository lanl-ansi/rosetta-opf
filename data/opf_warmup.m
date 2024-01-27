%
%   This matpower case was developed as part of the rosetta-opf project.
%
%   It is designed to be used for warming up the Julia JIT for solving similar
%   AC-OPF porblems. This network data is designed to cover a wide range of
%   model features in a minimal case file.  It is not representative of a
%   realistic power network.
%
%
function mpc = opf_warmup
mpc.version = '2';
mpc.baseMVA = 100.0;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 1	 150.0	  98.0	 0.1	  10.0	 1	    1.00000	    0.00000	 230.0	 1	    1.10000	    0.90000;
	2	 3	   0.0	   0.0	 0.0	   0.0	 1	    1.00000	    0.00000	 230.0	 1	    1.11000	    0.91000;
	3	 2	 200.0	 -52.0	 0.0	 -15.0	 1	    1.00000	    0.00000	 230.0	 1	    1.12000	    0.92000;
	4	 2	 300.0	 131.0	 0.0	   0.0	 1	    1.00000	    0.00000	 230.0	 1	    1.13000	    0.93000;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
mpc.gen = [
	3	 20.0	 0.0	 30.0	 -30.0	 1.0	 100.0	 1	 40.0	 0.0;
	3	 85.0	 0.0	 127.5	 -127.5	 1.0	 100.0	 1	 170.0	 0.0;
	4	 100.0	 0.0	 150.0	 -150.0	 1.0	 100.0	 1	 200.0	 0.0;
	2	 300.0	 0.0	 450.0	 -450.0	 1.0	 100.0	 1	 600.0	 0.0;
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	   1.00	   5.00	   0.00;
	2	 0.0	 0.0	 3	   0.30	  10.00	   0.00;
	2	 0.0	 0.0	 3	   0.00	  30.00	   0.00;
	2	 0.0	 0.0	 3	   0.00	  40.00	   0.00;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 2	 0.00281	 0.0281	 0.00712	 400.0	 400.0	 400.0	 1.1	  0.0	 1	 -30.0	 30.0;
	2	 3	 0.00108	 0.0108	 0.01852	 426.0	 426.0	 426.0	 0.0	  5.0	 1	 -30.0	 30.0;
	2	 4	 0.00297	 0.0297	 0.00674	 240.0	 240.0	 240.0	 0.0	  0.0	 1	 -30.0	 30.0;
	3	 4	 0.00297	 0.0297	 0.00674	 426.0	 426.0	 426.0	 0.9	 -5.0	 1	 -30.0	 30.0;
];

