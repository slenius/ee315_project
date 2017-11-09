clear all
close all

% part C - regeneration time constant
% hw5 p2 pt c

% given params
c = load_project_const();
p1_c = struct();
p1_c.v_fs = 2;

% derived constants
lsb = p1_c.v_fs / c.n_codes;

% update from captab sims
c_l = 17e-15;

% update from dc op point sim
gm_m0 = 151e-6;

tau = c_l / gm_m0;

fprintf('\npart 1.c\n');
fprintf('regeneration time constant = %0.3f ps\n', tau * 1e12);

% calculate the v_id_min from tau and pmeta
% review 6, pg 7
v_id_min = 0.5 * c.p_meta * lsb;

fprintf('\npart 1.d\n');
fprintf('minimum resolvable voltage = %0.3f nV\n', v_id_min * 1e9);

% calculate the number of taus we need to settle
n_taus = log(c.n_codes / c.p_meta);

% dac settling time is equal to the comparator time, and D=0.5, so factor
% of two for computing the clock period
t_clk = 2 * n_taus * tau;
fclk_max = 1 / t_clk;

% number of clocks it takes to convert a sample - from timing diagram pg 3
n_clks_per_conv = 2 + c.bits;

% from clock rate calculate sample rate
t_sample = t_clk * n_clks_per_conv;
f_sample = 1 / t_sample;

fprintf('maximum sample frequency = %0.3f MHz\n', f_sample/1e6);

