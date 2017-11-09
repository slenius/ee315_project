clear all
close all
addpath('../common')
% part C - regeneration time constant
% hw5 p2 pt c

% given params
c = load_project_const();
p1_c = struct();
p1_c.v_fs = 2;

% derived constants
lsb = p1_c.v_fs / c.n_codes;

% update from captab sims
%c_l = 17e-15;

% update from dc op point sim
gm_m0 = 151e-6;
c_l = 14.4e-15;

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


% make fft of ideal and nonideal quantizer
spectre = '/afs/ir.stanford.edu/class/ee/cadence/MMSIM10/tools/spectre/matlab';
old_path = path;
path(old_path, spectre);

% ideal quantizer
data_dir = '../sam_cadence/pt_1_tran_psf';

ideal_vod = cds_srr(data_dir, 'tran-tran', 'vod');
ideal_vid = cds_srr(data_dir, 'tran-tran', 'vid');


figure()
plot(ideal_vid.time, ideal_vid.V, 'b-*');
hold on
plot(ideal_vod.time, ideal_vod.V*2, 'r-*');
title('Noiseless Transient Simulation');
legend('vid', 'vod*2');
xlabel('Time');
ylabel('Voltage');
ylim([-1.1, 1.1]);

N = 64;

d = ideal_vod.V(1:N);

fs = 21e6;

[s_dbfs, f] = my_psd_dbfs(d, fs, N);
[~, idx] = max(s_dbfs);
input_freq = f(idx);
fprintf('Input Frequency = %f MHz\n', input_freq/1e6);

% compute snr
snr_fft = my_snr_fft(d, fs);

figure();
plot(f/1e6, s_dbfs)
title_s = sprintf('Noiselsess Transient Simulation - SNR = %0.2fdB', snr_fft);
title(title_s);
ylabel('dBFS');
xlabel('Frequency (MHz)')

fprintf('Noiseless SNR: %0.3fdB\n', snr(d, fs, 6));
fprintf('Noiseless SNDR: %0.3fdB\n', sinad(d, fs));


% ideal quantizer
data_dir = '../sam_cadence/pt_1_tran_noise_psf';
%data_dir = '../sam_cadence/simulation/project_starter/spectre/schematic/psf';

noisy_vod = cds_srr(data_dir, 'tran-tran', 'vod');
noisy_vid = cds_srr(data_dir, 'tran-tran', 'vid');

figure()
plot(noisy_vid.time, noisy_vid.V, 'b-*');
hold on
plot(noisy_vod.time, noisy_vod.V*2, 'r-*');
title('Noisy Transient Simulation');
legend('vid', 'vod*2');
xlabel('Time');
ylabel('Voltage');
ylim([-1.1, 1.1]);

d = noisy_vod.V(1:N);

% compute snr
snr_fft = my_snr_fft(d, fs);

% compute fft
[s_dbfs, f] = my_psd_dbfs(d, fs, N);
[val, idx] = max(s_dbfs);
input_freq = f(idx);
fprintf('Input Frequency = %f MHz\n', input_freq/1e6);

figure();
plot(f/1e6, s_dbfs)
title_s = sprintf('Noisy Transient Simulation - SNR = %0.2fdB', snr_fft);
title(title_s);
ylabel('dBFS');
xlabel('Frequency (MHz)')




