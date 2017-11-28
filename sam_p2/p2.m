clear all
close all
addpath('../common')


% given params
c = load_project_const();
p2_c = struct();
p2_c.v_fs = -0.59*2;


% ideal quantizer
data_dir = '../sam_cadence/pt_2_tran_psf';

ideal_vod = cds_srr(data_dir, 'tran-tran', 'vod');
ideal_vid = cds_srr(data_dir, 'tran-tran', 'vid');


figure()
plot(ideal_vid.time, ideal_vid.V, 'b-*');
hold on
plot(ideal_vod.time, ideal_vod.V * (p2_c.v_fs / 2) / 0.5, 'r-*');
title('Noiseless Transient Simulation');
legend('vid', 'vod (scaled)');
xlabel('Time');
ylabel('Voltage');
ylim([-1.1, 1.1]);

N = 64;

d = ideal_vod.V(1:N);

fs = 100e6;

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