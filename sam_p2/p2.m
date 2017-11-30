clear all
close all
addpath('../common')


% given params
c = load_project_const();
p2_c = struct();
p2_c.v_fs = -0.59*2;

% worst case async speedup
N = c.bits;
speedup_t = ((N-1) * log(3) + log(2) + N/2*(N+1)*log(2)) / (N*(N+1)*log(2));
speedup = 1/speedup_t;

% calculate sync speed
tau = 16.7e-12;
n_taus = 20.7;
regen_time = n_taus * tau;
fo4_time = 40e-12;
slew_time = 64e-12;
high_time = regen_time + 2 * fo4_time + slew_time;
bit_time = high_time * 2;
sync_time = bit_time * 12;
sync_speed = 1/sync_time;
async_speed = speedup * sync_speed;


% ideal quantizer

figure();

%data_dir = '../sam_cadence/pt_2_tran_psf';
data_dir = '../sam_cadence/pt_2_tran_noise_psf_moderate_32_64_130';


nyquist_vod = cds_srr(data_dir, 'tran-tran', 'vod');
nyquist_vid = cds_srr(data_dir, 'tran-tran', 'vid');

N = 64;

d = nyquist_vod.V(1:N);

fs = 130e6;

[s_dbfs, f] = my_psd_dbfs(d, fs, N);
[~, idx] = max(s_dbfs);
input_freq = f(idx);
fprintf('Sampling Frequency = %f MHz\n', fs/1e6);
fprintf('Input Frequency = %f MHz\n', input_freq/1e6);

% compute snr
snr_fft = my_snr_fft(d, fs);

%figure()
subplot(4,1,1);
plot(nyquist_vid.time, nyquist_vid.V, 'b-*');
hold on
plot(nyquist_vod.time, nyquist_vod.V * (p2_c.v_fs / 2) / 0.5, 'r-*');
title_s = sprintf('Transient Noise Sim - Fin=%0.1fMHz', input_freq/1e6);
title(title_s);
legend('vid', 'vod (scaled)');
xlabel('Time');
ylabel('Voltage');
ylim([-1.1, 1.1]);



subplot(4,1,2);
plot(f/1e6, s_dbfs)
title_s = sprintf('Transient Noise Sim - Fin=%0.1fMHz - SNR=%0.2fdB', input_freq/1e6, snr_fft);
title(title_s);
ylabel('dBFS');
xlabel('Frequency (MHz)')

fprintf('SNR: %0.3fdB\n', snr(d, fs, 6));
fprintf('SNDR: %0.3fdB\n', sinad(d, fs));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_dir = '../sam_cadence/pt_2_tran_noise_psf_moderate_3_64_130';

lowf_vod = cds_srr(data_dir, 'tran-tran', 'vod');
lowf_vid = cds_srr(data_dir, 'tran-tran', 'vid');

d = lowf_vod.V(1:N);


[s_dbfs, f] = my_psd_dbfs(d, fs, N);
[~, idx] = max(s_dbfs);
input_freq = f(idx);
fprintf('Sampling Frequency = %f MHz\n', fs/1e6);
fprintf('Input Frequency = %f MHz\n', input_freq/1e6);

% compute snr
snr_fft = my_snr_fft(d, fs);

subplot(4,1,3);
plot(lowf_vid.time, lowf_vid.V, 'b-*');
hold on
plot(lowf_vod.time, lowf_vod.V * (p2_c.v_fs / 2) / 0.5, 'r-*');
title_s = sprintf('Transient Noise Sim - Fin=%0.1fMHz', input_freq/1e6);
title(title_s);
legend('vid', 'vod (scaled)');
xlabel('Time');
ylabel('Voltage');
ylim([-1.1, 1.1]);



subplot(4,1,4);
plot(f/1e6, s_dbfs)
title_s = sprintf('Transient Noise Sim - Fin=%0.1fMHz - SNR=%0.2fdB', input_freq/1e6, snr_fft);
title(title_s);
ylabel('dBFS');
xlabel('Frequency (MHz)')

fprintf('SNR: %0.3fdB\n', snr(d, fs, 6));
fprintf('SNDR: %0.3fdB\n', sinad(d, fs));
