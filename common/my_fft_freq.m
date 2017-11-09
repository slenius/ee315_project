function f = my_fft_freq(fs, N)

f = (1:N/2) ./ (N/2);
f = f * fs/2;

end