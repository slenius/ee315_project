function [s_dbfs, f] = my_psd_dbfs(x, fs, N)
  
  s = abs(fft(x));
  s = s/max(s);
  s = s(1:(1+N/2));
  f = [0:1:(N/2)]*fs/N;
  
  s_dbfs = 20*log10(s);
  
end