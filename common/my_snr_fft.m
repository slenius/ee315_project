function snr_fft = my_snr_fft(x, fs)

N = length(x);

% make sure N is power of two
c = isreal(log2(N)) && rem(log2(N),1)==0;
assert(c, 'N is not a power of two!');
    
% compute SNR from FFT directly
s = abs(fft(x));
s = s/max(s);
s = s(1:N/2);
f = [0:1:(N/2-1)]*fs/N;
[~, i_fund] = max(s);
fin = f(i_fund);
i_3rd = find(f==fin*3);
i_5th = find(f==fin*5);

% compute the power by squaring
r = s.^2;

% remove the DC, fundamental and 3rd/5th harmonics
r_no_h = r;
r_no_h(1) = 0;
r_no_h(i_fund) = 0;
r_no_h(i_3rd) = 0;
r_no_h(i_5th) = 0;

a = r(4);
b = sum(r_no_h);

snr_fft = 10*log10(a/b);

end