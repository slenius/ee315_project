function [ sigma_n ] = erf_noise( vdin, p )

sigma_n=1/(erfinv(2*p-1)*(sqrt(2))/(vdin));
end

