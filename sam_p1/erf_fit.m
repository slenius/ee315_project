function fit = erf_fit( x,X,Y )
% This function is called by lsqnonlin.
% x is a vector which contains the coefficients of the
% equation.  X and Y are the option data sets that were
% passed to lsqnonlin.

p=x(1);
vdin=x(2);
sigma_n=x(3);


fit=(1+erf(vdin/(sqrt(2)*sigma_n)))/(2)-p;

end

