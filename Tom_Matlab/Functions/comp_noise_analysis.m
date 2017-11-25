function [ comparator_noise ] = comp_noise_analysis( vdin,p, num_sims )

p=p/num_sims;
initial=[1];
fun=@(sigma,vdinval)((1+erf((vdinval)/(sqrt(2)*sigma)))/(2));
sigma=lsqcurvefit(fun,initial,vdin,p);
figure
stem(vdin*1e3,p,'r','LineWidth',2)
hold on
vdinfit=[-2.5e-3:1e-5:2.5e-3];
plot(vdinfit*1e3,fun(sigma,vdinfit),'k','LineWidth',2)
xlim([-1.2 1.2])
legend('Simulation','Erf Fit','Location','northwest')
xlabel('vid [mV]')
ylabel('Probability')
title('Probability vs Dif Input Voltage')
set(gca,'FontSize',14)
text(-1.1,0.82,sprintf('sigma_n = %4.2f uVRMS',sigma*1e6),'FontSize',14)
comparator_noise=sigma^2;



end

