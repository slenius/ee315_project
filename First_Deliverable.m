% EE315 Project
% First Deliverable
clc
clear all
close all
cd First_Deliverable/
load('Gamma.mat');
cd Inverter_Design/
%% Part A - Inverter Design



fprintf('**** Part A ****\n')
fprintf('Inverter Design\n\n')

% Design an inverter chain to buffer the comparator outputs and drive a 
% load capacitance of 50 fF. Use a fan-out (FO) of 4 for each stage. 
% Instantiate this inverter chain in ?hw_comparator_testbench? schematic
% (in Project_Lib2 library). Use this inverter chain together with the 
% comparator in all following parts.

% Optimal number of inverters will be written as: 
% N=log_4(C_load/C_in)
%
C_load=50e-15;
% Need to find C_in for minimum sized transistor. Get this value from
% simulation
load('Inverter_Design/CompleteSimulation.mat')
load('Inverter_Design/SimulationSingle.mat')
Vdd=1.2;
C_test=[3e-16:0.2e-16:6e-16];
for i=1:length(C_test)
[ rise_time(i), fall_time(i) ] = rise_fall_time( CompleteSimulation(:,[2*i-1:2*i]), Vdd );
end
mean_time=(rise_time+fall_time)./2;
[rise_time_inv, fall_time_inv]=rise_fall_time( SimulationSingle, Vdd );
mean_time_inv=(rise_time_inv+fall_time_inv)/2;
[~,cap_index]=min(abs(mean_time-mean_time_inv));
C_in=C_test(cap_index);
fprintf('Inverter Input Cap: %4.2ffF\n',C_in*1e15)
N=(log(C_load/C_in))/log(4);
fprintf('Number of inverters in FO4 chain: %i\n',round(N))

% figure
% plot(Sim.out_2f(:,1)*1e12,Sim.out_2f(:,2))
%% Thermal Noise Voltage - Hand Analysis
fprintf('\n\n**** Part E ****\n')
fprintf('Thermal Noise Hand Analysis\n\n')
% constants
k=1.38064852e-23;
T=300; %K

% collect values from simulation
gm_n=143.3e-6;
gm_p=151.2e-6;
gds_n=12.15e-6;
gds_p=4.852e-6;
C_L=1.47e-14;


% Derive gamma from simulation, using Vdc biased -> 0.6V and Vds 1.2V
% Gamma will be beyond the flicker corner

gamma=Gamma(end,2);
figure
loglog(Gamma(:,1),Gamma(:,2),'r','LineWidth',2)
hold on
text(1e10,1.3,strcat({'\downarrow \gamma \approx '},{num2str(gamma)}),'fontsize',14)
grid on
xlabel('Frequecy [Hz]')
ylabel('\gamma')
title('Noise Simulation to Extract Gamma')
set(gca,'fontsize',14)
gamma=Gamma(end,2);



% Calculate Gain
A_v=-(gm_n)/(gds_n+gds_p-gm_p);

% Calculate total output noise
Req=1/(gds_n+gds_p);% r_on||r_op
V_noise_tot_out=2*((k*T*gamma*(gm_n+gm_p)*Req)/C_L);
% Calculate total input noise
V_noise_tot_in=V_noise_tot_out/(abs(A_v)^2);

fprintf('Total output referred noise voltage: %4.2f mVRMS\n',sqrt(V_noise_tot_out*1e6))
fprintf('Total input referred  noise voltage: %4.2f mVRMS\n',sqrt(V_noise_tot_in*1e6))

%% Thermal Noise Voltage - Simulation
fprintf('\n\n**** Part F ****\n')
fprintf('Thermal Noise Simulation\n\n')
V_noise_tot_out_sim=0.001026832; %V rms
fprintf('Total output referred noise voltage from simulation: %4.2fmV RMS\n',V_noise_tot_out_sim*1e3)
fprintf('Percent error wrt hand analysis: %4.2f%%\n',(V_noise_tot_out_sim-sqrt(V_noise_tot_out))/(sqrt(V_noise_tot_out))*100);
V_noise_tot_in_sim=(V_noise_tot_out_sim/A_v);
fprintf('Total input referred noise voltage from simulation: %4.2fmV RMS\n',V_noise_tot_in_sim*1e3)
fprintf('Percent error wrt hand analysis: %4.2f%%\n',(V_noise_tot_in_sim-sqrt(V_noise_tot_in))/(sqrt(V_noise_tot_in))*100);



%% Part G - Transient Noise Simulations
tau=95.3e-12;
FMAX=10/(2*pi*tau);
vdin=[-1e-3 -200e-6 200e-6 1e-3];
p=[82 420 630 920]/1000;
initial=[1];
fun=@(sigma,vdinval)((1+erf((vdinval)/(sqrt(2)*sigma)))/(2));
options = optimset('Display','off');
sigma=lsqcurvefit(fun,initial,vdin,p,[],[],options);
figure
stem(vdin*1e3,p,'r','LineWidth',2)
hold on
vdinfit=[-2.5e-3:1e-5:2.5e-3];
plot(vdinfit*1e3,fun(sigma,vdinfit),'k','LineWidth',2)
xlim([-2.5 2.5])
legend('Simulation','Erf Fit')
xlabel('vid [mV]')
ylabel('Probability')
title('Probability vs Dif Input Voltage')
set(gca,'FontSize',14)
text(-2,0.9,strcat({'\sigma_n = '},{num2str(sigma*1e3)},{' mV_{RMS}'}),'FontSize',14)
comparator_noise=sigma^2;
fprintf('Noise compared to Hand Analysis and Simulation: %4.2f%% and %4.2f%%\n',(sigma-(sqrt(V_noise_tot_in)))/(sqrt(V_noise_tot_in))*100,(sigma-(V_noise_tot_in_sim))/(V_noise_tot_in_sim)*100)
%% Part H - DAC Cap and Sampling Noise
fprintf('\n\n**** Part H ****\n')
unit_cap=2.5e-15;
multipliers=[1 1 2 4 8 16 32 64 128];
CDAC_Single=sum(multipliers*unit_cap);
fprintf('C_DAC for a single DAC: %4.2fpF or %iCu\n',CDAC_Single*1e12,CDAC_Single/unit_cap)
fprintf('C_DAC for total DAC: %4.2fpF or %iCu\n',CDAC_Single*1e12*2,CDAC_Single/unit_cap*2)
total_sampling_noise=(2*k*T/CDAC_Single);
fprintf('Total sampling noise: %4.2f mVRMS\n',sqrt(total_sampling_noise)*1e3)

%% Part I - Hand Calculated SNR
fprintf('\n\n**** Part I ****\n')
V_FS=[-1 1];
Vpeak=max(V_FS);
Vrms=Vpeak/sqrt(2);


V_ref=1;
delta=2*V_ref/(2^10);
quant_noise=delta^2/12;

fprintf('Quantization noise: %4.2f mVRMS\n',sqrt(quant_noise)*1e3)
total_noise=quant_noise+total_sampling_noise+comparator_noise;
fprintf('Total noise: %4.2f mVRMS\n',sqrt(total_noise)*1e3)
V_signal=Vrms;
SNR=10*log10(V_signal^2/total_noise);
SQNR=10*log10(V_signal^2/quant_noise);
fprintf('SNR and SQNR: %4.2fdB and %4.2fdB\n',SNR,SQNR) 
SQNR_ideal=6.02*10+1.76;

