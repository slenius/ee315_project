clear all
clc
close all
load('/afs/ir.stanford.edu/users/t/o/tomflo/ee315/opus/Simulation_Results/comparator_DC.mat');
load('90nch.mat');
load('90pch.mat');
% Noise requirement
clearvars W_input
SNR_spec_target=55;
SNR_spec=SNR_spec_target+5;
k=1.38e-23;
T=300;

Vdd=1.2;

% design such that sigma noise = LSB
B=10;
V_FS=2;
LSB=V_FS/(2^B);
Psig=(0.5*V_FS)^2/2;
P_noise=Psig/(10^(SNR_spec/10));
sigma_test=sqrt(P_noise);
fprintf('Necessary noise spec: %4.2f uVRMS\n',sigma_test*1e6);
fprintf('LSB: %4.2f uV\n',LSB*1e6)

P_meta_spec=10^-7;
Vin_min=2*P_meta_spec*LSB;
fprintf('Req min resolvable input voltage: %4.2f pV\n',Vin_min*1e12)

% Input noise from comparator, designed in previous script
comp_noise=500e-6; %Vrms
% extract parastic capacitance from comparator input pair
comp_input_cap=Minput.cgg;
% unit capacitance for DAC
DAC_unit_cap=2.5e-15;
C_par=DAC_unit_cap+comp_input_cap;
kT_C_noise=k*T/C_par;
kT_C_noise_rms=sqrt(kT_C_noise);
fprintf('Input kT/C noise: %4.2f uVrms\n',sqrt(kT_C_noise)*1e6)
if sqrt(kT_C_noise)<sigma_test
fprintf('Noise spec met with %4.2f uVrms to spare\n',(sigma_test-kT_C_noise_rms)*1e6)
else
    fprintf('Noise spec not met! Over by %4.2f uVrms\n',(kT_C_noise_rms-sigma_test)*1e6)

end


%%