% EE315 Project
% First Deliverable
cd First_Deliverable/Inverter_Design
%% Part A - Inverter Design

clc
clear all
close all

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
fprintf('Inverter Input Cap: %4.2ffF\n',C_in*1e15w)
N=(log(C_load/C_in))/log(4);
fprintf('Number of inverters in FO4 chain: %i\n',round(N))

% figure
% plot(Sim.out_2f(:,1)*1e12,Sim.out_2f(:,2))
%% Thermal Noise Voltage - Hand Analysis
fprintf('\n\n**** Part E ****\n')
fprintf('Thermal Noise Hand Analysis\n\n')
% constants
k=1.38064852e-23;
T=25+273; %K

% collect values from simulation
gm_n=143.3e-6;
gm_p=151.2e-6;
gds_n=12.15e-6;
gds_p=4.852e-6;
C_L=1.47e-14;



% Calculate Gain
A_v=-(gm_n)/(gds_n+gds_p-gm_p);

% Calculate total output noise
gamma=0.25; %<----- Needs to be derived from simulation!!!
Req=1/(gds_n+gds_p);% r_on||r_op
p_noise=4*k*T*gamma/gm_p;
n_noise=gamma*gm_n*Req*k*T/C_L;
V_noise_tot_out=(p_noise+n_noise)*2
% Calculate total input noise
V_noise_tot_in=V_noise_tot_out/(abs(A_v)^2);

fprintf('Total output referred noise voltage: %4.2fmV RMS\n',sqrt(V_noise_tot_out*1e6))
fprintf('Total input referred  noise voltage: %4.2fmV RMS\n',sqrt(V_noise_tot_in*1e6))

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

vdin=[-1e-3 -200e-6 200e-6 1e-3];
p=[0.2 0.4 0.6 0.8];
initial=[1 1 1];
x=lsqnonlin(@erf_fit,initial,[],[],[],vdin,p)

