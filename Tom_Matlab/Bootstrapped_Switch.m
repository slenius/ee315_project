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


%% Sizing the input transistor
% Use Kazim thesis. Extract the input resistance for different widths
load('90nch.mat');
load('90pch.mat');
clc
close all
clearvars W R_on
W=logspace(log10(180e-9),log10(500e-6),100);
L=90e-9;
VGS=1.2;
for i=1:length(W)
    R_on(i)=1./(lookup(nch,'GDS_W','W',W(i),'VGS',VGS,'L',L)*W(i));
    Cgd(i)=lookup(nch,'CGD_W','W',W(i),'VGS',VGS,'L',L)*W(i); 
end

% Calculate total charge injection fro extracted caps and hold cap
Vin=1.2;
C_hold=DAC_unit_cap;
Q=abs((-(Vin+Vdd)*C_hold*Cgd)./(C_hold+Cgd));
V_charge=Q/C_hold;



%% Create figure to find optimum switch width
Font_Size=12;
figure
subplot(1,3,3)
% find second derivative of product
first_deriv=diff(R_on.*Q);
% find maximum in first deriv
[~,index_first_deriv_max]=max(-first_deriv);
opt_width=W(index_first_deriv_max+1);
loglog(W(2:end)*1e6,-first_deriv,'b','LineWidth',2)
hold on
loglog(opt_width*1e6,-first_deriv(index_first_deriv_max),'pk','MarkerSize',12,'MarkerFaceColor','y')
grid on
xlabel('W [\mum]')
title('Negative first derivative of product.')
set(gca,'Fontsize',Font_Size)
subplot(1,3,1)
[AX,R_on_plot,Qinj_plot]=plotyy(W*1e6,R_on*1e-3,W*1e6,Q,'loglog','loglog')
hold on
set(AX,{'ycolor'},{'k';'r'},{'Fontsize'},{Font_Size;Font_Size})
set(R_on_plot,'LineWidth',2,'color','k')
set(Qinj_plot,'LineWidth',2,'color','r')
ylabel(AX(1),'On Resistance [k\Omega]')
ylabel(AX(2),'Charge injected [C]')
xlabel('W [\mum]')
title('R_{on} and Q_{inj} for Varied Transistor Width')

subplot(1,3,2)
loglog(W*1e6,R_on.*Q,'g','LineWidth',2)
hold on
loglog(opt_width*1e6,R_on(index_first_deriv_max+1)*Q(index_first_deriv_max+1),'pk','MarkerSize',12,'MarkerFaceColor','y')

xlabel('W [\mum]')
ylabel('Product R_{on} and Q_{inj}')
grid on
set(gca,'FontSize',Font_Size)
title('Product of R_{on} and Q_{inj}')

%% Size switch and extract capacitance
M10.W=opt_width;
M10.L=L;
M10.CGG=lookup(nch,'CGG_W','W',M10.W,'L',M10.L,'VGS',1.2)*M10.W


%% Size the bootstrapped cap
% this cap needs to be sufficiently large to prevent punchthrough effects
% from significanlty reducing the voltage seen at M10 gate, creating signal
% dependency

% The signal dependence can be seen from lectures 5-6 p 38, with the
% residual signal dependence proportional to Cpar/(Cboot + Cpar). Because
% of this, we want to minimize the parastic capactiance seen by the
% bootstrap cap and increase the bootstrap cap to limit the signal
% dependence.

% As a first pass, size the cap such that it is 10x larger than the M10
% gate capacitance. 

Cboot=5*M10.CGG

%%

Vfs=2;
Vref=Vfs/2;
B=10;
Psig=(0.5*Vfs)^2/2
Delta=Vfs/(2^B)
Pquant=(Delta^2)/12;
SNR_spec=60;
gamma=SNR_spec/10;
k=1.38e-23;
T=300;
% Use equation from Kazim thesis
Chold_min=(k*T)/((Psig)/(10^gamma)-Pquant);
fprintf('Necessary min hold cap size: %4.2f fF\n',Chold_min*1e15)

%%
C=2.5e-12;
SNR=10*log10((alpha)/((beta)/(12)+(k*T/C)))


%% Size the charging branches
% we need to ensure that the the current through the capacitance is
% sufficient to fully charge Cboot back to Vdd in worst case. This speed is
% determined by the speed of the comparator. 

% From the previous design, we determined that the maximum sampling speed
% was 100MHz. This means that the clock is running at ~12GHz. 
fclock=e9;
t_clock=1/fclock;


%% Next, size the PMOS transistor to limit the transit time 

% find the on resistance for the pmos
clearvars R_on Cgd
for i=1:length(W)
    R_on(i)=1./(lookup(pch,'GDS_W','W',W(i),'VGS',VGS,'L',L)*W(i));
    CDD(i)=lookup(pch,'CDD_W','W',W(i),'VGS',VGS,'L',L)*W(i); 
end


figure
subplot(1,2,1)
[AX,R_on_plot,Qinj_plot]=plotyy(W*1e6,R_on*1e-3,W*1e6,CDD,'loglog','loglog')
hold on
set(AX,{'ycolor'},{'k';'r'},{'Fontsize'},{Font_Size;Font_Size})
set(R_on_plot,'LineWidth',2,'color','k')
set(Qinj_plot,'LineWidth',2,'color','r')
ylabel(AX(1),'On Resistance [k\Omega]')
ylabel(AX(2),'Charge injected [C]')
xlabel('W [\mum]')
title('R_{on} and Q_{inj} for Varied Transistor Width')


