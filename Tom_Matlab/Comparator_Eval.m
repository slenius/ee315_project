clear all
clc
addpath('Functions');
%% Calculate the regeneration time constant
%load('/afs/ir.stanford.edu/users/t/o/tomflo/ee315/opus/Simulation_Results/comp_initial_sim_results');

%% Sizing the comparator
clc
clear all
close all
load('90nch.mat');
% Designing the comparator starting from the noise requirement
clearvars W_input
SNR_spec_target=55;
SNR_spec=SNR_spec_target+5;
Vdd=1.2;

% design such that sigma noise = LSB
B=10;
V_FS=2;
LSB=V_FS/(2^B);
Psig=(0.5*V_FS)^2/2;
P_noise=Psig/(10^(SNR_spec/10));
sigma_test=sqrt(P_noise);
fprintf('Necessary noise spec: %4.2f uVRMS\n',sigma_test*1e6);

% use equations from Razavi and lecture for input referred noise voltage
% and gain for strong arm latch
% vin^2=(8*gamma*k*T)/(Av*Cp)
% Av=gm1*Vtn/ICM

% First find the necessary gain for some design choice gm_ID
% choose GM_ID
GM_ID=[12:0.1:19];
for i=1:length(GM_ID); %reasonable speed/power tradeoff
    Vthn=250e-3;
    Av=GM_ID(i)*Vthn;
    %fprintf('Necessary gain: %4.2f \n',Av)
    
    % extract the capacitance value from vin^2 equation
    gamma=0.83706;
    k=1.38064852e-23;
    T=300;
    Cp(i)=(8*gamma*k*T)/(sigma_test^2*Av);
    %fprintf('Minimum capacitance to meet noise spec: %4.2ffF\n',Cp*1e15)
    
    % Use GM_ID to pull out some transistor sizes
    GM_CDD= lookup(nch, 'GM_CDD', 'GM_ID', GM_ID(i),'VDS',Vdd);
    gm(i)=GM_CDD*Cp(i);
    ID(i)=gm(i)/GM_ID(i);
    W_ID=lookup(nch,'W_ID','GM_ID',GM_ID(i));
    W_input(i)=ID(i)*W_ID;
end
figure
subplot(3,1,1)
plot(GM_ID,W_input*1e6)
title('GM_ID vs W')
subplot(3,1,2)
plot(GM_ID,gm,'r')
title('GM_ID vs gm')

subplot(3,1,3)
plot(GM_ID,ID*1e6,'b')
title('GM_ID vs ID')

%%

clc
clear all
close all
load('90nch.mat');
% Designing the comparator starting from the noise requirement
clearvars W_input
SNR_spec_target=55;
SNR_spec=SNR_spec_target+5;

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

% use equations from Razavi and lecture for input referred noise voltage
% and gain for strong arm latch
% vin^2=(8*gamma*k*T)/(Av*Cp)
% Av=gm1*Vtn/ICM

% First find the necessary gain for some design choice gm_ID
% choose GM_ID
GM_ID=15;
Vthn=250e-3;
Av=GM_ID*Vthn;
%fprintf('Necessary gain: %4.2f \n',Av)

% extract the capacitance value from vin^2 equation
gamma=0.83706;
k=1.38064852e-23;
T=300;
Cp=(8*gamma*k*T)/(sigma_test^2*Av);
fprintf('Minimum capacitance to meet noise spec: %4.2ffF\n',Cp*1e15)

% pull Vgs value from GM_ID=2/Vov
Vs1=400e-3; % approximations
Vov=2/GM_ID;
Vg1=Vov+Vs1+Vthn;
Vgs1=Vg1-Vs1;
fprintf('Necessary common mode biasing: %4.2f mV\n',Vg1*1e3)
%% plot GM_ID vs VGS curve
close all
clearvars GM_ID_nch GM_ID_pch ID_nch ID_pch
Vgs_test=[0.01:0.0075:1.2];
% extract values from our process
for i=1:length(Vgs_test)
   GM_ID_nch(i)=lookup(nch,'GM_ID','VGS',Vgs_test(i));
   GM_ID_pch(i)=lookup(pch,'GM_ID','VGS',Vgs_test(i));
   ID_nch(i)=lookup(nch,'ID','VGS',Vgs_test(i),'VDS',0.4);
   ID_pch(i)=lookup(pch,'ID','VGS',Vgs_test(i),'VDS',0.4);
end

Font_Size=12;
figure
subplot(1,3,1)
[AX,GM_ID_plot,ID_plot]=plotyy(Vgs_test,GM_ID_nch,Vgs_test,ID_nch*1e6,'plot','semilogy');
set(AX,{'ycolor'},{'k';'r'},{'Fontsize'},{Font_Size;Font_Size}) 
set(GM_ID_plot,'LineWidth',2,'color','k')
set(ID_plot,'LineWidth',2,'color','r')
ylabel(AX(1),'g_m/I_D')
ylabel(AX(2),'I_D [\mu A]')
xlabel('V_{gs} [mV]')
title('nch')
xlim([min(Vgs_test) max(Vgs_test)])


subplot(1,3,2)
[AX,GM_ID_plot,ID_plot]=plotyy(Vgs_test,GM_ID_pch,Vgs_test,ID_pch*1e6,'plot','semilogy');
set(AX,{'ycolor'},{'k';'r'},{'Fontsize'},{Font_Size;Font_Size}) 
set(GM_ID_plot,'LineWidth',2,'color','k')
set(ID_plot,'LineWidth',2,'color','r')
ylabel(AX(1),'g_m/I_D')
ylabel(AX(2),'I_D [\mu A]')
xlabel('V_{sg} [mV]')
title('pch')
xlim([min(Vgs_test) max(Vgs_test)])


W_L=[2:4:10];
subplot(1,3,3)
for i=1:length(W_L)
semilogx(ID_nch./W_L(i)*1e6,GM_ID_nch,'-','LineWidth',1.5);
legend_entry{i}=strcat('W_n/L=',num2str(W_L(i)));
hold on
end
hold on
for i=1:length(W_L)
semilogx(ID_pch./W_L(i)*1e6,GM_ID_pch,'-.','LineWidth',1.5);
legend_entry{i+length(W_L)}=strcat('W_p/L=',num2str(W_L(i)));
hold on
end
legend(legend_entry)
xlabel('ID/(W/L) [\mu A]')
ylabel('g_m/I_D')
title('g_m/I_D vs I_d/(W/L)')
grid minor
set(gca,'Fontsize',Font_Size)
xlim([min(ID_nch./W_L(i)*1e6) max(ID_nch./W_L(i)*1e6)])

%hold on
%semilogx(ID_pch./W_L,GM_ID_pch,'-.r','LineWidth',2)


%% Use curves to design transistors
clc
% For now, assume everything is minimum length
% select a reasonable GM_ID from curves for the input transistor pair
GM_ID_input=10;
[~,nch_index]=min(abs(GM_ID_nch-GM_ID_input))
Vgs1=Vgs_test(nch_index);
W_L=2; % keep minimum!
% find index for GM_ID vs ID curve

I_D=ID_nch(nch_index)
I_bias=2*I_D
% find GM_ID for the tail
GM_ID_tail=lookup(nch,'GM_ID','ID',I_bias,'VGS',0.6)
% find the Vds of the tail transistor
Vds_tail=lookup(nch,'ID_GDS','GM_ID',GM_ID_tail)
% Find common mode biasing for Min
V_CM_in=Vgs1+Vds_tail
