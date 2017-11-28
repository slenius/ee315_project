clear all
clc
addpath('Functions');

%% plot GM_ID vs VGS curve
close all
load('90nch.mat');
load('90pch.mat');
Vgs_test=[0.01:0.0075:1.2];
min_W=180e-9;
% extract values from our process
for i=1:length(Vgs_test)
    GM_ID_nch(i)=lookup(nch,'GM_ID','VGS',Vgs_test(i),'W',min_W);
    GM_ID_pch(i)=lookup(pch,'GM_ID','VGS',Vgs_test(i),'W',min_W);
    ID_W_nch(i)=lookup(nch,'ID_W','VGS',Vgs_test(i),'W',min_W);
    ID_W_pch(i)=lookup(pch,'ID_W','VGS',Vgs_test(i),'W',min_W);
end

Font_Size=12;
figure
subplot(1,3,1)
[AX,GM_ID_plot,ID_plot]=plotyy(Vgs_test,GM_ID_nch,Vgs_test,ID_W_nch*1e6*min_W,'plot','semilogy');
set(AX,{'ycolor'},{'k';'r'},{'Fontsize'},{Font_Size;Font_Size})
set(GM_ID_plot,'LineWidth',2,'color','k')
set(ID_plot,'LineWidth',2,'color','r')
ylabel(AX(1),'g_m/I_D')
ylabel(AX(2),'I_D [\mu A]')
xlabel('V_{gs} [mV]')
title('nch')
xlim([min(Vgs_test) max(Vgs_test)])


subplot(1,3,2)
[AX,GM_ID_plot,ID_plot]=plotyy(Vgs_test,GM_ID_pch,Vgs_test,ID_W_pch*min_W*1e6,'plot','semilogy');
set(AX,{'ycolor'},{'k';'r'},{'Fontsize'},{Font_Size;Font_Size})
set(GM_ID_plot,'LineWidth',2,'color','k')
set(ID_plot,'LineWidth',2,'color','r')
ylabel(AX(1),'g_m/I_D')
ylabel(AX(2),'I_D [\mu A]')
xlabel('V_{sg} [mV]')
title('pch')
xlim([min(Vgs_test) max(Vgs_test)])


%W_L=logspace(0.301,3,6);
W_L_lookup_nch=min_W/nch.L(1);
W_L_lookup_pch=min_W/pch.L(1);
ID_W_L_n=(ID_W_nch*min_W)/W_L_lookup_nch;
ID_W_L_p=(ID_W_pch*min_W)/W_L_lookup_pch;

subplot(1,3,3)
for i=1:length(W_L_lookup_nch)
    semilogx(ID_W_L_n,GM_ID_nch,'-','LineWidth',1.5);
    legend_entry{i}=strcat('W_n/L=',num2str(round(W_L_lookup_nch(i))));
    hold on
end
hold on
for i=1:length(W_L_lookup_pch)
    semilogx(ID_W_L_p,GM_ID_pch,'-.','LineWidth',1.5);
    legend_entry{i+length(W_L_lookup_nch)}=strcat('W_p/L=',num2str(round(W_L_lookup_pch(i))));
    hold on
end
legend(legend_entry)
xlabel('ID/(W/L) [\mu A]')
ylabel('g_m/I_D')
title('g_m/I_D vs I_d/(W/L)')
grid minor
set(gca,'Fontsize',Font_Size)
%xlim([min(ID_nch./W_L_lookup(i))*1e6 max(ID_nch./W_L_lookup(1))*1e6])

%hold on
%semilogx(ID_pch./W_L,GM_ID_pch,'-.r','LineWidth',2)

%% Using noise spec to find minimum necessary capacitance at node P,Q

clc
close all
%load('90nch.mat');
%load('90pch.mat');
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

P_meta_spec=10^-7;
Vin_min=2*P_meta_spec*LSB;
fprintf('Req min resolvable input voltage: %4.2f pV\n',Vin_min*1e12)

% use equations from Razavi and lecture for input referred noise voltage
% and gain for strong arm latch
% vin^2=(8*gamma*k*T)/(Av*Cp)
% Av=gm1*Vtn/ICM

%% First find the necessary gain for some design choice gm_ID
clc
close all
% choose GM_ID
clearvars M1
M1.GM_ID=13 ;
for i=1:length(M1.GM_ID)
Vthn=400e-3;
Av=M1.GM_ID(i)*Vthn;
%fprintf('Necessary gain: %4.2f \n',Av)
% extract the capacitance value from vin^2 equation
gamma=0.83706;
k=1.38064852e-23;
T=300;
Cp(i)=(8*gamma*k*T)/(sigma_test^2*Av);
end
fprintf('Minimum capacitance to meet noise spec: %4.2ffF\n',Cp*1e15)




% Use GM_ID to design transistors
% For now, assume everything is minimum length
%%%% DESIGN CHOICE %%%%%%
% budget how much capacitance will be allowed at input
M1_cap_budget=1/10;
M1.L=90e-9;
W=logspace(-7,-4,100);
clearvars CDD
for i=1:length(W)
  M1.CDD(i)=lookup(nch,'CDD_W','GM_ID',M1.GM_ID,'W',W(i),'L',M1.L)*W(i);
end
figure
semilogx(W,M1.CDD*1e15,[min(W) max(W)],[Cp*M1_cap_budget Cp*M1_cap_budget]*1e15)

% find width that meets spec
[~,width_index]=min(abs(M1.CDD-Cp*M1_cap_budget));
M1.W=W(width_index);
M1.CDD=M1.CDD(width_index);
fprintf('Input transistor sizing: %4.2f um\n',M1.W*1e6)
M1.ID=lookup(nch,'ID_W','GM_ID',M1.GM_ID,'L',M1.L)*M1.W;
fprintf('Necessary bias current: %4.2f uA\n',M1.ID*1e6*2)
Vthn=lookup(nch,'VT_ID','GM_ID',M1.GM_ID,'W',M1.W,'L',M1.L,'ID',M1.ID)*M1.ID;
M1.VGS=lookupVGS(nch,'GM_ID',M1.GM_ID,'L',M1.L);
% now that the input transistor is sized, let's size the tail current
% source!
% this transistor just needs to be ensured to drain enough current when
% switched high. This shouldn't be a problem when the device is in triode
% region, as VGS will be very high. To make sure of this function, we size
% the transistor to be large as to reduce the VDS necessary for the
% subsequent biasing of the input common mode
clearvars M7
M7.L=90e-9;
M7.W=M1.W;
M7.ID=M1.ID*2;
%M7.GDS=
M7.GM_ID=lookup(nch,'GM_ID','VGS',Vdd,'W',M7.W,'L',M7.L);


% size the remainder of the circuit














%%

M7.L=0.09e-6;
M7.W=M1.W;
M7.ID=M1.ID*2;
VDS_test=logspace(-3, 1,50);
VGS_test=linspace(0.4,1.2,9);
for i=1:length(VDS_test)
    for j=1:length(VGS_test)
        M7.GM_ID(i,j)=lookup(nch,'GM_ID','VDS',VDS_test(i),'ID',M7.ID,'VGS',VGS_test(j),'L',M7.L);
    end
end
figure
for j=1:length(VGS_test)
plot(VDS_test,M7.GM_ID(:,j),'LineWidth',3);
legend_entry{j}=strcat('VGS=',num2str(VGS_test(j)));
hold on
end
xlabel('VDS [V]')
ylabel('GM_ID')
legend(legend_entry)
set(gca,'fontsize',14)
%close all
M7.VGS=0.7;
[~,VGS_index]=min(abs(M7.VGS-VGS_test));
M7.VDS=0.1;
M1.VCM=M1.VGS+M7.VDS
[~,VDS_index]=min(abs(M7.VDS-VDS_test));
M7.GM_ID=M7.GM_ID(VDS_index,VGS_index);
M7.W=lookup(nch,'W_ID','GM_ID',M7.GM_ID,'VDS',M7.VDS,'L',M7.L)*M7.ID


%%

M7.Vds=1/(lookup(nch,'GDS_ID','GM_ID',M7.GM_ID,'ID',M7.ID,'L',M7.L));
%
M7.W=lookup(nch,'W_ID','GM_ID',M7.GM_ID,'ID',M7.ID,'L',M7.L)*M7.ID;
fprintf('Width of tail transistor: %4.2fum \n',M7.W*1e6)
M7.gds=lookup(nch,'GDS_ID','GM_ID',M7.GM_ID,'ID',M7.ID,'W',M7.W,'L',M7.L)*M7.ID;
M7.Vds=M7.ID/M7.gds
M7.Vgs=lookupVGS(nch,'GM_ID',M7.GM_ID,'VDS',M7.Vds,'VSB',0,'L',M7.L);
fprintf('Tail transistor biasing: %4.2f V\n',M7.Vgs)

%% Calculate the timing of comparator
load('/afs/ir.stanford.edu/users/t/o/tomflo/ee315/opus/Simulation_Results/comp_initial_sim_results.mat');
N_bit=10;
clc
num_clock_sims=1;
[ comp_dec_time, reg_time, recharge_time ] = Comparator_Max_Freq( comp_initial_sim_results.StrongArmPos, comp_initial_sim_results.StrongArmNeg, comp_initial_sim_results.vclock, num_clock_sims );
fprintf('***** Simulations for minimum input voltage %4.2f pV ******\n\n',Vin_min*1e12)
fprintf('Comparator decision time: %4.2f ps\n',comp_dec_time*1e12)
fprintf('Regeneration time: %4.2f ps\n',reg_time*1e12)
fprintf('Recharge time: %4.2f ps\n',recharge_time*1e12)
if reg_time>recharge_time
    fprintf('Regeneration time is %4.2f ps slower than recharge time\n',(reg_time-recharge_time)*1e12)
else 
    fprintf('Regeneration time is %4.2f ps faster than recharge time\n',(-reg_time+recharge_time)*1e12)

end
fprintf('Max clock frequency: %4.2f GHz\n',1/((max([comp_dec_time recharge_time]))*2)*1e-9)
fprintf('Max synchronous sample frequency: %4.2f MHz\n',1/((max([comp_dec_time recharge_time]))*2*(N_bit+2))*1e-6)
fprintf('Approximate asynchronous max frequency: > %4.2f MHz\n',(1/((max([comp_dec_time recharge_time]))*2*(N_bit+2))*1e-6)*2)

%% Find input referred noise offset
vdin=[-1e-3 -600e-6 -200e-6 200e-6 600e-6 1e-3];
p=[5 19 83 146 184 198 ];
num_sims=200;
[ comparator_noise ] = comp_noise_analysis( vdin,p, num_sims );


