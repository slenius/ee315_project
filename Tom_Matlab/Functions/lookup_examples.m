% Usage examples for function "lookup"
% Boris Murmann
% Stanford University
%clearvars;
close all;

% Load data tables (this loads the structures nch and pch)
%load 90nch.mat;
%load 90pch.mat;
dev = nch;

% Plot drain characteristics for different VGS at minimum L (default value)
vgs = 0:0.1:max(nch.VGS);
id = lookup(dev, 'ID', 'VGS', vgs, 'VDS', dev.VDS);
figure;
plot(nch.VDS, id')
xlabel('V_D_S [V]')
ylabel('I_D [A]')
%%
% Plot gm_ID vs. VGS for different L
gm_ID1 = lookup(dev, 'GM_ID', 'VGS', dev.VGS, 'L', dev.L);
figure;
plot(dev.VGS, gm_ID1)
ylabel('g_m/I_D [S/A]')
xlabel('V_G_S [V]')

%% gm/ID vector for remaining plots
gm_ID = 3:0.2:min(max(gm_ID1'));

% Plot fT against gm_ID for different L and VDS
wt = lookup(dev, 'GM_CGG', 'GM_ID', gm_ID, 'L', dev.L);
figure;
semilogy(gm_ID, wt/2/pi)
xlabel('g_m/I_D [S/A]')
ylabel('f_T [Hz]')
legend('min')

%%
% Plot ID/W against gm_ID for different L
% Note that VDS is not specified here; it then defaults to max(nch.VDS)/2
id_w = lookup(dev, 'ID_W', 'GM_ID', gm_ID, 'L', dev.L');
figure;
semilogy(gm_ID, id_w)
xlabel('g_m/I_D [S/A]')
ylabel('I_D/W [A/m]')

% Plot gm/gds against gm_id at two different L and default VDS
lmin = min(dev.L);
gm_gds = lookup(dev, 'GM_GDS', 'GM_ID', gm_ID, 'L', [lmin 2*lmin]);
figure;
plot(gm_ID, gm_gds)
xlabel('g_m/I_D [S/A]')
ylabel('g_m/g_d_s')

% Plot gamma factor versus gm/ID
kBoltz = 1.38064852e-23;
gamma = lookup(dev, 'STH_GM', 'GM_ID', gm_ID)/4/kBoltz/dev.TEMP;
figure;
plot(gm_ID, gamma, [0 max(gm_ID)], 2/3*[1 1], 'k--');
xlabel('g_m/I_D [S/A]')
ylabel('\gamma')

% Find thermal noise factor in triode region (should be 1)
gamma_triode = lookup(dev, 'STH_GDS', 'VGS', max(dev.VGS), 'VDS', 0)/4/kBoltz/dev.TEMP

% Given gm_ID, VDS, VSB and L, find VGS
vgs1 = lookupVGS(dev, 'GM_ID', [10 12], 'VDS', 0.6, 'VSB', 0.3, 'L', 0.2)

% Given gm_ID and everything else at default, find VGS
vgs2 = lookupVGS(dev, 'GM_ID', 10)


%%
clc
CDD_W=lookup(nch,'CDD_W','GM_ID',10);
W=50e-6;
CGG=lookup(nch,'CGG_W','GM_ID',10)*W
CDD=CDD_W*W

