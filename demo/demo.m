% Program that calculates the buoyancy work required to bring a water parcel
% from different depths (zint) to a reference depth (default = -10 m).
% The ocean is considered to be at rest and the parcel at zint is in hydrostatic equilibrium.
% The vertical discretization of the vertical profile may be non-equidistant.

clearvars; close all; clc
addpath('../scripts')

%% 1. Load data
% These Argo profiles had a value interpolated to -10 m.

iargo=3;
if iargo==1; load data/D5903264_520.mat; end % cada 2 m
if iargo==2; load data/D1900039_112.mat; end % cada 10, 20 y 50 m
if iargo==3; load data/D1900044_079.mat; end % varía en todo el perfil

rho = flip(rho);
z = flip(z);

%% 2. Compute the Work done by buoyancy

[WB, z_wb] = buoyancy_potential_work(rho, z);


%% 3. Figure

close all
figure(1)
tf = tiledlayout(1,1);
ax1 = axes(tf);
plot(ax1,rho,z,'-b','linewidth',3); grid on; xlabel('\sigma_0 (kg\cdotm^{-3})'); ylabel('Depth (m)');
xlim([min(rho) max(rho)]); ylim([min(z) max(z)])
ax1.XColor = 'b';
ax1.YColor = 'b';

ax2 = axes(tf);
plot(ax2,WB,z_wb,'-r','linewidth',3); xlabel('WB (J\cdotm^{-3})'); ylabel('Depth (m)');
xlim([min(WB) max(WB)]); ylim([min(z) max(z)])
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.XColor = 'r';
ax2.YColor = 'r';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';
