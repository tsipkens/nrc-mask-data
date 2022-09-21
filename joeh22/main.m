

addpath cmap autils;

clear;
close all;
clc;

% Colormap used for some plots.
cm = inferno(300);
cm = [0, 1, 0; cm];  % add green at the beginning of the colormap

% Read data and format.
s = tools.json_read('sipkens22.json');  % read as Matlab structure
s = db.fillnan(s);
d = s(1).dm;

% Get mass-mobility information (requires autils).
prop = massmob.init('salt');  % get mass-mobility parameters

%== INITIAL PROCESSING ===================================================%
% Compute integrated quantities.
n_ni = pfe.npfe([s.smps1], [s.smps2]);  % NPFE
m_ni = pfe.mpfe_ni([s.smps1], [s.smps2], [s.dm], prop);  % MPFE (numerical integration)
m_hc = pfe.mpfe_hc([s.smps1], [s.smps2], [s.dm], prop);  % MPFE (Hatch-Choate)

% Compute scattering-based PFE.
ri = 1.54 + 0j;  % particle refractive index (salt)
int = mie.get_intensity(532e-9, [s.dm] .* 1e-9, ri, [], 45/180*pi);  % use Mie theory (requires autils)
sca_45 = pfe.scapfe_ni([s.smps1], [s.smps2], int);  % use numerical integration to compute PFE

% Get SPFE and PFE at dm = 100 nm.
[sp, ssp] = pfe.spfe([s.smps1], [s.smps2], ...
    [s.se1] ./ sqrt([s.N_SMPS1_xN]), [s.se2] ./ sqrt([s.N_SMPS2_xN]), 1);
pfe_100 = pfe.spfe_d([s.smps1], [s.smps2], [s.dm], 100);

[~, ssp1] = pfe.spfe([s.smps1], [s.smps2], ...
    [s.se1], [s.se2], 1);  % SPFE with different uncertainty intervals than above
dmpps = tools.get_mpps(sp, [s.dm], ssp1);  % get most-penetrating particle size (MPPS)

%{
% Alternative MPPS calc.
[dg_up, sg_up] = get_geo([s.smps1], [s.dm]);  % get geometric moments (upstream)
[dg_down, sg_down] = get_geo([s.smps2], [s.dm]);  % get geometric moments (downstream)
dmpps = tools.hatch2mpps(dg_up, sg_up, dg_down, sg_down);
%}

color = max(log10(min(dmpps, 500)), 1);  % for coloring points below (capped on either end)

% Convert to penetration and copy back to structure.
for ii=1:length(s)
    s(ii).m_hc = 1 - m_hc(ii);
    s(ii).n_ni = 1 - n_ni(ii);
    s(ii).sca = 1 - sca_45(:,ii);
    s(ii).sp = 1 - sp(:,ii);
end


%%
% FIG 2: Plot NPFE v. MPFE. 
figure(2);
scatter(1 - n_ni, 1 - m_hc, 20, color, 'filled');
hold on;
plot([0.01,1], [0.01, 1], 'k');
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');
xlim([1e-3, 1.5]);  ylim([3e-4, 1.5]);
title('NPFE v. MPFE');
xlabel('Number-based pen.');
ylabel('Mass-based pen.');

colormap(cm);
colorbar;

%%
% FIG 3: Compare mass- and scattering-based PFEs.
figure(3);
scatter(1 - m_hc, 1 - sca_45, 20, color, 'filled');
hold on;
plot([0.01,1], [0.01, 1], 'k');
hold off

xlim([1e-3, 1.5]);  ylim([3e-4, 1.5]);
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');
title('MPFE v. Scattering-based PFE');
xlabel('Mass-based pen.');
ylabel('Scattering-based pen.')

colormap(cm);
colorbar;

%%
% FIG 4: SPFE for subset of data with known properties.
figure(4);
fgl = 3:4;  % Greenline
fh = color > 2.1;  % flag high MPPS
fl = color < 1.8;  % flag low MPPS
m1 = median(real(log(1 - sp(:,fh))), 2, 'omitnan');
s1l = prctile(real(log(1 - sp(:,fh))), 25, 2);
s1u = prctile(real(log(1 - sp(:,fh))), 75, 2);
m2 = median(real(log(1 - sp(:,fl))), 2, 'omitnan');
s2l = prctile(real(log(1 - sp(:,fl))), 25, 2);
s2u = prctile(real(log(1 - sp(:,fl))), 75, 2);
mgl = median(real(log(1 - sp(:,fgl))), 2, 'omitnan');

plot(d, exp(m1), 'k');
hold on;
plot(d, exp(s1l), '--');
plot(d, exp(s1u), '--');
plot(d, exp(m2), 'k');
plot(d, exp(s2l), '--');
plot(d, exp(s2u), '--');
plot(d, exp(mgl), 'k-', 'LineWidth', 1.5);  % Greenline
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'normal', 'YDir', 'reverse');
xlim([22, 320]);
title('SPFE for data subset');
legend({'dMPPS > 125 nm','err','err', 'dMPPS < 63 nm','err','err', 'Greenline'});

%%
% FIG 5: SPFE @ 100 nm v. MPFE.
figure(5);
scatter(1 - m_hc, 1 - pfe_100, 20, color, 'filled');
hold on;
plot([0.01,1], [0.01, 1], 'k');
hold off;

% Get correlation between MPFE and SPFE @ 100 nm.
A = 1 - [pfe_100', m_hc'];
A(any(isnan(A), 2), :) = [];
A(any(isinf(A), 2), :) = [];
R2 = corrcoef(A(:, 1), A(:, 2))

set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');
xlim([1e-3, 1.5]);
ylim([3e-4, 1.5]);
title('SPFE @ 100 v. MPFE');
ylabel('Pen. @ 100 nm');
xlabel('Mass-based pen.');

colormap(cm);
colorbar;


%%
% FIG 6: Compare SMPS-based and CPC NPFE.
figure(6);
plot(1 - n_ni, [s.N_CPC2_xAvg] ./ [s.N_CPC1_xAvg], '.');
hold on;
plot([0.01,1], [0.01, 1], 'k');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');

% Get correlation between SMPS-based and CPC NPFE.
A = 1 - [n_ni', [s.N_CPC2_xAvg]' ./ [s.N_CPC1_xAvg]'];
A(any(isnan(A), 2), :) = [];
A(any(isinf(A), 2), :) = [];
R1 = corrcoef(A(:, 1), A(:, 2))

title('CPC v. SMPS NPFE');
ylabel('SMPS number-based pen.');
xlabel('CPC number-based pen.');


%%
% FIG 7: TSI v. PFEMS
figure(7);

s_f = db.combine(s, 'lot', 'm_hc');  % average over the lots

errorbar([s_f.m_hc], [s_f.TSI1], ...
    [s_f.m_hc_min] - [s_f.m_hc], [s_f.m_hc_max] - [s_f.m_hc], 'o', 'horizontal');
hold on;
errorbar([s_f.m_hc], [s_f.TSI2], ...
    [s_f.m_hc_min] - [s_f.m_hc], [s_f.m_hc_max] - [s_f.m_hc], 'o', 'horizontal');
plot([0.01,1], [0.01, 1], '-k');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');

% For Greenline samples (add as a star).
s_green = db.filter(s_f, 'id', 'Green*');  % filter for Greenline
hold on;
plot([s_green.m_hc], [s_green.TSI2], 'kp', 'MarkerSize', 10, 'LineWidth', 2);  % add this point
hold off;

% Compute correlation between quantities. 
A = [[s_f.m_hc]', [s_f.TSI1]'; [s_green.m_hc]', [s_green.TSI2]'];
A(any(isnan(A), 2), :) = [];
A(any(isinf(A), 2), :) = [];
R3 = corrcoef(A(:, 1), A(:, 2))

title('PFEMS MPFE v. TSIs');
xlabel('PFEMS mass-based pen.');
ylabel('TSI mass-based pen.');

