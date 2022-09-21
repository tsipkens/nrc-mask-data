

addpath cmap autils;

clear;
close all;
clc;

% Colormap used for some plots.
cm = inferno(300);
cm = [0, 1, 0; cm];  % add green at the beginning of the colormap

% Read data and format.
s = tools.json_read('sipkens22.json');  % read as Matlab structure
s = s(3:4);  % SELECT only Greenline, then match analysis from main.m
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


%== FIG. 4: Size-resolved penetration ================%
mgl = median(real(log(1 - sp)), 2, 'omitnan');

figure(4);
plot(d, exp(mgl), 'k-', 'LineWidth', 1.5);
hold on;
plot(d, 1 - sp);
hold off;
set(gca, 'XScale', 'log', 'YScale', 'linear', ...
'XDir', 'normal', 'YDir', 'normal');
xlim([22, 320]);


%%
%-- Fitting --%
flag_x = and([s(1).dm] > 22, [s(1).dm] < 320);  % assume dms are the same

dmx = [s(1).dm(flag_x)];
fitfun = @(x) exp(x(3)) .* normpdf(log(dmx), x(1), x(2));
xfun = @(x) fitfun(x) - exp(mgl(flag_x));

x0 = [log(155), log(2), 1];
x1 = lsqnonlin(xfun, x0);
dmpps = exp(x1(1))
sp = exp(x1(2))

% Add to Fig. 4.
figure(4);
hold on;
plot(dmx, fitfun(x1), 'g--');
hold off;
