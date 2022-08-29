

addpath cmap autils;

clear;
close all;
clc;


% Do some data analysis. 
s = tools.json_read('sipkens22.json');
s = db.fillnan(s);
d = s(1).dm;

% s = s([359, 134, 371, 643, 649, 297, 622, 643]);


prop = massmob.init('salt');  % get mass-mobility parameters
prop_alt = massmob.init('zet', 2.97, 'rho100', 1940);

% Compute integrated quantities.
[dg_up, sg_up] = get_geo([s.smps1], [s.dm]);
[dg_down, sg_down] = get_geo([s.smps2], [s.dm]);
n_ni = pfe.npfe([s.smps1], [s.smps2]);
m_ni = pfe.mpfe_ni([s.smps1], [s.smps2], [s.dm], prop);
m_hc = pfe.mpfe_hc([s.smps1], [s.smps2], [s.dm], prop);
m_hc2 = pfe.mpfe_hc([s.smps1], [s.smps2], [s.dm], prop_alt);

ri = 1.54 + 0j;  % particle refractive index (salt)
int = mie.get_intensity(532e-9, [s.dm] .* 1e-9, ri, [], 45/180*pi);
sca_45 = pfe.scapfe_ni([s.smps1], [s.smps2], int);

[sp, ssp] = pfe.spfe([s.smps1], [s.smps2], ...
    [s.se1] ./ sqrt([s.N_SMPS1_xN]), [s.se2] ./ sqrt([s.N_SMPS2_xN]), 1);
pfe_100 = pfe.spfe_d([s.smps1], [s.smps2], [s.dm], 100);

[~, ssp1] = pfe.spfe([s.smps1], [s.smps2], ...
    [s.se1], [s.se2], 1);
dmpps = tools.get_mpps(sp, [s.dm], ssp1);
% dmpps = tools.hatch2mpps(dg_up, sg_up, dg_down, sg_down);


for ii=1:length(s)
    s(ii).m_hc = 1 - m_hc(ii);
    s(ii).n_ni = 1 - n_ni(ii);
    s(ii).sca = 1 - sca_45(:,ii);
    s(ii).sp = 1 - sp(:,ii);
end


% Percent error relative to zet = 3.
pe_zet = (-m_hc2 + m_hc) ./ (1 - m_hc);
figure(1);
subplot(3, 1, [2,3])
plot(1 - n_ni, 1 - m_hc, '.');
hold on;
plot(1 - n_ni, 1 - m_hc2, 'o', 'MarkerSize', 4);
plot([0.01,1], [0.01, 1]);
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');
xlabel('Number-based pen.');
ylabel('Mass-based pen.');
lim_x = xlim();

subplot(3, 1, 1);
plot(1 - m_hc, 1 - pe_zet, '.');
set(gca, 'XScale', 'log', 'XDir', 'reverse');
title('Effect of mass-mobility exponent');
ylabel('Fractional difference');
xlim(lim_x);


%%
% Plot NPFE v. MPFE. 
figure(2);
cl = max(log10(min(dmpps, 500)), 1);
scatter(1 - n_ni, 1 - m_hc, 20, cl, 'filled');
hold on;
% plot(1 - n_ni, 1 - m_ni, 'ok', 'MarkerSize', 4)
% plot(1 - n_ni, 1 - m_hc2, 'ok', 'MarkerSize', 4);  % for alternate zet
plot([0.01,1], [0.01, 1], 'k');
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');
xlim([1e-3, 1.5]);
ylim([3e-4, 1.5]);
title('NPFE v. MPFE');
xlabel('Number-based pen.');
ylabel('Mass-based pen.');

cm = inferno(300);
cm = [0,1,0;cm];
colormap(cm);
colorbar;


% Compare mass- and scattering-based PFEs.
figure(3);
scatter(1 - m_hc, 1 - sca_45, 20, cl, 'filled');
hold on;
plot([0.01,1], [0.01, 1], 'k');
hold off

xlim([1e-3, 1.5]);
ylim([3e-4, 1.5]);
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');
title('Scattering and mass-based PFEs');

cm = inferno(300);
cm = [0,1,0;cm];
colormap(cm);
colorbar;


%%
figure(4);
f1 = cl > 2.1;
f2 = cl < 1.8;
m1 = median(real(log(1 - sp(:,f1))), 2, 'omitnan');
s1l = prctile(real(log(1 - sp(:,f1))), 25, 2);
s1u = prctile(real(log(1 - sp(:,f1))), 75, 2);
m2 = median(real(log(1 - sp(:,f2))), 2, 'omitnan');
s2l = prctile(real(log(1 - sp(:,f2))), 25, 2);
s2u = prctile(real(log(1 - sp(:,f2))), 75, 2);

plot(d, exp(m1), 'k');
hold on;
plot(d, exp(s1l), '--');
plot(d, exp(s1u), '--');
plot(d, exp(m2), 'k');
plot(d, exp(s2l), '--');
plot(d, exp(s2u), '--');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'normal', 'YDir', 'reverse');
xlim([22, 320]);


figure(6);
scatter(1 - m_hc, 1 - pfe_100, 20, cl, 'filled');
hold on;
plot([0.01,1], [0.01, 1], 'k');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');

A = 1 - [pfe_100', m_hc'];
A(any(isnan(A), 2), :) = [];
A(any(isinf(A), 2), :) = [];
R2 = corrcoef(A(:, 1), A(:, 2))

xlim([1e-3, 1.5]);
ylim([3e-4, 1.5]);
title('PFE @ 100 v. MPFE');
ylabel('Pen. @ 100 nm');
xlabel('Mass-based pen.');

cm = inferno(300);
cm = [0,1,0;cm];
colormap(cm);
colorbar;


%%

figure(5);
plot(1 - n_ni, [s.N_CPC2_xAvg] ./ [s.N_CPC1_xAvg], '.');
hold on;
plot([0.01,1], [0.01, 1], 'k');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');

A = 1 - [n_ni', [s.N_CPC2_xAvg]' ./ [s.N_CPC1_xAvg]'];
A(any(isnan(A), 2), :) = [];
A(any(isinf(A), 2), :) = [];
R1 = corrcoef(A(:, 1), A(:, 2))

title('CPC v. SMPS NPFE');
xlabel('SMPS number-based pen.');
xlabel('CPC number-based pen.');



%%
% TSI v. PFEMS
sin_f = db.combine(s, 'lot', 'm_hc');

figure(7);

errorbar([sin_f.m_hc], [sin_f.TSI1], ...
    [sin_f.m_hc_min] - [sin_f.m_hc], [sin_f.m_hc_max] - [sin_f.m_hc], '.', 'horizontal');
hold on;
errorbar([sin_f.m_hc], [sin_f.TSI2], ...
    [sin_f.m_hc_min] - [sin_f.m_hc], [sin_f.m_hc_max] - [sin_f.m_hc], '.', 'horizontal');
% plot(m_hc, [sin.TSI1], '.b');
% plot(m_hc, [sin.TSI2], '.r');
plot([0.01,1], [0.01, 1], '-k');
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');

sin_f2 = db.filter(sin_f, 'id', 'Green*');
hold on;
% plot([sin_f.m_hc], [sin_f.TSI1], 'or');
% plot([sin_f.m_hc], [sin_f.TSI2], 'ob');
plot([sin_f2.m_hc], [sin_f2.TSI2], 'kp', 'MarkerSize', 10);
hold off;

A = [[sin_f.m_hc]', [sin_f.TSI1]'; [sin_f2.m_hc]', [sin_f2.TSI2]'];
A(any(isnan(A), 2), :) = [];
A(any(isinf(A), 2), :) = [];
R3 = corrcoef(A(:, 1), A(:, 2))

title('PFEMS MPFE v. TSIs');
xlabel('PFEMS mass-based pen.');
ylabel('TSI mass-based pen.');


%{
%-- Weight/charge state specific plots -----%
figure(10);
clf;
vec1 = {'H', 'MH', 'ML', 'L'};
vec2 = {'H', 'M', 'N'};
sym = {'o', '<', '>', 's'};
cl2 = {'g', 'r', 'b'};
hold on;
for ii=1:length(vec1)
    for jj=1:length(vec2)
        sij = db.filter(s, 'lot', ['R.', vec1{ii}, '.', vec2{jj}, '*']);
        plot([sij.n_ni], [sij.m_hc], [sym{ii}, cl2{jj}]);
    end
end
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XDir', 'reverse', 'YDir', 'reverse');
hold on;
plot([0.1, 1], [0.1, 1], 'k-');
hold off;


figure(11);
clf;
hold on;
for ii=1:length(vec2)
    sij = db.filter(s, 'lot', ['R.*.', vec2{ii}, '.*']);
    fd = and(~isnan(prctile([sij.sp], 25, 2)), ...
        ~isnan(prctile([sij.sp], 75, 2)));
    
    p25 = prctile([sij.sp], 25, 2);
    p75 = prctile([sij.sp], 75, 2);
    area([sij(1).dm(fd), fliplr(sij(1).dm(fd))], ...
        [p25(fd), fliplr(p75(fd) - p25(fd))]);
    
    sij = db.filter(s, 'lot', ['R.*.', vec2{ii}, '.*']);
    plot(sij(1).dm, nanmean([sij.sp], 2), [cl2{ii}, '-']);
end
sgl = db.filter(s, 'lot', 'Greenline');
sgl = sgl([2,3]);
plot(sij(1).dm, nanmean([sgl.sp], 2), 'k-')
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log', 'YDir', 'reverse');
ylim([4e-2, 1]);
xlim([22, 320]);
%-----------------------------------------------------%
%}
