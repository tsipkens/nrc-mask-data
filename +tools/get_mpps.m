
% GET_MPPS  Compute the most penetrating particle size from SRFE curve.
%  
%  NOTE: Designed for use with SMPS measurments, where results are reliable
%  between dmin and dmax mobility diameter.
%  
%  AUTHOR: Timothy Sipkens, 2021-07-21

function [dmpps, varmpps, pmpps] = get_mpps(eta, d, s)

disp('Estimating MPPS:');
tools.textbar([0,size(eta, 2)]);

dmax = 250;
dmin = 25;

P = 1 - eta;

% Loop over sets of size-resolved measurements.
for ii=1:size(eta, 2)

    di = d(:, ii);
    Pi = P(:, ii);
    si = s(:, ii);

    mu0 = Pi;
    sig0 = si;

    Pi = mu0;
    si = sig0;

    %-{
    % Remove faulty values.
    % Set so as to be very certain that it is not these values.
    f_null = or(or(Pi == 0, isnan(Pi)), si > Pi);
    si(isnan(Pi)) = NaN;
    Pi(isnan(Pi)) = NaN;

    diffm = (Pi - Pi');
    diffs = sqrt(si .^ 2 + (si .^ 2)');

    % Compute probabilities of bins being MPPS.
    pr_ij = 0.5 .* (1 - erf(diffm ./ (sqrt(2) .* diffs))) + eps;
    pr_i = -sum(log(pr_ij), 2, 'omitnan')';

    % Eliminated nulls flagged above.
    pr_i(f_null) = -inf;

    % Add filter to ignore often stray 
    % measurements below dmin and above dmax. 
    pr_i(di < dmin) = -inf;
    pr_i(di > dmax) = -inf;

    % Account for probability peak/PDF width, 
    % as a shortcut around full numerical integration.
    pr_i = pr_i - log(si)';
    pr_i = pr_i - max(pr_i);  % scale distribution
    pr_i = exp(pr_i);  % apply exponential
    
    % Filter based on local fluctuations.
    pr_i = pr_i .* (movstd(Pi', 6) < 500);

    % Scale pr_ij.
    pr_i = pr_i ./ nansum(pr_i);
    pr_i(isnan(pr_i)) = 0;
    %}

    dmpps(ii) = exp(nansum(pr_i .* log(di)'));
    gvarmpps(ii) = exp(sum(pr_i .* (log(di)' - log(dmpps(ii))) .^ 2));
    varmpps(:, ii) = [(gvarmpps(ii) - 1) .* ...
        exp(2 .* log(dmpps(ii)) + log(gvarmpps(ii))), 0];

    [~,idx] = min(abs(log(di) - log(dmpps(ii))));
    pmpps(ii) = min(Pi(idx), 100);  % penetration at MPPS
    
    %-{
    % Alternative intervals.
    % Get interval around dmpps when penetration is similar.
    t0 = find((Pi + si) < pmpps(ii));
    if any(t0 < idx)
        t1 = di(max(t0(t0 < idx)));
    else
        t1 = 1;  % goes to min.
    end
    if any(t0 > idx)
        t2 = di(min(t0(t0 > idx)));
    else
        t2 = 1e3;  % goes to max.
    end
    varmpps(:, ii) = abs([t1, t2] - dmpps(ii)) .^ 2;
    %}

    %{
    %== DIAGNOSTIC PLOT ==============================================%
    figure(gcf);
    plot(di, Pi, '.-');
    hold on;
    plot(di, Pi - si, '.-');
    plot(di, Pi + si, '.-');
    [~, im] = min((di - dmpps(ii)) .^ 2, [], 'omitnan');
    plot(dmpps(ii), Pi(im), 'ok', 'MarkerSize', 30);
    errorbar(dmpps(ii), Pi(im), ...
        sqrt(varmpps(1, ii)), sqrt(varmpps(2, ii)), 'horizontal', 'k');
    hold off;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    ylim([-inf, 1.2]);
    pause(0.5);
    %=================================================================%
    %}
    
    if dmpps(ii) == 1  % MPPS not found, set to min. - 5
        dmpps(ii) = dmin - 2;
    end
    if dmpps(ii) == dmin  % MPPS not found, set to min. - 5
        dmpps(ii) = dmin - 2;
    end
    
    tools.textbar([ii,size(eta, 2)]);
end

end
