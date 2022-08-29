
% HATCH2MPPS  Convert GMD and GSD (e.g., from Hatch-Choate) to MPPS.
%  Uses upstream and downstream quantities. 
%  
%  AUTHOR: Timothy Sipkens, 2021-07-2-

function [mpps, sg_p] = hatch2mpps(dg_up, sg_up, dg_down, sg_down)

sg_up = log(sg_up);  % convert to std. dev. in log-space
sg_down = log(sg_down);

dg_up = log(dg_up);  % convert to std. dev. in log-space
dg_down = log(dg_down);

sg_p = 1 ./ (1 ./ (sg_up .^ 2) - 1 ./ (sg_down .^ 2));  % penetration GSD

mpps = sg_p .* (dg_up ./ sg_up .^ 2 - dg_down ./ sg_down .^ 2);
mpps = exp(mpps);

sg_p = sqrt(sg_p);

% Simplified variant as above suffers from noise.
% Currently overwrites above. 
% Assumes the penetration has a very wide GSD, such that 
% SG_UP and SG_DOWN are effectively the same.
mpps = exp((dg_up + dg_down) ./ 2);
mpps(mpps > 100) = min(mpps) - 3;

end

