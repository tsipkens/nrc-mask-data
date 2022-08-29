
% COMBINE  Combine repeat entries in a struct into a single entry.
%   Replaces existing entry with a mean and adds a field with a std. dev.
%   
%   AUTHOR: Timothy Sipkens

function [s] = combine(s, fiel, fiel2)

n2 = find(strcmp(fiel2, fields(s)));
fiel3 = [fiel2, '_std'];
fiel3a = [fiel2, '_min'];
fiel3b = [fiel2, '_max'];

ii = 1;
while ii <= length(s)
    fl = strcmp(s(ii).(fiel), {s((ii+1):length(s)).(fiel)});

    if any(fl)
        fl_idx = find(fl) + ii;
        
        st = nanstd([s([ii, fl_idx]).(fiel2)]);
        mi = nanmin([s([ii, fl_idx]).(fiel2)]);
        ma = nanmax([s([ii, fl_idx]).(fiel2)]);
        s(ii).(fiel2) = nanmean([s([ii, fl_idx]).(fiel2)]);
        
        s(ii).(fiel3) = st;
        s(ii).(fiel3a) = mi;
        s(ii).(fiel3b) = ma;
        
        s(fl_idx) = [];
    else
        s(ii).(fiel3) = NaN;
        s(ii).(fiel3a) = NaN;
        s(ii).(fiel3b) = NaN;
    end
    
    ii = ii+1;
end

end

