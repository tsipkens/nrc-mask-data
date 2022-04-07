
% FILLNAN  Fills empty entires in structures with NaN.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-15

function s = fillnan(s)

fiels = fields(s);

for ii=1:length(fiels)
    if isnumeric(s(1).(fiels{ii}))
        for jj=1:length(s)
            if isempty(s(jj).(fiels{ii}))
                s(jj).(fiels{ii}) = NaN;
            end
        end
    end
end

end
