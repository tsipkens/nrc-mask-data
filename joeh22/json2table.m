
% JSON2TABLE  Convert JSON formatted data to table.
%  
%  Also removes any fields containing rays for CSV export.
%  
%  AUTHOR: Timothy Sipkens, 2022-04-07

function s = json2table(fi, fo, fiels)

if ~exist('fo', 'var'); fo = []; end
if ~exist('fiels', 'var'); fiels = []; end

s = tools.json_read(fi);

% Remove fields not in fiels.
fiel = fields(s);
if isempty(fiels); fiels = fiel; end
for ii=1:length(fiel)
    if ~any(strcmp(fiel{ii}, fiels))
        s = rmfield(s, fiel{ii});
    end
end

% Remove fields containing arrays.
fiel = fields(s);
for ii=1:length(fiel)
    a = s.(fiel{ii});
    w = whos('a');

    if ~all(w.size == 1)
        s = rmfield(s, fiel{ii});
    end
end

s = struct2table(s);

if ~isempty(fo)
    writetable(s, fo);
end

end
