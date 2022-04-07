
function [s, t] = json_read(fn)

if ~iscell(fn); fn = {fn}; end

for ii=1:length(fn)
    fid = fopen(fn{ii});
    
    t1 = fscanf(fid, '%s');
    if ii==1; t = t1;
    else; t = [t(1:(end-1)), ', ', t1(2:end)];
    end
    
    fclose(fid);
end

s = jsondecode(t);

end