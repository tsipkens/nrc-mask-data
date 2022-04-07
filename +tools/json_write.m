

function [] = json_write(fn, s)

t = jsonencode(s);

fid = fopen(fn, 'w');
fprintf(fid, t);
fclose(fid);

end
