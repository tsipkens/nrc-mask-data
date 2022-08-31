
% JSON_READ  A function to read in JSON files in Matlab.
%  
%  S = json_read(FN) reads the JSON file with file name, FN, given as a
%  string. Outputs a Matlab structure containing the data.
%  
%  AUTHOR: Timothy Sipkens, 2022-08-31

function s = json_read(fn)

fid = fopen(fn, 'r');

txt = fscanf(fid, '%s');
s = jsondecode(txt);

fclose(fid);

end
