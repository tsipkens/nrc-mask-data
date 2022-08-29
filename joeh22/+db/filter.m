
% FILTER  Filter the structure, by field and value.
%  
%  AUTHOR: Timothy Sipkens, 2022-03-15

function s = filter(s, fiel, val)

flag = regexp({s.(fiel)}, regexptranslate('wildcard', val));

for ii=1:length(flag)
    if isempty(flag{ii})
        flag{ii} = 0;
    end
end

flag = cell2mat(flag);
flag = (flag == 1);

s = s(flag);

end
