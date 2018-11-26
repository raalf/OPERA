clc
clear

fcn = 'J5'
str = fileread([fcn, '/', fcn, '_in.txt']);

k = regexp(str,'[\r]');
header = str(1:k(1))
body = str(k(1):end);

k = strfind(header, ':') + 1;
var_name = [];
prev = false;
j = 1;
while true
    valid = ~strcmp(header(k),'-') & ~strcmp(header(k),'>') & ~strcmp(header(k),' ') & ~strcmp(header(k),'~') & ~strcmp(header(k),',');
    if valid && k == length(header)
        disp('Big trouble in Little China')
    end
    
    if valid && ~prev
        var_name = [];
        var_name = [var_name header(k)];
        prev = true;
    elseif valid && prev && k < length(header)
        var_name = [var_name header(k)];
        prev = true;
    elseif (~valid && prev) || (k == length(header) && prev)
        var_name = regexprep(var_name, '__', '_');
        names{j} = var_name;
        j = j + 1;
        prev = false;
    end
    
    if k == length(header)
        break;
    end
    k = k + 1;
end
tmp_exp = names(end:-2:2);
var_exp = names(end-1:-2:1);

expression = [{'__','/','*','\^'} tmp_exp];
replace = [{'_','./','.*','.^'} var_exp];
body = regexprep(body, expression, replace)

body = strrep(body,char(10),'');  % remove LF characters

fid = fopen([fcn, '/', fcn, '_out.txt'],'wt');
fprintf(fid, body);
fclose(fid);