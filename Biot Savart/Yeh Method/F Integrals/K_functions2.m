clc
clear

fn = 0:7;

for i = 1:length(fn)

fcn = ['K',num2str(fn(i))]; 

str = fileread(['K Integrals/', fcn, '_in.txt']);

k = regexp(str,'[\r]');
tmp_exp = [];
var_exp = [];
body = str;

expression = [{'__','/','*','\^','~ ','~;','~)','~,'} tmp_exp];
replace = [{'_','./','.*','.^',' ',';',')',','} var_exp];
body = regexprep(body, expression, replace);

body = strrep(body,char(10),'');  % remove LF characters
temp = splitlines(body);
temp = temp(~cellfun('isempty',temp));
lastline = temp{end};
k = strfind(lastline,'=');
tnum = lastline(1:k(1)-2);

fcnheader = sprintf('function I = fcn%s_new(S, T, u, alpha, F_1, F_2, tol)\n',fcn);

I = sprintf('\nI = %s;\n',tnum);
fcnfooter = sprintf('\nend\n');

fid = fopen(['G:\GIT\opera\K Functions\fcn', fcn, '_new.m'],'wt');
% fid = fopen(['C:\Users\travi\OneDrive\Desktop\GIT\opera\fcn', fcn, '.m'],'wt');
fprintf(fid, [fcnheader body I fcnfooter]);
fclose(fid);

% delete(['K Integrals/', fcn, '_in.txt'])

end








