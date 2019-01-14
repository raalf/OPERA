clc
clear

letter = 'T'

fn = 0:7;

for i = 1:length(fn)

fcn = ['K',num2str(fn(i)),letter]; 

str = fileread(['K Integrals/', fcn, '_in.txt']);

k = regexp(str,'[\r]');
tmp_exp = [];
var_exp = [];
body = str;

expression = [{'__','/','*','\^'} tmp_exp];
replace = [{'_','./','.*','.^'} var_exp];
body = regexprep(body, expression, replace);

body = strrep(body,char(10),'');  % remove LF characters
temp = splitlines(body);
temp = temp(~cellfun('isempty',temp));
lastline = temp{end};
k = strfind(lastline,'=');
if length(k) > 1
   disp('Problem in last line'); 
end
tnum = lastline(1:k(1)-2);

fcnheader = sprintf('function I = fcn%s(S, T, u, alpha, F, tol)\n',fcn);

I = sprintf('\nI = %s(:,:,2) - %s(:,:,1);\n',tnum, tnum);
fcnfooter = sprintf('\nend\n');

fid = fopen(['G:\GIT\opera\K Functions\fcn', fcn, '.m'],'wt');
% fid = fopen(['C:\Users\travi\OneDrive\Desktop\GIT\opera\fcn', fcn, '.m'],'wt');
fprintf(fid, [fcnheader body I fcnfooter]);
fclose(fid);

% delete(['K Integrals/', fcn, '_in.txt'])

end








