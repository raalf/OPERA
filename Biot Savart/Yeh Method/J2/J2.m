clc
clear

str = fileread('J2_in.txt');

%Warning, the following variable name replacements were made: x__m~ -> cg, xi__1~ -> cg1, y__m~ -> cg3, z__m~ -> cg5
tmp_exp = {'cg3','cg5','cg1','cg'};
var_exp = {'y_m','z_m','xi_1','x_m'};

expression = [{'__','/','*','\^'} tmp_exp];
replace = [{'_','./','.*','.^'} var_exp];
str = regexprep(str, expression, replace)
str = strrep(str,char(10),'');  % remove LF characters

fid = fopen('J2_out.txt','wt');
fprintf(fid, str);
fclose(fid);