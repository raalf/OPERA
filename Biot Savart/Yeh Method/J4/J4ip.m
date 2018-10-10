clc
clear

str = fileread('J4ip_in.txt');
% Warning, the following variable name replacements were made:...
% C~ -> cg, D__LE~ -> cg1, D__TE~ -> cg3, E~ -> cg5, x__m~ -> cg7, xi__1~ -> cg9, xi__3~ -> cg11, y__m~ -> cg13

tmp_exp = {'cg13','cg11','cg5','cg1','cg3',  'cg7','cg9','cg'};
var_exp = {'y_m','xi_3',   'E',  'D_LE','D_TE', 'x_m','xi_1','C'};

expression = [{'__','/','*','\^'} tmp_exp];
replace = [{'_','./','.*','.^'} var_exp];
str = regexprep(str, expression, replace)

str = strrep(str,char(10),'');  % remove LF characters

fid = fopen('J4ip_out.txt','wt');
fprintf(fid, str);
fclose(fid);