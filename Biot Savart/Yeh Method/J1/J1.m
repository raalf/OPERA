clc
clear

str = fileread('J1_in.txt');

% Warning, the following variable name replacements were made: x__m~ -> cg, xi__1~ -> cg1, xi__3~ -> cg3, y__m~ -> cg5, z__m~ -> cg7
tmp_exp = {'cg15','cg13','cg11','cg5','cg1','cg3',  'cg7','cg9','cg'};
var_exp = {'z_m','y_m','xi_3',   'E',  'D_LE','D_TE', 'x_m','xi_1','C'};

expression = [{'__','/','*','\^'} tmp_exp];
replace = [{'_','./','.*','.^'} var_exp];
str = regexprep(str, expression, replace)

str = strrep(str,char(10),'');  % remove LF characters

fid = fopen('J1_out.txt','wt');
fprintf(fid, str);
fclose(fid);