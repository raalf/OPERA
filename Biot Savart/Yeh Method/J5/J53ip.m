clc
clear

str = fileread('J53ip_in.txt');

% D__LE~ -> cg, D__TE~ -> cg1, E~ -> cg3, eta__1~ -> cg5, eta__2~ -> cg7, eta__3~ -> cg9, x__m~ -> cg11, xi__1~ -> cg13, y__m~ -> cg15
tmp_exp = {'cg15','cg13','cg11','cg9','cg7','cg5','cg3','cg1','cg'};
var_exp = {'y_m','xi_1','x_m','eta_3','eta_2','eta_1','E','D_TE','D_LE'};

expression = [{'__','/','*','\^','CodeGeneration:-'} tmp_exp];
replace = [{'_','./','.*','.^',''} var_exp];
str = regexprep(str, expression, replace)

str = strrep(str,char(10),'');  % remove LF characters

fid = fopen('J53ip_out.txt','wt');
fprintf(fid, str);
fclose(fid);