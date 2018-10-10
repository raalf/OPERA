clc
clear

str = fileread('J51ip_in.txt');
% C~ -> cg, D__LE~ -> cg1, D__TE~ -> cg3, E~ -> cg5, eta__1~ -> cg7, eta__2~ -> cg9, eta__3~ -> cg11, x__m~ -> cg13, xi__1~ -> cg15, y__m~ -> cg17
tmp_exp = {'cg17','cg15','cg13','cg11','cg9','cg7','cg5','cg3','cg1','cg'};
var_exp = {'y_m','xi_1','x_m','eta_3','eta_2','eta_1','E','D_TE','D_LE','C'};

expression = [{'__','/','*','\^','CodeGeneration:-'} tmp_exp];
replace = [{'_','./','.*','.^',''} var_exp];
str = regexprep(str, expression, replace)

str = strrep(str,char(10),'');  % remove LF characters

fid = fopen('J51ip_out.txt','wt');
fprintf(fid, str);
fclose(fid);