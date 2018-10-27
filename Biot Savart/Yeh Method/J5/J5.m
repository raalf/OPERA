clc
clear

str = fileread('J5_in.txt');

% C~ -> cg, D__LE~ -> cg1, D__TE~ -> cg3, E~ -> cg5, x__m~ -> cg7, xi__1~ -> cg9, xi__3~ -> cg11, y__m~ -> cg13, z__m~ -> cg15
tmp_exp = {'cg15','cg13','cg11','cg5','cg1','cg3',  'cg7','cg9','cg'};
var_exp = {'z_m','y_m','xi_3',   'E',  'D_LE','D_TE', 'x_m','xi_1','C'};

expression = [{'__','/','*','\^'} tmp_exp];
replace = [{'_','./','.*','.^'} var_exp];
str = regexprep(str, expression, replace)

str = strrep(str,char(10),'');  % remove LF characters

fid = fopen('J5_out.txt','wt');
fprintf(fid, str);
fclose(fid);