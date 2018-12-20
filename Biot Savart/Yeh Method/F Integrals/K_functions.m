clc
clear

% fcn = 'K0'; tnum = 't254';
% fcn = 'K1'; tnum = 't236';
% fcn = 'K2'; tnum = 't254';
% fcn = 'K3'; tnum = 't225';
% fcn = 'K4'; tnum = 't218';
% fcn = 'K5'; tnum = 't254';
% fcn = 'K6'; tnum = 't381';
% fcn = 'K7'; tnum = 't390';

% fcn = 'K0ip'; tnum = 't86';
% fcn = 'K1ip'; tnum = 't67';
% fcn = 'K2ip'; tnum = 't50';
% fcn = 'K3ip'; tnum = 't37';
% fcn = 'K4ip'; tnum = 't17';
% fcn = 'K5ip'; tnum = 't17';
% fcn = 'K6ip'; tnum = 't40';
% fcn = 'K7ip'; tnum = 't55';

% fid = fopen(['K Integrals/', fcn, '_in.txt'],'wt');
% fclose(fid);

str = fileread(['K Integrals/', fcn, '_in.txt']);

k = regexp(str,'[\r]');
tmp_exp = [];
var_exp = [];
body = str;

expression = [{'__','/','*','\^'} tmp_exp];
replace = [{'_','./','.*','.^'} var_exp];
body = regexprep(body, expression, replace);

body = strrep(body,char(10),'');  % remove LF characters

fcnheader = sprintf('function I = fcn%s(S, t, u, alpha, F, tol)\n',fcn);

corr = sprintf('\nF(abs(F(:,:,1)) < tol,:,1) = sign(F(abs(F(:,:,1)) < tol,:,2)).*tol;\nF(abs(F(:,:,2)) < tol,:,2) = sign(F(abs(F(:,:,2)) < tol,:,1)).*tol;\n');

I = sprintf('\nI = %s(:,:,2) - %s(:,:,1);\n',tnum, tnum);
fcnfooter = sprintf('\nend\n');


fid = fopen(['G:\GIT\opera\fcn', fcn, '.m'],'wt');
fprintf(fid, [fcnheader corr body I fcnfooter]);
fclose(fid);











