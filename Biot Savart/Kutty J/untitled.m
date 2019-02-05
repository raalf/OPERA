clc
clear

%%
str = fileread('u.txt');
str = splitlines(str);
strout = fcnM2M(str);

fid = fopen('u_out.txt','wt');
fprintf(fid, sprintf('%s\n',strout{:}));
fclose(fid);

%%
str = fileread('v.txt');
str = splitlines(str);
strout = fcnM2M(str);

fid = fopen('v_out.txt','wt');
fprintf(fid, sprintf('%s\n',strout{:}));
fclose(fid);

%%
str = fileread('w.txt');
str = splitlines(str);
strout = fcnM2M(str);

fid = fopen('w_out.txt','wt');
fprintf(fid, sprintf('%s\n',strout{:}));
fclose(fid);
