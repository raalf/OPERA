clc
clear

str = fileread('F_ints.txt');
str = splitlines(str);

fname = 'F_ints_calls.txt';
delete(fname);
for i = 1:length(str) - 1
   
    % Clean normally
    body = fcnM2M(str{i});
    
    % Clean specifically for this style
    body = regexprep(body,'K_','K'); % K_1 => K1
    idx = strfind(body,'K1'); % Adding K0 to last terms
    body = [body(1:idx + 1), ' + ( ', body(idx + 2:end-3), ' ).*K0', body(end-2:end-1),', 2);'];
    idx = strfind(body,'=');
    num = sscanf(body(1:idx),'F%g'); % F1 => F(:,1)
    body = ['F(:,',num2str(num),') = sum(', body(idx+1:end)];
    f{i} = body;
    
    fid = fopen(fname,'at');
    fprintf(fid, '%s\n',f{i});
    fclose(fid);
end