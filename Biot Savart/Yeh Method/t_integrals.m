clc
clear


fid = fopen('D_call.txt','wt');

for i = 1:3
    for j = 1:6
        
        str = fileread(['D_Ints/D', num2str(i), num2str(j), '.txt']);
        body = fcnM2M(str);
        body = ['t', num2str(i), num2str(j), ' = ', body(5:end)];
        fprintf(fid, body);        
    end
end



fclose(fid);






































