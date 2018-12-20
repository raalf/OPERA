clc
clear

str{1} = '(2 * L * S + L) * F ^ 3 + (6 * N * S - 3 * N) * F ^ 2 + (-3 * L * alpha + 6 * L * u) * F + N * alpha + 2 * u * N';
str{2} = '(2 * L * S + L) * F ^ 4 + (6 * N * S - 3 * N - (2 * L * S + L) * x_m) * F ^ 3 + (-3 * L * alpha + 6 * L * u - (6 * N * S - 3 * N) * x_m) * F ^ 2 + (N * alpha + 2 * u * N - (-3 * L * alpha + 6 * L * u) * x_m) * F - (N * alpha + 2 * u * N) * x_m';
str{3} = '(2 * L * S + L) * F ^ 5 + (6 * N * S - 3 * N - 2 * (2 * L * S + L) * x_m) * F ^ 4 + (-3 * L * alpha + 6 * L * u - 2 * (6 * N * S - 3 * N) * x_m + (2 * L * S + L) * x_m ^ 2) * F ^ 3 + (N * alpha + 2 * u * N - 2 * (-3 * L * alpha + 6 * L * u) * x_m + (6 * N * S - 3 * N) * x_m ^ 2) * F ^ 2 + (-2 * (N * alpha + 2 * u * N) * x_m + (-3 * L * alpha + 6 * L * u) * x_m ^ 2) * F + (N * alpha + 2 * u * N) * x_m ^ 2';
str{4} = 'F ^ 4 + (2 * L * S * y_m + L * y_m) * F ^ 3 + (6 * N * S * y_m - 3 * N * y_m + 2 * alpha) * F ^ 2 + (-3 * L * alpha * y_m + 6 * L * u * y_m) * F + N * y_m * alpha + 2 * y_m * u * N + alpha ^ 2';
str{5} = '(L * S - L) * F ^ 5 + (3 * N * S - 3 * N + 2 * y_m) * F ^ 4 + (2 * L * S * y_m ^ 2 + L * S * alpha + L * y_m ^ 2 - 4 * L * alpha + 3 * L * u) * F ^ 3 + (6 * N * S * y_m ^ 2 + 3 * N * S * alpha - 3 * N * y_m ^ 2 - 4 * N * alpha + N * u + 4 * y_m * alpha) * F ^ 2 + (-3 * L * alpha * y_m ^ 2 + 6 * L * u * y_m ^ 2 - 3 * L * alpha ^ 2 + 3 * L * alpha * u) * F + N * alpha * y_m ^ 2 + 2 * N * u * y_m ^ 2 - N * alpha ^ 2 + N * alpha * u + 2 * y_m * alpha ^ 2';
str{6} = 'F ^ 5 + (2 * L * S * y_m + y_m * L - x_m) * F ^ 4 + (6 * N * y_m * S - 3 * N * y_m + 2 * alpha - (2 * L * S * y_m + y_m * L) * x_m) * F ^ 3 + (-3 * L * alpha * y_m + 6 * L * u * y_m - (6 * N * y_m * S - 3 * N * y_m + 2 * alpha) * x_m) * F ^ 2 + (N * y_m * alpha + 2 * y_m * u * N + alpha ^ 2 - (-3 * L * alpha * y_m + 6 * L * u * y_m) * x_m) * F - (N * y_m * alpha + 2 * y_m * u * N + alpha ^ 2) * x_m';
str{7} = '(2 * L * S + L) * F ^ 6 + (6 * N * S - 3 * N - 3 * (2 * L * S + L) * x_m) * F ^ 5 + (-3 * L * alpha + 6 * L * u - 3 * (6 * N * S - 3 * N) * x_m + 3 * (2 * L * S + L) * x_m ^ 2) * F ^ 4 + (N * alpha + 2 * u * N - 3 * (-3 * L * alpha + 6 * L * u) * x_m + 3 * (6 * N * S - 3 * N) * x_m ^ 2 - (2 * L * S + L) * x_m ^ 3) * F ^ 3 + (-3 * (N * alpha + 2 * u * N) * x_m + 3 * (-3 * L * alpha + 6 * L * u) * x_m ^ 2 - (6 * N * S - 3 * N) * x_m ^ 3) * F ^ 2 + (3 * (N * alpha + 2 * u * N) * x_m ^ 2 - (-3 * L * alpha + 6 * L * u) * x_m ^ 3) * F - (N * alpha + 2 * u * N) * x_m ^ 3';
str{8} = '(3 * S - 1) * F ^ 6 + (3 * L * S * y_m + 6 * L * N - 3 * L * y_m) * F ^ 5 + (9 * N * S * y_m - 9 * N * y_m + 6 * S * alpha + 3 * y_m ^ 2 - 3 * alpha + 3 * u) * F ^ 4 + (2 * L * S * y_m ^ 3 + 3 * L * S * alpha * y_m + L * y_m ^ 3 + 12 * L * N * alpha - 12 * L * alpha * y_m + 9 * L * u * y_m) * F ^ 3 + (6 * N * S * y_m ^ 3 + 9 * N * S * alpha * y_m - 3 * N * y_m ^ 3 - 12 * N * alpha * y_m + 3 * N * u * y_m + 3 * S * alpha ^ 2 + 6 * alpha * y_m ^ 2 - 3 * alpha ^ 2 + 6 * alpha * u) * F ^ 2 + (-3 * L * alpha * y_m ^ 3 + 6 * L * u * y_m ^ 3 + 6 * L * N * alpha ^ 2 - 9 * L * alpha ^ 2 * y_m + 9 * L * alpha * u * y_m) * F + N * alpha * y_m ^ 3 + 2 * N * u * y_m ^ 3 - 3 * N * alpha ^ 2 * y_m + 3 * N * alpha * u * y_m + 3 * alpha ^ 2 * y_m ^ 2 - alpha ^ 3 + 3 * alpha ^ 2 * u';
str{9} = 'F ^ 6 + (2 * L * S * y_m + y_m * L - 2 * x_m) * F ^ 5 + (x_m ^ 2 - 2 * x_m * (2 * L * S * y_m + y_m * L) + 6 * N * y_m * S - 3 * N * y_m + 2 * alpha) * F ^ 4 + (x_m ^ 2 * (2 * L * S * y_m + y_m * L) - 2 * x_m * (6 * N * y_m * S - 3 * N * y_m + 2 * alpha) - 3 * L * alpha * y_m + 6 * L * u * y_m) * F ^ 3 + (x_m ^ 2 * (6 * N * y_m * S - 3 * N * y_m + 2 * alpha) - 2 * x_m * (-3 * L * alpha * y_m + 6 * L * u * y_m) + N * y_m * alpha + 2 * y_m * u * N + alpha ^ 2) * F ^ 2 + (x_m ^ 2 * (-3 * L * alpha * y_m + 6 * L * u * y_m) - 2 * x_m * (N * y_m * alpha + 2 * y_m * u * N + alpha ^ 2)) * F + x_m ^ 2 * (N * y_m * alpha + 2 * y_m * u * N + alpha ^ 2)';
str{10} = '(L * S - L) * F ^ 6 + (-x_m * (L * S - L) + 3 * N * S - 3 * N + 2 * y_m) * F ^ 5 + (-x_m * (3 * N * S - 3 * N + 2 * y_m) + 2 * L * S * y_m ^ 2 + L * S * alpha + L * y_m ^ 2 - 4 * L * alpha + 3 * L * u) * F ^ 4 + (-x_m * (2 * L * S * y_m ^ 2 + L * S * alpha + L * y_m ^ 2 - 4 * L * alpha + 3 * L * u) + 6 * N * S * y_m ^ 2 + 3 * N * S * alpha - 3 * N * y_m ^ 2 - 4 * N * alpha + N * u + 4 * y_m * alpha) * F ^ 3 + (-x_m * (6 * N * S * y_m ^ 2 + 3 * N * S * alpha - 3 * N * y_m ^ 2 - 4 * N * alpha + N * u + 4 * y_m * alpha) - 3 * L * alpha * y_m ^ 2 + 6 * L * u * y_m ^ 2 - 3 * L * alpha ^ 2 + 3 * L * alpha * u) * F ^ 2 + (-x_m * (-3 * L * alpha * y_m ^ 2 + 6 * L * u * y_m ^ 2 - 3 * L * alpha ^ 2 + 3 * L * alpha * u) + N * alpha * y_m ^ 2 + 2 * N * u * y_m ^ 2 - N * alpha ^ 2 + N * alpha * u + 2 * y_m * alpha ^ 2) * F - x_m * (N * alpha * y_m ^ 2 + 2 * N * u * y_m ^ 2 - N * alpha ^ 2 + N * alpha * u + 2 * y_m * alpha ^ 2)';
str{11} = '(2 * L * S + L) * F ^ 7 + (-4 * x_m * (2 * L * S + L) + 6 * N * S - 3 * N) * F ^ 6 + (6 * x_m ^ 2 * (2 * L * S + L) - 4 * x_m * (6 * N * S - 3 * N) - 3 * L * alpha + 6 * L * u) * F ^ 5 + (-4 * x_m ^ 3 * (2 * L * S + L) + 6 * x_m ^ 2 * (6 * N * S - 3 * N) - 4 * x_m * (-3 * L * alpha + 6 * L * u) + N * alpha + 2 * u * N) * F ^ 4 + (x_m ^ 4 * (2 * L * S + L) - 4 * x_m ^ 3 * (6 * N * S - 3 * N) + 6 * x_m ^ 2 * (-3 * L * alpha + 6 * L * u) - 4 * x_m * (N * alpha + 2 * u * N)) * F ^ 3 + (x_m ^ 4 * (6 * N * S - 3 * N) - 4 * x_m ^ 3 * (-3 * L * alpha + 6 * L * u) + 6 * x_m ^ 2 * (N * alpha + 2 * u * N)) * F ^ 2 + (x_m ^ 4 * (-3 * L * alpha + 6 * L * u) - 4 * x_m ^ 3 * (N * alpha + 2 * u * N)) * F + x_m ^ 4 * (N * alpha + 2 * u * N)';
str{12} = '(L * S - L) * F ^ 7 + (3 * N * S - 3 * N + 2 * y_m - 2 * (L * S - L) * x_m) * F ^ 6 + (2 * L * S * y_m ^ 2 + L * S * alpha + L * y_m ^ 2 - 4 * L * alpha + 3 * L * u - 2 * (3 * N * S - 3 * N + 2 * y_m) * x_m + (L * S - L) * x_m ^ 2) * F ^ 5 + (6 * N * S * y_m ^ 2 + 3 * N * S * alpha - 3 * N * y_m ^ 2 - 4 * N * alpha + u * N + 4 * alpha * y_m - 2 * (2 * L * S * y_m ^ 2 + L * S * alpha + L * y_m ^ 2 - 4 * L * alpha + 3 * L * u) * x_m + (3 * N * S - 3 * N + 2 * y_m) * x_m ^ 2) * F ^ 4 + (-3 * L * alpha * y_m ^ 2 + 6 * L * u * y_m ^ 2 - 3 * L * alpha ^ 2 + 3 * L * alpha * u - 2 * (6 * N * S * y_m ^ 2 + 3 * N * S * alpha - 3 * N * y_m ^ 2 - 4 * N * alpha + u * N + 4 * alpha * y_m) * x_m + (2 * L * S * y_m ^ 2 + L * S * alpha + L * y_m ^ 2 - 4 * L * alpha + 3 * L * u) * x_m ^ 2) * F ^ 3 + (N * alpha * y_m ^ 2 + 2 * N * u * y_m ^ 2 - N * alpha ^ 2 + N * alpha * u + 2 * y_m * alpha ^ 2 - 2 * (-3 * L * alpha * y_m ^ 2 + 6 * L * u * y_m ^ 2 - 3 * L * alpha ^ 2 + 3 * L * alpha * u) * x_m + (6 * N * S * y_m ^ 2 + 3 * N * S * alpha - 3 * N * y_m ^ 2 - 4 * N * alpha + u * N + 4 * alpha * y_m) * x_m ^ 2) * F ^ 2 + (-2 * (N * alpha * y_m ^ 2 + 2 * N * u * y_m ^ 2 - N * alpha ^ 2 + N * alpha * u + 2 * y_m * alpha ^ 2) * x_m + (-3 * L * alpha * y_m ^ 2 + 6 * L * u * y_m ^ 2 - 3 * L * alpha ^ 2 + 3 * L * alpha * u) * x_m ^ 2) * F + (N * alpha * y_m ^ 2 + 2 * N * u * y_m ^ 2 - N * alpha ^ 2 + N * alpha * u + 2 * y_m * alpha ^ 2) * x_m ^ 2';
str{13} = '(4 * L * S - L) * F ^ 7 + (12 * N * S - 12 * y_m * S - 9 * N + 4 * y_m) * F ^ 6 + (-6 * L * S * y_m ^ 2 - 24 * L * N * y_m + 8 * L * S * alpha + 6 * L * y_m ^ 2 - 11 * L * alpha + 12 * L * u) * F ^ 5 + (-18 * N * y_m ^ 2 * S + 24 * N * alpha * S + 18 * N * y_m ^ 2 - 24 * S * alpha * y_m - 4 * y_m ^ 3 - 19 * N * alpha + 4 * N * u + 12 * y_m * alpha - 12 * u * y_m) * F ^ 4 + (-2 * L * S * y_m ^ 4 - 6 * L * S * alpha * y_m ^ 2 - L * y_m ^ 4 - 48 * L * N * alpha * y_m + 4 * L * S * alpha ^ 2 + 24 * L * alpha * y_m ^ 2 - 18 * L * u * y_m ^ 2 - 19 * L * alpha ^ 2 + 24 * L * alpha * u) * F ^ 3 + (-6 * N * y_m ^ 4 * S - 18 * N * y_m ^ 2 * alpha * S + 3 * N * y_m ^ 4 + 12 * N * alpha ^ 2 * S + 24 * N * y_m ^ 2 * alpha - 6 * N * y_m ^ 2 * u - 12 * S * alpha ^ 2 * y_m - 8 * alpha * y_m ^ 3 - 11 * N * alpha ^ 2 + 8 * N * alpha * u + 12 * y_m * alpha ^ 2 - 24 * alpha * u * y_m) * F ^ 2 + (3 * L * alpha * y_m ^ 4 - 6 * L * u * y_m ^ 4 - 24 * L * N * alpha ^ 2 * y_m + 18 * L * alpha ^ 2 * y_m ^ 2 - 18 * L * alpha * u * y_m ^ 2 - 9 * L * alpha ^ 3 + 12 * L * alpha ^ 2 * u) * F - N * alpha * y_m ^ 4 - 2 * N * u * y_m ^ 4 + 6 * N * alpha ^ 2 * y_m ^ 2 - 6 * N * alpha * u * y_m ^ 2 - 4 * alpha ^ 2 * y_m ^ 3 - N * alpha ^ 3 + 4 * N * alpha ^ 2 * u + 4 * y_m * alpha ^ 3 - 12 * alpha ^ 2 * u * y_m';

df(1,:) = [1 -1];
df(2,:) = [-1 1];
df(3,:) = [1 -1];
df(4,:) = [1 -1];
df(5,:) = [1 -1];
df(6,:) = [-1 1];
df(7,:) = [-1 1];
df(8,:) = [1 -1];
df(9,:) = [1 -1];
df(10,:) = [-1 1];
df(11,:) = [1 -1];
df(12,:) = [1 -1];
df(13,:) = [-1 1];

fid = fopen('FCalls.txt', 'wt');
fclose(fid);
for i = 1:size(str,2)
    body = fcnM2M(str{i});
    body = regexprep(body, {'-','+'}, {' - ', ' + '});
    body = regexprep(body, {'F.\^'}, {'K'});
    idx = strfind(body,'F'); 
    body = [body(1:idx + 1), ' + ( ', body(idx + 2:end), ' ).*K0'];
    body = regexprep(body, {'F'}, {'K1'});
    cmnd = ['F(:,', num2str(i), ') = sum([', num2str(df(i,1)),'*(1/3) ', num2str(df(i,2)), '*(1/3)].*(', body, '), 2);'];
    fid = fopen('FCalls.txt', 'a+');
    fprintf(fid, '%s\n', cmnd);
    fclose(fid);
end


body = fcnM2M(body);
body = regexprep(body, {'-','+'}, {' - ', ' + '});
body = regexprep(body, {'F.\^'}, {'K'});








