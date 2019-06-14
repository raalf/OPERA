function body = fcnM2M(body)
expression = {'__','/','*','\^','~'};
replace = {'_','./','.*','.^',''};
body = regexprep(body, expression, replace);

expression = {'C','D_LE','x_m','y_m','z_m','eta_2','eta_3','xi_2','xi_3'};
replace = {'C(vecBI)','D_LE(vecBI)','x_m(vecBI)','y_m(vecBI)','z_m(vecBI)','eta_2(vecBI)','eta_3(vecBI)','xi_2(vecBI)','xi_3(vecBI)'};
body = regexprep(body, expression, replace);

body = strrep(body,char(10),'');  % remove LF characters
% body = body(find(~isspace(body)));
end