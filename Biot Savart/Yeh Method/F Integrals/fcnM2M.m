function body = fcnM2M(body)
expression = {'__','/','*','\^'};
replace = {'_','./','.*','.^'};
body = regexprep(body, expression, replace);
body = strrep(body,char(10),'');  % remove LF characters
body = body(find(~isspace(body)));
end