function body = fcnM2M(body)
expression = {'__','/','*','\^','~','csgn'};
replace = {'_','./','.*','.^','','sign'};
body = regexprep(body, expression, replace);
body = strrep(body,char(10),'');  % remove LF characters
body = body(find(~isspace(body)));
end