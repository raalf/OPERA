clc
clear

% letter = {'A','B','C','D','E','F','G','H'}
letter = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T'};

fn = 0:7;

for i = 1:length(fn)

fcn = ['K',num2str(fn(i))]; 

fcnheader = sprintf('function I = fcn%s(S, t, u, alpha, F, tol, idx)',fcn);

% fcnheader = sprintf('function I = fcn%s(S, t, u, alpha, F, tol,',fcn);
% 
% for jj = 1:length(letter)
% 
%     if jj < length(letter)
%     fcnheader = [fcnheader, ' idx_', letter{jj}, ','];
%     else
%     fcnheader = [fcnheader, ' idx_', letter{jj}, ')'];  
%     end
%     
% end

bdy1 = sprintf('\nlen = size(S,1);\nI = nan(len,2);\n');

bdy2 = '';
for jj = 1:length(letter)
%     bdy2 = [bdy2; sprintf('\nI(idx_%s) = fcnK%d%s(S(idx_%s), t(idx_%s), u(idx_%s), alpha(idx_%s), reshape(F(idx_%s(:,:,[1,1])),[],1,2), tol);', letter{jj}, fn(i), letter{jj}, letter{jj}, letter{jj}, letter{jj}, letter{jj}, letter{jj})];
    bdy2{jj} = sprintf('\nidx_tmp = logical(idx(:,:,%d));\nI(idx_tmp) = fcnK%d%s(S(idx_tmp), t(idx_tmp), u(idx_tmp), alpha(idx_tmp), reshape(F(idx_tmp(:,:,[1,1])),[],1,2), tol);', jj, fn(i), letter{jj});
end
% bdy2 = convertCharsToStrings(bdy2');
footer = sprintf('\n\nend\n');

whole_thing = strjoin(convertCharsToStrings([fcnheader bdy1 bdy2 footer]))

fid = fopen(['G:\GIT\opera\K Functions\fcn', fcn, '.m'],'wt');
% fid = fopen(['C:\Users\travi\OneDrive\Desktop\GIT\opera\fcn', fcn, '.m'],'wt');

fprintf(fid, whole_thing);
fclose(fid);

end








