clc
clear

xlims = 5;
ylims = 5;
zlims = 5;

x = -xlims:0.1:xlims;
y = -ylims:0.1:ylims;
z = -zlims:0.1:zlims;

[X,Y,Z] = meshgrid(x,y,z);

fpg = unique([reshape(X,[],1) reshape(Y,[],1) reshape(Z,[],1)],'rows');

len = length(fpg(:,1));

dvenum = 