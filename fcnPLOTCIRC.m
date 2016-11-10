function [hFig1] = fcnPLOTCIRC(hFig1, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, vecUINF)
%FCNPLOTCIRC Summary of this function goes here
%   Detailed explanation goes here

set(0,'CurrentFigure',hFig1);

hold on

count = 1;
for i = 0:0.05:1
    for j = 1-i:-0.05:0
        lambda(count,1) = i;
        lambda(count,2) = j;
        lambda(count,3) = 1-i-j;
        count = count + 1;
    end
end

for i = 1:valNELE
    
    x1 = matPLEX(1,1,i);
    x2 = matPLEX(2,1,i);
    x3 = matPLEX(3,1,i);
    
    y1 = matPLEX(1,2,i);
    y2 = matPLEX(2,2,i);
    y3 = matPLEX(3,2,i);
    
    z1 = matPLEX(1,3,i);
    z2 = matPLEX(2,3,i);
    z3 = matPLEX(3,3,i);
    
    eta = (x1.*lambda(:,1) + x2.*lambda(:,2) + x3.*lambda(:,3));
    xsi = (y1.*lambda(:,1) + y2.*lambda(:,2) + y3.*lambda(:,3));
    
    circ = matCOEFF(i,1).*(eta.^2) + matCOEFF(i,2).*eta + matCOEFF(i,3).*(xsi.^2) ...
        + matCOEFF(i,4).*xsi + matCOEFF(i,5);
    
    vort = 2.*matCOEFF(i,1).*eta + matCOEFF(i,2) + 2.*matCOEFF(i,3).*xsi + matCOEFF(i,4);
   
    vort2 = matCOEFF(i,1).*eta + matCOEFF(i,2);
    vort3 = matCOEFF(i,3).*xsi + matCOEFF(i,4);
    
    len = length(eta);
    
    % Global
%     etaxsi = fcnTOGLOB(repmat(i,len,1), [eta xsi circ], matDVE, matDVECT, matVLST);
%     DT = delaunay(etaxsi(:,1), etaxsi(:,2));   
%     trisurf(DT, etaxsi(:,1), etaxsi(:,2), etaxsi(:,3) - matCENTER(i,3),'EdgeColor','r','FaceColor','r','FaceAlpha',0.5,'EdgeAlpha',0.5)
    
    etaxsi = fcnTOGLOB(repmat(i,len,1), [eta xsi zeros(len,1)], matDVE, matDVECT, matVLST);
    DT = delaunay(etaxsi(:,1), etaxsi(:,2));   
%     trisurf(DT, etaxsi(:,1), etaxsi(:,2), circ + matCENTER(i,3),'EdgeColor','r','FaceColor','r','FaceAlpha',0.5,'EdgeAlpha',0.5)
    trisurf(DT, etaxsi(:,1), etaxsi(:,2), vort + matCENTER(i,3),'EdgeColor','b','FaceColor','b','FaceAlpha',0.5,'EdgeAlpha',0.5)
%     trisurf(DT, etaxsi(:,1), etaxsi(:,2), vort2,'EdgeColor','g','FaceColor','g','FaceAlpha',0.5,'EdgeAlpha',0.5)
%     trisurf(DT, etaxsi(:,1), etaxsi(:,2), vort3,'EdgeColor','m','FaceColor','m','FaceAlpha',0.5,'EdgeAlpha',0.5)
    
end

end

