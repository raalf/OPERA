function [hFig1] = fcnGIF(valTIMESTEP, valNELE, matDVE, matVLST, matCENTER, matELST, matDVECT, matPLEX, matCOEFF, matUINF, matROTANG, ...
                    valWNELE, matWDVE, matWVLST, matWCENTER, matWELST, matWDVECT, valWSIZE, valPRESTEPS, matWVGRID, case_num)

hFig1 = fcnPLOTBODY(0, matDVE, valNELE, matVLST, matELST, matDVECT, matCENTER, matPLEX, matCOEFF, matUINF, matROTANG, [], 'opengl');
hFig1 = fcnPLOTWAKE(0, hFig1, matWDVE, valWNELE, matWVLST, matWELST, matWDVECT, matWCENTER, valWSIZE, valPRESTEPS, matWVGRID);
view([-45 20]);

frame = getframe(hFig1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

if ~exist('GIF', 'dir')
   mkdir('GIF');
end

% Write to the GIF File
gif_str = ['GIF/output_',num2str(case_num),'.gif'];

if valTIMESTEP == 1
    imwrite(imind,cm, gif_str,'gif', 'Loopcount',inf);
else
    imwrite(imind,cm, gif_str,'gif','WriteMode','append');
end