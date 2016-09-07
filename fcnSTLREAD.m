function [temp, isBinary] = fcnSTLREAD(stlfile)

%{
% fgetl(fp); % Skipping header line

% count = 1;
% while ~feof(fp)
%     try
%         % Skip "facet normal" words
%         fscanf(fp,'%s',2);
%
%         % Normal
%         temp(count,:,4) = fscanf(fp,'%f',3);
%         fgetl(fp); % Skip to new line
%
%         % Skip "outer loop" line
%         fgetl(fp);
%
%         % Skip "vertex" words
%         fscanf(fp,'%s',1);
%
%         % Vertex 1 points
%         temp(count,:,1) = fscanf(fp,'%f',3);
%         fgetl(fp); % Skip to new line
%
%         % Skip "vertex" words
%         fscanf(fp,'%s',1);
%
%         % Vertex 2 points
%         temp(count,:,2) = fscanf(fp,'%f',3);
%         fgetl(fp); % Skip to new line
%
%         % Skip "vertex" words
%         fscanf(fp,'%s',1);
%
%         % Vertex 3 points
%         temp(count,:,3) = fscanf(fp,'%f',3);
%         fgetl(fp); % Skip to new line
%
%         % Skip "endloop" line
%         fgetl(fp);
%
%         % Skip "endfacet" line
%         fgetl(fp);
%
%         count = count + 1;
%     catch
%         break;
%     end
% end
% fclose(fp);
%}

if nargin == 0
    warning('No STL-file is specified');
end

try
    % Try to read an STL ASCII file
    temp = stlAread(stlfile);
    isBinary=false;
catch
    
    try
        % Try to read an STL binary file
        temp = stlBread(stlfile);
        isBinary=true;
        
    catch
        error('File could not be read!')
    end
    
end


function temp = stlAread(stlfile)
% Reads an STL ASCII file

fid=fopen(stlfile,'r');

fileTitle=sscanf(fgetl(fid),'%*s %s');

vnum=0;
fclr=0;
testASCII=true;
lineCount=0;

while feof(fid) == 0
    
    stlLine=fgetl(fid);
    keyWord=sscanf(stlLine,'%s');
    
    if strncmpi(keyWord,'c',1) == 1;
        fclr = sscanf(stlLine,'%*s %f %f %f');
    elseif strncmpi(keyWord,'v',1) == 1;
        vnum = vnum+1;
        vertex(:,vnum) = sscanf(stlLine,'%*s %f %f %f');
        clr(:,vnum) = fclr;
    elseif testASCII
        lineCount=lineCount+1;
        
        if lineCount>20
            if vnum>2
                testASCII=false;
            else
                error('File is not an STL ASCII file!')
            end
        end
    end
end

temp(:,1:3,1) = [vertex(1,1:3:end)' vertex(2,1:3:end)' vertex(3,1:3:end)'];
temp(:,1:3,2) = [vertex(1,2:3:end)' vertex(2,2:3:end)' vertex(3,2:3:end)'];
temp(:,1:3,3) = [vertex(1,3:3:end)' vertex(2,3:3:end)' vertex(3,3:3:end)'];

fclose(fid);


function temp = stlBread(stlfile)
% Reads an STL binary file

fid = fopen(stlfile,'r');

fileTitle = fread(fid,80,'uchar=>schar');
fnum = fread(fid,1,'int32');

temp = zeros(fnum,3,3);

FVCD = uint8(zeros(3,fnum));

for i=1:fnum
    
%     if mod(fnum/100, i) == 0
%        perc = 100*(i/fnum);
%        disp(perc);
%     end
    
    normal = fread(fid,3,'float32');
    vertex1 = fread(fid,3,'float32');
    vertex2 = fread(fid,3,'float32');
    vertex3 = fread(fid,3,'float32');
    clr = fread(fid,1,'uint16');
    
    if bitget(clr,16) == 1
        rd = bitshift(bitand(65535,clr),-10);
        grn = bitshift(bitand(2047,clr),-5);
        bl = bitand(63,clr);
        FVCD(:,i) = [rd;grn;bl];
    end
    
    temp(i,1:3,1) = [vertex1(1) vertex1(2) vertex1(3)];
    temp(i,1:3,2) = [vertex2(1) vertex2(2) vertex2(3)];
    temp(i,1:3,3) = [vertex3(1) vertex3(2) vertex3(3)];
    
end

fclose(fid);
