function [points, RO, KR, FH, mat, lager, lasten] = read_file(~, file)

%-----------------------------------------------------------
%Read Geometry Data from file
%Output: points,Rohre,mat,lager,lasten
%Output Rohre: RO,KR, (?)
%-----------------------------------------------------------

% clc;

% %-----------------------------------------------------------
% %Read Materialparameter
% %-----------------------------------------------------------
% mat = read_in3('MAT',4,file);
%
% %-----------------------------------------------------------
% %Read Points
% %-----------------------------------------------------------
% points = read_in3('Knoten','Elemente',4,file);
%
%
% %-----------------------------------------------------------
% %Read Rohre
% %-----------------------------------------------------------
% %-----------------------------------------------------------
% %Read RO Rohre
% %-----------------------------------------------------------
% RO = read_in('Elemente','Lager','RO',3,file);
% %-----------------------------------------------------------
% %Read KR Rohre
% %-----------------------------------------------------------
% KR = read_in('Elemente','Lager','KR',4,file);
%
%
% %-----------------------------------------------------------
% %Read DIR Lager
% %-----------------------------------------------------------
% lager = read_in('Lager','Lasten','DIR',13,file);
%
%
% %-----------------------------------------------------------
% %Read NEU Lasten
% %-----------------------------------------------------------
% lasten = read_in('Lasten','Ende','NEU',7,file);

%-----------------------------------------------------------
%Read Materialparameter
%-----------------------------------------------------------
mat = read_in3('MAT',12,file);
mat = mat(:,2:12);

%-----------------------------------------------------------
%Read Points
%-----------------------------------------------------------
points = read_in3('KNO',4,file);
points = points(:,2:4);

%-----------------------------------------------------------
%Read RO Rohre
%-----------------------------------------------------------
RO = read_in3('RO',4,file);
%-----------------------------------------------------------
%Read KR Rohre
%-----------------------------------------------------------
KR = read_in3('KR',5,file);
%-----------------------------------------------------------
%Read FH Federhänger
%-----------------------------------------------------------
FH = read_in3('FH',3,file);

%-----------------------------------------------------------
%Read DIR Lager
%-----------------------------------------------------------
lager = read_in3('DIR',15,file);


%-----------------------------------------------------------
%Read NEU Lasten
%-----------------------------------------------------------
lasten = read_in3('NEU',8,file);

    function output = read_in3(keyword,count,file)
        
        %-----------------------------------------------------------
        %Read_in mit Keyword über gesamtes File
        %-----------------------------------------------------------
        
        fid = fopen(file);
        k = 1;
        
        while ~feof(fid)
            
            
            str = fscanf(fid, '%s', 1);
            
            
            if (strcmp(str,'#'))
                str_temp = fscanf(fid, '%s', 1);
                continue;
            end
            
            
            if (strcmp(str,keyword))
                
                for i=1:count
                    val = fscanf(fid, '%g', 1);
                    
                    if isempty(val)
                        
                        output(k) = 0;
                        k=k+1;
                        
                        continue;
                    end
                    
                    output(k) = val;
                    k=k+1;
                    
                end
            end
            
        end
        
        
        %output
        
        fclose(fid);
        
        if exist('output')
            
            fprintf('Read_in: \t Anzahl: %g \t Keyword: %s \n',length(output)/count,keyword);
            output =  reshape(output,count,length(output)/count)';
            
        else
            output = [];
            fprintf('Keyword %s does not exist\n',keyword);
        end
        
        
    end
end

