function split_RO(this)

%Zerlegung lokal oder global
%Liefert p_new, RO/KR_new

%KR_raw_new = KR_raw;
p = this.Points;


%globale Verfeinerung
if this.RO_factor_global > 1
    
   this.RO_raw(:,4) = this.RO_factor_global;
   
  
   
   RO_raw_new = zeros(this.RO_factor_global * size(this.RO_raw,1),3);
   %RO_raw_new(:,1) = RO_mat;
   p_new = zeros(size(p,1) + size(this.RO_raw,1)*(this.RO_factor_global-1),3);
    
end

%lokaleVerfeinerung
if this.RO_factor_global <= 1
    
   %this.RO_raw(:,4) = this.RO_factor_global;
   
   RO_count_add = 0;
   for j=1:size(this.RO_raw,1)
      if this.RO_raw(j,4) >= 1
      RO_count_add = RO_count_add + this.RO_raw(j,4) - 1; 
      end
   end
   
   %RO_count_add
   
   
   RO_raw_new = zeros(size(this.RO_raw,1)+RO_count_add,3);
   %RO_raw_new(:,1) = RO_mat;
   p_new = zeros(size(p,1) +RO_count_add,3);
    

end



%this.RO_raw


p_new(1:size(p,1),:) = p;
RO_count = 1;
point_count = 1;


for i = 1:size(this.RO_raw,1)

    RO_factor = this.RO_raw(i,4);
    RO_mat_act = this.RO_raw(i,1);
    
    if RO_factor <= 1
        
        %Füge keine neuen Punkte ein, füge altes Rohr wieder ein 
        
                RO_raw_new(RO_count,2) = this.RO_raw(i,2);
                RO_raw_new(RO_count,3) = this.RO_raw(i,3);
                %setze altes Material
                RO_raw_new(RO_count,1) = this.RO_raw(i,1);
                RO_count = RO_count+1;
        
        continue;
    end
    
    %RO_factor
    
    %P_start = this.RO_raw(i,2)
    %P_end = this.RO_raw(i,3)
    
    %x_coord_start = p(this.RO_raw(i,2),1)
    
    lx = p(this.RO_raw(i,3),1) - p(this.RO_raw(i,2),1);
    ly = p(this.RO_raw(i,3),2) - p(this.RO_raw(i,2),2);
    lz = p(this.RO_raw(i,3),3) - p(this.RO_raw(i,2),3);
    dx = lx/RO_factor;
    dy = ly/RO_factor;
    dz = lz/RO_factor;
    
    %Füge neue Punkte ein
    
    for j = 1:RO_factor-1
        
            %size(p,1)
        
          p_new(size(p,1)+point_count,1) =  p(this.RO_raw(i,2),1) + j*dx   ;  
          p_new(size(p,1)+point_count,2) =  p(this.RO_raw(i,2),2) + j*dy   ;  
          p_new(size(p,1)+point_count,3) =  p(this.RO_raw(i,2),3) + j*dz   ;
         
          
          
          
          
          P_start = this.RO_raw(i,2);
          
          
          P_act = size(p,1)+point_count;
          
          %P_act
          
          point_count = point_count+1;
          
          P_end = this.RO_raw(i,3);
          
          %Füge Anfangsrohr ein
          if (j==1)
                RO_raw_new(RO_count,2) = P_start;
                RO_raw_new(RO_count,3) = P_act;
                RO_raw_new(RO_count,1) = RO_mat_act;
                RO_count = RO_count+1;
             
          end
          
          %Füge Zwischenrohr(e) ein
          %j
          if(j~=1 )
                RO_raw_new(RO_count,2) = P_old;
                RO_raw_new(RO_count,3) = P_act;
                RO_raw_new(RO_count,1) = RO_mat_act;
                RO_count = RO_count+1;
                %RO_count
          end
          
          %Füge Endrohr ein
          if (j==RO_factor-1)
                RO_raw_new(RO_count,2) = P_act;
                RO_raw_new(RO_count,3) = P_end;
                RO_raw_new(RO_count,1) = RO_mat_act;
                RO_count = RO_count+1;
          end
          
          P_old = P_act;
        
    end
    
    
   
%     %Füge neue Rohre ein
%     for j = 1:RO_factor
%        
%         %Startpunkt
%         RO_raw_new(j+(i-1)*(RO_factor),2) = size(p,1)+j+(i-1)*(RO_factor-1);
%         %Endpunkt
%         RO_raw_new(j+(i-1)*(RO_factor),3) = size(p,1)+j+1+(i-1)*(RO_factor-1);
%         
%     end
    
    
    
end

%if this.RO_factor_global ~= 1
%RO_raw_new

% fprintf('\nSize RO_new: %d\n',size(RO_raw_new,1));
% fprintf('Size P_new: %d\n',size(p_new,1));
% 
% fprintf('\nSize RO: %d\n',size(this.RO_raw,1));
% fprintf('Size P: %d\n',size(p,1));

%this.RO_raw = RO_raw_new;
%p=p_new;

%if KR_factor_global ~= 1
%KR_raw_new

fprintf('\nNeue Größen\n');
fprintf('Size RO_new: %d\n',size(RO_raw_new,1));
fprintf('Size P_new: %d\n',size(p_new,1));

fprintf('\nAlte Größen\n');
fprintf('Size RO: %d\n',size(this.RO_raw,1));
fprintf('Size P: %d\n',size(p,1));

this.Points = p_new;
this.RO_raw = RO_raw_new;
% p_new2
% KR_raw_new

% fprintf('\nSize RO_new: %d\n',size(RO_raw_new,1));
% fprintf('Size P_new: %d\n',size(p_new,1));
% 
% fprintf('\nSize RO: %d\n',size(this.RO_raw,1));
% fprintf('Size P: %d\n',size(p,1));

