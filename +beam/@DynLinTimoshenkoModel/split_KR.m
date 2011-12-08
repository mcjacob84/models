function split_KR(this)

%Zerlegung lokal oder global
%Liefert p_new, RO/KR_new

p = this.Points;

if (size(this.KR_raw,1) == 0)
    return;
end

%Zerlege KR Rohre

%globale Verfeinerung
if this.KR_factor_global > 1
    
    this.KR_raw(:,5) = this.KR_factor_global;
    
    KR_raw_new = zeros(this.KR_factor_global * size(this.KR_raw,1),5);
    %KR_raw_new(:,1) = KR_mat;
    p_new = zeros(size(p,1) + size(this.KR_raw,1)*(this.KR_factor_global-1),3);
    
    KR_raw_new(:,5) = this.KR_factor_global;
    
end

%lokaleVerfeinerung
if this.KR_factor_global <= 1
    
    %this.KR_raw(:,4) = this.KR_factor_global;
    
    KR_count_add = 0;
    for j=1:size(this.KR_raw,1)
        if this.KR_raw(j,5) >= 1
            KR_count_add = KR_count_add + this.KR_raw(j,5) - 1;
        end
    end
    
    KR_raw_new = zeros(size(this.KR_raw,1)+KR_count_add,5);
    %KR_raw_new(:,1) = KR_mat;
    p_new = zeros(size(p,1) +KR_count_add,3);
end

p_new(1:size(p,1),:) = p;
KR_count = 1;
point_count = 1;
for i = 1:size(this.KR_raw,1)
    
    KR_factor = this.KR_raw(i,5);
    KR_mat_act = this.KR_raw(i,1);
    
    if KR_factor <= 1
        
        %Füge keine neuen Punkte ein, füge altes Rohr wieder ein
        
        KR_raw_new(KR_count,2) = this.KR_raw(i,2);
        KR_raw_new(KR_count,3) = this.KR_raw(i,3);
        KR_raw_new(KR_count,4) = this.KR_raw(i,4);
        KR_raw_new(KR_count,5) = 1;
        %setze altes Material
        KR_raw_new(KR_count,1) = this.KR_raw(i,1);
        KR_count = KR_count+1;
        
        continue;
    end

    %Punktkoordinaten der Zwischenpunkte
    
    % Anfangs- und Endpunkt
    P_start = [this.KR_raw(i,2) this.KR_raw(i,3)];
    P_end = this.KR_raw(i, 4);

    KR_factor_act = this.KR_raw(i,5);
    s = ( 0 : (pi/2)/KR_factor_act : pi/2 );
    
    e_x = ( p_new(P_start(1),:) - p_new(P_end,:) )';
    rad = norm(e_x);
    e_x = e_x / norm(e_x);
    e_y = (p(P_start(2),:) - p(P_end,:))' /norm(e_x);
    e_y = e_y - (e_x'*e_y) * e_x;
    e_y = e_y / norm(e_y);
    e_z = [e_x(2)*e_y(3) - e_x(3)*e_y(2);
        e_x(3)*e_y(1) - e_x(1)*e_y(3);
        e_x(1)*e_y(2) - e_x(2)*e_y(1)];
    e_z = e_z / norm(e_z);
    T = [e_x e_y e_z];

    x = rad * cos(s);
    y = rad * sin(s);
    z = 0*s;
    COR = T * [x; y; z];
    COR2 = [( p_new(P_end, 1) + COR(1,:))' ( p_new(P_end, 2) + COR(2,:))' ( p_new(P_end, 3) + COR(3,:))'];
    
    for j = 1:KR_factor-1
        
        p_new(size(p,1)+point_count,1) =  COR2(j+1,1);
        p_new(size(p,1)+point_count,2) =  COR2(j+1,2);
        p_new(size(p,1)+point_count,3) =  COR2(j+1,3);
        
        P_act = size(p,1)+point_count;
        
        point_count = point_count+1;
        
        %Füge Anfangsrohr ein
        if (j==1)
            KR_raw_new(KR_count,2) = P_start(1);
            KR_raw_new(KR_count,3) = P_act;
            
            KR_raw_new(KR_count,4) = P_end;
            
            KR_raw_new(KR_count,5) = KR_factor_act;
            
            KR_raw_new(KR_count,1) = KR_mat_act;
            KR_count = KR_count+1;
            
        end
        
        %Füge Zwischenrohr(e) ein
        %j
        if(j~=1 )
            KR_raw_new(KR_count,2) = P_old;
            KR_raw_new(KR_count,3) = P_act;
            
            KR_raw_new(KR_count,4) = P_end;
            
            KR_raw_new(KR_count,5) = KR_factor_act;
            
            KR_raw_new(KR_count,1) = KR_mat_act;
            KR_count = KR_count+1;
            %KR_count
        end
        
        %Füge Endrohr ein
        if (j==KR_factor-1)
            KR_raw_new(KR_count,2) = P_act;
            KR_raw_new(KR_count,3) = P_start(2);
            
            KR_raw_new(KR_count,4) = P_end;
            
            KR_raw_new(KR_count,5) = KR_factor_act;
            
            KR_raw_new(KR_count,1) = KR_mat_act;
            KR_count = KR_count+1;
        end
        
        P_old = P_act;
    end
end

fprintf('\nNeue Größen\n');
fprintf('Size KR_new: %d\n',size(KR_raw_new,1));
fprintf('Size P_new: %d\n',size(p_new,1));


fprintf('\nAlte Größen\n');
fprintf('Size KR: %d\n',size(this.KR_raw,1));
fprintf('Size P: %d\n',size(p,1));

this.Points = p_new;
this.KR_raw = KR_raw_new;
