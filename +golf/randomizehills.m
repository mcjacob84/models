[h,r] = generatehills(size(hills,1),max_hillheight,max_hillradius,min_hillradius);
hills(:,3) = h; hills(:,4) = r;
% Lochdaten anpassen
HZ = green(HX,HY);
plotgreen();
% Erzeugt neue höhen für das Gras