 if wi
                    iap = x(4*m+1:5*m);
                    bar = x(5*m+1:6*m);
                    yai = x(6*m+1:7*m);
                    xab = x(7*m+1:end);
 end
                
   % iap
                fx(4*m+1:5*m) = -(rc(3)+rc(4))*ya.*iap - rc(8)*iap + rc(15) + rc(14)*yai;
                % bar
                fx(5*m+1:6*m) = -rc(11)*xa.*bar + rc(18)*xab - rc(19)*bar + rc(19);
                % yai
                fx(6*m+1:7*m) = rc(3)*ya.*iap -(rc(14)+rc(7))*yai;
                % xab
                fx(7*m+1:8*m) = rc(11)*xa.*bar-(rc(18)+rc(13))*xab;