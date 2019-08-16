function [xy,xnaca,ynaca,xp,yp] = xylagran(i,n,nn,xnaca,ynaca,xp,yp,lnaca)
    
    is=1;
    ie=lnaca;
    switch nn
        
        case 1
        [AAA,ii]=min((xnaca(is:ie-1,n)-xp(i,n+1)).*(xnaca(is+1:ie,n)-xp(i,n+1)));
        ii=ii-1;
        maxjj=2;
        ip=min(max(ii,is+1-1),ie-maxjj-1);
        xy=0;
        alag(:)=xp(i,n+1)-xnaca(ip-1+1:ip+maxjj+1,n);
        for jj=-1:maxjj
            blag(:)=xnaca(ip+jj+1,n)-xnaca(ip-1+1:ip+maxjj+1,n);
            ao=1;
            bo=1;
            for ii=-1:maxjj
                if(ii~=jj)
                    ao=ao*alag(ii+2);
                    bo=bo*blag(ii+2);
                end
            end
            xy=xy+ao*ynaca(ip+jj+1,n)/bo;
        end
        
        case 2
        ii=3;
        ao=abs(xp(i,n+1)-xnaca(is,n));
        bo=abs(xp(i,n+1)-xnaca(ie,n));
        if(ao<bo)
            [AAA,ip]=min((ynaca(is:is+ii-1,n)-yp(i,n+1)).*...
                      (ynaca(is+1:is+ii,n)-yp(i,n+1)));
            ip=ip-1;
            
        else
            [AAA,ip]=min((ynaca(ie-ii+1:ie,n)-yp(i,n+1)).*...
                      (ynaca(ie-ii:ie-1,n)-yp(i,n+1)));
                  ip=ip+ie-ii-1;
        end
        
        xlag(is+1:ie+1)=xnaca(:,n);
        xlag(is)=xnaca(is+1,3-n);
        xlag([ie+2,ie+3])=xnaca([ie-1,ie-2],3-n);
        ylag(is+1:ie+1)=ynaca(:,n);
        ylag(is)=ynaca(is+1,3-n);
        ylag([ie+2,ie+3])=ynaca([ie-1,ie-2],3-n);
        
        xy=0;
        ipn=2;
        alag(:)=yp(i,n+1)-ylag(ip-1+ipn:ip+2+ipn);
        for jj=-1:2
            blag(:)=ylag(ip+jj+ipn)-ylag(ip-1+ipn:ip+2+ipn);
            ao=1;
            bo=1;
            for ii=-1:2
                if(ii~=jj)
                    ao=ao*alag(ii+2);
                    bo=bo*blag(ii+2);
                end
            end
            xy=xy+ao*xlag(ip+jj+ipn)/bo;
        end
        
    end

end