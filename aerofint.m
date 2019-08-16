function [xp,yp,xnaca,ynaca] = aerofint(i,ii,n,nn,shg,xnaca,ynaca,xp,yp,lnaca)
    err=1.0e-4;
    fshg0=0.8;
    fshg1=(2-fshg0)/fshg0;
    switch nn
        case 1
            ra0=fshg0*(xp(ii,n+1)-xp(2*ii-i,n+1));
            xp(i,n+1)=xp(ii,n+1)+ra0;
            [useless,xnaca,ynaca,xp,yp]=xylagran(i,n,1,xnaca,ynaca,xp,yp,lnaca);
            yp(i,n+1) = useless;
        case 2
            ra0=(-1)^n*shg*fshg0;
            yp(i,n+1)=yp(ii,n+1)+ra0;
            [useless,xnaca,ynaca,xp,yp]=xylagran(i,n,2,xnaca,ynaca,xp,yp,lnaca);
            xp(i,n+1) = useless;
    end
    ra1=fshg1*ra0;
    ra2=errfn(i,ii,n,shg,xp,yp);
    while (abs(ra2)>err)
        switch nn
            case 1
                xp(i,n+1)=xp(ii,n+1)+ra1;
                [useless,xnaca,ynaca,xp,yp]=xylagran(i,n,1,xnaca,ynaca,xp,yp,lnaca);
                yp(i,n+1) = useless;
            case 2
                yp(i,n+1)=yp(ii,n+1)+ra1;
                [useless,xnaca,ynaca,xp,yp]=xylagran(i,n,2,xnaca,ynaca,xp,yp,lnaca);
                xp(i,n+1) = useless;
        end
        ra3=errfn(i,ii,n,shg,xp,yp);
        res=ra1-ra3*(ra1-ra0)/(ra3-ra2);
        ra0=ra1;
        ra1=res;
        ra2=ra3;
    end

end