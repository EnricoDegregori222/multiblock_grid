%% Input
filename_input =input('Please enter the input filename:', 's');
% filename_input = strcat('input/', filename_input);
run(strcat('input/', filename_input));

%% Airfoil section grid
xp = zeros(lxit+1,4);
yp = xp;
xq = zeros(lett+1,4);
yq = xq;
pxi = zeros(lxit+1,1);
qet = zeros(lett+1,1);
qeto = zeros(lett+1,2);
aaa = zeros(length(yq),1);
bbb = zeros(length(xp),1);
ccc = pxi;
ddd = qet;

xy = dlmread(strcat('aerofoils/', filename_input));
lnaca = length(xy)/2;
xnaca(:,1) = xy(1:lnaca,1);
xnaca(:,2) = xy(lnaca+1:end,1);
ynaca(:,1) = xy(1:lnaca,2);
ynaca(:,2) = xy(lnaca+1:end,2);
figure(1)
plot(xnaca(:,1),ynaca(:,1), '-ok'); hold on; grid on;
plot(xnaca(:,2),ynaca(:,2), '-ok');

shs=smg;
she=shs;

is=1;
ie=lnaca;
 for m=0:nref
     for i=is
        res=0.5*(xnaca(i,1)+xnaca(i,2));
        xnaca(i,1)=res;
        xnaca(i,2)=res;
        res=0.5*(ynaca(i,1)+ynaca(i,2));
        ynaca(i,1)=res;
        ynaca(i,2)=res;
     end
        cha(:)=[-1,4,10,4,-1]/16.0;
        xf(:,:)=xnaca(:,:);
        yf(:,:)=ynaca(:,:);
     for n=1:2
     for i=is+2:ie-2
        xf(i,n)=sum(cha(:).*xnaca(i-2:i+2,n));
        yf(i,n)=sum(cha(:).*ynaca(i-2:i+2,n));
     end
        xf(is,n)=sum(cha(:).*[xnaca(is+2:-1:is+1,3-n);xnaca(is:is+2,n)]);
        yf(is,n)=sum(cha(:).*[ynaca(is+2:-1:is+1,3-n);ynaca(is:is+2,n)]);
        xf(is+1,n)=sum(cha(:).*[xnaca(is+1,3-n);xnaca(is:is+3,n)]);
        yf(is+1,n)=sum(cha(:).*[ynaca(is+1,3-n);ynaca(is:is+3,n)]);
     end
     for n=1:2
         for i=is:ie-2
            xnaca(i,n)=xf(i,n);
            ynaca(i,n)=yf(i,n);
         end
     end
 end
 
 xb = 0.5*(xnaca(is,1)+xnaca(is,2));
 yb = 0.5*(ynaca(is,1)+ynaca(is,2));
 xc(1) = xnaca(ie,1);
 yc(1) = ynaca(ie,1);
 xc(2) = xnaca(ie,2);
 yc(2) = ynaca(ie,2);
 
figure(1)
plot(xnaca(:,1),ynaca(:,1), '-or'); hold on; grid on;
plot(xnaca(:,2),ynaca(:,2), '-or');
title('Compare original and smoothed airfoil');
figure(2)
plot(xnaca(:,1),ynaca(:,1), '-or'); hold on; grid on;
plot(xnaca(:,2),ynaca(:,2), '-or');


 for n=1:2
    xp(lxis1+1,n+1)=xb;
    yp(lxis1+1,n+1)=yb;
    xp(lxie1+1,n+1)=xc(n);
    yp(lxie1+1,n+1)=yc(n);
    p=polyfit(xnaca(end-1:end,n),ynaca(end-1:end,n),1);
    xp(lxie1,n+1) = xnaca(end,n)-(smg*1.);
    yp(lxie1,n+1) = polyval(p,xnaca(end,n)-smg);
    ll=n_initial_points; % "ll" must be equal to or larger than 4.
    for i=lxis1+1+1:lxis1+ll+1
         if (i<=lxis1+2+1)
            [xp,yp,xnaca,ynaca] = aerofint(i,i-1,n,2,shs,xnaca,ynaca,xp,yp,lnaca);
        else
            [xp,yp,xnaca,ynaca] = aerofint(i,i-1,n,1,shs,xnaca,ynaca,xp,yp,lnaca);
        end
    end
    for i=lxie1-2+1:-1:lxie1-ll+1
        if(i>=lxie1-1+1)
            [xp,yp,xnaca,ynaca] = aerofint(i,i+1,n,2,she,xnaca,ynaca,xp,yp,lnaca);
        else
            [xp,yp,xnaca,ynaca] = aerofint(i,i+1,n,1,she,xnaca,ynaca,xp,yp,lnaca);
        end
    end
    tmpa=sum(xp(lxis1+ll-4+1:lxis1+ll+1,n+1).*[3.0,-16.0,36.0,-48.0,25.0]')/12.0;
    tmpb=sum(xp(lxie1-ll+4+1:-1:lxie1-ll+1,n+1).*[3.0,-16.0,36.0,-48.0,25.0]')/12.0;
    ip=lxis1+ll+1;
    im=lxi1-2*ll;
    [bbb,ccc] = gridf(xp(:,n+1),pxi,xp(lxis1+ll+1,n+1),xp(lxie1-ll+1,n+1),tmpa,-tmpb,lxit,im,ip);
    xp(:,n+1) = bbb;
    pxi=ccc;
    for i=lxis1+ll+1+1:lxie1-ll-1+1
        [useless,xnaca,ynaca,xp,yp] = xylagran(i,n,1,xnaca,ynaca,xp,yp,lnaca);
        yp(i,n+1)= useless;
    end
 end

figure(2)
plot(xp(lxis1+1:lxie1+1,2),yp(lxis1+1:lxie1+1,2), '-o'); hold on; grid on;
plot(xp(lxis1+1:lxie1+1,3),yp(lxis1+1:lxie1+1,3), '-o');
title('Compare smoothed and meshed airfoil');