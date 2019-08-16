%% Aerofoil section
grid_aero

if gridsave==1
    filename_output = input('Please enter the grid output filename:', 's');
end

lximb = [lxi0,lxi1,lxi2,lxi0,lxi1,lxi2,lxi2];
letmb = [let0,let0,let0,let0,let0,let0,let1];
lzemb = [lze0,lze0,lze0,lze0,lze0,lze0,lze0];
if gridtec==1
    for zz=1:mbk
        xxxx{zz} = zeros(lximb(zz)+1, letmb(zz)+1,lzemb(zz)+1);
        yyyy{zz} = zeros(lximb(zz)+1, letmb(zz)+1,lzemb(zz)+1);
    end
end
%% Spanwise extension
xcLE = xc;
ycLE = yc;
xpo(lxis1+1:lxie1+1,2:3) = xp(lxis1+1:lxie1+1,2:3);
ypo(lxis1+1:lxie1+1,2:3) = yp(lxis1+1:lxie1+1,2:3);
kplot = 1;
for k=0+1:lze0+1
    zs(k) = span*((lze0+1-k) / (lze0) - 0.5);
    xa=-doml0;
    xd=doml1-szth1;
    xe=doml1;
    ya=-domh;
    yd=domh;
    fctr=2*pi/wlew;
    shswle=shs*sqrt(1+(fctr*wlea*cos(fctr*(zs(k)-zs(0+1))))^2);

%% Wavy leading edge and bump
tmp=1-wlea*sin(2*pi*(zs(k)-zs(0+1))/wlew)
xcm = mean(xcLE);
ycm = mean(ycLE);
for n=1:2
    for i=lxis1+1:lxie1+1
        xp(i,n+1)=tmp*(xpo(i,n+1)-xcm)+xcm;
        yp(i,n+1)=tmp*(ypo(i,n+1)-ycm)+ycm;
    end
    xc(n) = xp(lxie1+1,n+1);
    yc(n) = yp(lxie1+1,n+1);
end

xb=0.5*(xp(lxis1+1,1+1)+xp(lxis1+1,2+1)); 
yb=0.5*(yp(lxis1+1,1+1)+yp(lxis1+1,2+1));

%% Horizontal interface
sho=szth0 / nsponge;

for n=1:2
ip=0+1;
im=lxi0-lxi01;
[bbb,ccc] = gridf(xp(:,n+1),pxi,xa,xb-doml01,sho,free,lxit,im,ip);
xp(:,n+1) = bbb;
pxi = ccc;

ip = lxi0-lxi01+1;
im = lxi01;
[bbb,ccc] = gridf(xp(:,n+1),pxi,xb-doml01,xb,pxi(ip-1),shswle,lxit,im,ip);
xp(:,n+1) = bbb;
pxi = ccc;

if k==0+1
    [AAA,lh]=min(abs(xa+szth0-xp(1:lxi0,n+1)));
    lh=lh-1;
    lxisz=lh;
end

ip=lxis2+1;
im=lxi2-lxisz;
[bbb,ccc] = gridf(xp(:,n+1),pxi,xc(n),xd,she,free,lxit,im,ip);
xp(:,n+1) = bbb;
pxi = ccc;

ip=ip+im;
im=lxisz;
[bbb,ccc] = gridf(xp(:,n+1),pxi,xd,xe,pxi(ip-1),free,lxit,im,ip);
xp(:,n+1) = bbb;
pxi=ccc;

tmpa=pxi(lxisz+1);
tmpb=pxi(lxit-lxisz-1);
for m=1:2
    switch m
        case 1
            is=1;
            ie=lxie0+1;
            ii=is;
            xo=xb-xa;
            yo=yb;
            ygo=0;
        case 2
            is=lxis2+1;
            ie=lxit+1;
            ii=ie;
            xo=xe-xc(n);
            yo=yc(n);
            ygo=-yc(n)+(2-n)*spy;
    end
    am=yo/sin(pi/2)^2;
    for i=is:ie
        gf = (-1)^(m+1)*(xp(i,n+1)-xp(ii,n+1))/xo;
        yp(i,n+1) = am*sin(pi/2*gf)^2 - ygo*cos(pi/2*gf);
        yp(i,n+1) = yo;
    end
end
end

if k==1
    figure(3)
    plot(xp(:,2),yp(:,2), 'k', 'LineWidth', 2); grid on; hold on;
    plot(xp(:,3),yp(:,3), 'k', 'LineWidth', 2);
end

%% Top and bottom boundaries
fctr=3;
sho=szth0 / nsponge;

for n=0:3:3
    switch n
        case 0
            nn=1;
            yv=ya;
        case 3
            nn=2;
            skew=-skew;
            yv=yd;
    end
ip=0+1;
im=lxi0-lxi01;
[bbb, ccc] =  gridf(xp(:,n+1),pxi,xa-skew,xb-skew-spx-doml01,sho,free,lxit,im,ip);
xp(:,n+1) = bbb;
pxi = ccc;

ip = lxi0-lxi01+1;
im = lxi01;
[bbb,ccc] = gridf(xp(:,n+1),pxi,xb-skew-spx-doml01,xb-skew-spx,pxi(ip-1),fctr*shs,lxit,im,ip);
xp(:,n+1) = bbb;
pxi = ccc;

ip=lxis1+1;
im=lxi1;
[bbb,ccc] =  gridf(xp(:,n+1),pxi,xb-skew-spx,xc(nn)-skew+0*spx,fctr*shs,1/fctr*fctr*she,lxit,im,ip);
xp(:,n+1) = bbb;
pxi = ccc;

ip=lxis2+1;
im=lxi2-lxisz;
[bbb,ccc] =  gridf(xp(:,n+1),pxi,xc(nn)-skew+0*spx,xd-skew,1/fctr*fctr*she,free,lxit,im,ip);
xp(:,n+1) = bbb;
pxi = ccc;

ip=ip+im;
im=lxisz;
[bbb,ccc] =  gridf(xp(:,n+1),pxi,xd-skew,xe-skew,pxi(ip-1),free,lxit,im,ip);
xp(:,n+1) = bbb;
pxi = ccc;

skew=-skew;

for i=1:lxit+1
    yp(i,n+1)=yv;
end
end

if k==1
    figure(3)
    plot(xp(:,1), yp(:,1), 'k', 'LineWidth', 2);
    plot(xp(:,4), yp(:,4), 'k', 'LineWidth', 2);
end

%% Vertical interface
fctr=smgvr*sqrt(0.5);
ya = -domh1;
yd = domh1;
yq1 = yq;
for m=1:2
    switch m
        case 1
            yo1=yb;
            yo2=yb;
            sho=fctr*shs;
        case 2
            yo1=yc(2);
            yo2=yc(1);
            sho=fctr*she*sqrt(2);
    end
    jj=ptvr*let01;
    ra0=fracvr*domh1;
    
    jp=let01-jj+1;
    jm=jj;
    [aaa,ddd] = gridf(yq(:,m+1),qet,yo2-ra0,yo2,free,sho,lett,jm,jp);
    yq(:,m+1) = aaa;
    qet = ddd;
    
    jp=0+1;
    jm=let01-jj;
    [aaa,ddd] = gridf(yq(:,m+1),qet,ya,yo2-ra0,free,qet(jp+jm),lett,jm,jp);
    yq(:,m+1) = aaa;
    qet = ddd;
    
    jp=lets11+1;
    jm=jj;
    [aaa,ddd] = gridf(yq(:,m+1),qet,yo1,yo1+ra0,sho,free,lett,jm,jp);
    yq(:,m+1) = aaa;
    qet = ddd;
    
    jp=jp+jm;
    jm=let01-jj;
    [aaa,ddd] = gridf(yq(:,m+1),qet,yo1+ra0,yd,qet(jp-1),free,lett,jm,jp);
    yq(:,m+1) = aaa;
    qet = ddd;
    
    jp = 0+1;
    jm = let0-let01;
    [aaa,ddd] = gridf(yq1(:,m+1),qet,-domh,ya,free,qet(1),lett,jm,jp);
    yq1(:,m+1) = aaa;
    qet1 = ddd;
    
    jp = lett+1 - (let0-let01);
    jm = let0-let01;
    [aaa,ddd] = gridf(yq1(:,m+1),qet,yd,domh,qet(lets11+1+jj+let01-jj),free,lett,jm,jp);
    yq1(:,m+1) = aaa;
    qet1 = ddd;
    
    qeto(:,m)=qet(:);
end

for n=0:2:2
    switch n
        case 0
            js=1;
            je=lete01+1;
        case 2 
            js=lets11+1;
            je=lett1+1;
    end
    for m=1:2
        res=skew/domh;
        switch m
            case 1
                i=lxis1+1;
                pp=-fix(n/2);
                qq=(1-fix(n/2));
            case 2
                i=lxie1+1;
                pp=fix(n/2);
                qq=(fix(n/2)-1);
        end
        tmpa=xp(i,n+1+1)-xp(i,n+1);
        if n==0
            tmpb = yq(let01+1,m+1) - yq(1,m+1);
        else
            tmpb = -(yq(lets11+1,m+1) - yq(lets11+1+let01,m+1));
        end
        ra1=pp;
        ra2=(3*tmpa-(2*pp+qq)*tmpb)/tmpb^2;
        ra3=((pp+qq)*tmpb-2*tmpa)/tmpb^3;
        
        for j=js:je
            tmp=yq(j,m+1)-yp(i,n+1);
            if n == 0
                tmp = yq(j,m+1)-yq(1,m+1);
            end
            xq(j,m+1)=xp(i,n+1)+ra1*tmp+ra2*tmp^2+ra3*tmp^3;
        end
        if n==0
            jk=js;
        else
            jk=je;
        end
        xq1(m+1) = xq(jk,m+1);
    end
end

xq(:,3) = xc(1)*ones(length(yq(:,3)),1);

yq2 = yq;
xq2 = xq;
for m=1:2
    yq(1:let0-let01,m+1) = yq1(1:let0-let01,m+1);
    yq(let0-let01+1:let0+1,m+1) = yq2(1:let01+1,m+1);
    yq(let0+2:let0+2+let01,m+1) = yq2(lets11+1:lets11+1+let01,m+1);
    yq(lett+1-(let0-let01)+1:lett+1,m+1) = yq1(lett+1 - (let0-let01)+1:lett+1,m+1);
    
    xq(1:let0-let01,m+1) = xq1(m+1)*ones(let0-let01,1);
    xq(let0-let01+1:let0+1,m+1) = xq2(1:let01+1,m+1);
    xq(let0+2:let0+2+let01,m+1) = xq2(lets11+1:lets11+1+let01,m+1);
    xq(lett+1-(let0-let01)+1:lett+1,m+1) = xq1(m+1)*ones(let0-let01,1);
end

xq(:,3) = xc(1)*ones(length(yq(:,3)),1);

    jp=1;
    jm=let1;
    qetr=zeros(jp,1);
    yy7=qetr;
    fctr=smgvr*sqrt(0.5);
    sho=fctr*she*sqrt(2);
    [aaa,ddd] = gridf(yy7,qetr,yc(1),yc(2),sho,sho,lett,jm,jp);
    yy7 = aaa;
    qetr = ddd;
    xx7 = xc(1) * ones(length(yy7),1);
    xqr(:,3) = [xq(1:let0+1,3);xx7;xq(let0+2:end,3)];
    yqr(:,3) = [yq(1:let0+1,3);yy7';yq(let0+2:end,3)];

if k==1
    figure(3)
    plot(xq(:,2), yq(:,2), 'k', 'LineWidth', 2);
    plot(xq(:,3), yq(:,3), 'k', 'LineWidth', 2);
end

%% Left and right boundaries
fctr=2*smgvr*sqrt(0.5);
for m=0:3:3
    switch m
        case 0
            yo1=yb;
            yo2=yb;
            sho=fctr*shs;
        case 3 
            yo1=yc(2);
            yo2=yc(1)-spy;
            sho=fctr*she*sqrt(2)/2;
    end
    jj=ptvr*let01;
    ra0=fracvr*domh1;
    
    jp=let01-jj+1;
    jm=jj;
    [aaa,ddd] = gridf(yq(:,m+1),qet,yo2-ra0,yo2,free,sho,lett,jm,jp);
    yq(:,m+1) = aaa;
    qet = ddd;
    
    jp=0+1;
    jm=let01-jj; 
    [aaa,ddd] = gridf(yq(:,m+1),qet,ya,yo2-ra0,free,qet(jp+jm),lett,jm,jp);
    yq(:,m+1) = aaa;
    qet = ddd;
    
    jp=lets11+1; 
    jm=jj; 
    [aaa,ddd] = gridf(yq(:,m+1),qet,yo1,yo1+ra0,sho,free,lett,jm,jp);
    yq(:,m+1) = aaa;
    qet = ddd;
    
    jp=jp+jm;
    jm=let01-jj;
    [aaa,ddd] = gridf(yq(:,m+1),qet,yo1+ra0,yd,qet(jp-1),free,lett,jm,jp);
    yq(:,m+1) = aaa;
    qet = ddd;
    
    jp = 0+1;
    jm = let0-let01;
    [aaa,ddd] = gridf(yq1(:,m+1),qet,-domh,ya,free,qet(1),lett,jm,jp);
    yq1(:,m+1) = aaa;
    qet1 = ddd;
    
    jp = lett+1 - (let0-let01);
    jm = let0-let01;
    [aaa,ddd] = gridf(yq1(:,m+1),qet,yd,domh,qet(lets11+1+jj+let01-jj),free,lett,jm,jp);
    yq1(:,m+1) = aaa;
    qet1 = ddd;
    
end

fctr=2*skew/(yd-ya);
for j=1:lett+1
    xq(j,1)=xa-skew+fctr*(yq(j,1)-ya);
    xq(j,4)=xe-skew+fctr*(yq(j,4)-ya);
end

yq2=yq;
for m=0:3:3
    yq(1:let0-let01,m+1) = yq1(1:let0-let01,m+1);
    yq(let0-let01+1:let0+1,m+1) = yq2(1:let01+1,m+1);
    yq(let0+2:let0+2+let01,m+1) = yq2(lets11+1:lets11+1+let01,m+1);
    yq(lett+1-(let0-let01)+1:lett+1,m+1) = yq1(lett+1-(let0-let01)+1:lett+1,m+1);
end

jp=1;
jm=let1;
qetr=zeros(jp,1);
fctr=2*smgvr*sqrt(0.5)/2;
sho=fctr*she*sqrt(2);
[aaa,ddd] = gridf(yy7,qetr,yc(1),yc(2),sho,sho,lett,jm,jp);
yy7 = aaa;
qetr = ddd;

xx7 = xe * ones(length(yy7),1);
xqr(:,4) = [xq(1:let0+1,4);xx7;xq(let0+2:end,4)];
yqr(:,4) = [yq(1:let0+1,4);yy7';yq(let0+2:end,4)];

if k==1
    figure(3)
    plot(xq(:,1), yq(:,1), 'k', 'LineWidth', 2);
    plot(xq(:,4), yq(:,4), 'k', 'LineWidth', 2);
    grid off;
end

%% Grid output
for zz=1:7
    switch zz
        case 1
            n=0; m=0; js=1; je=lete0+1; is=1; ie=lxie0+1; jsr=js; jer=je;
            xleft = xq; yleft = yq; xright = xq; yright = yq;
            xxx = zeros(ie-is+1, je-js+1); yyy = zeros(ie-is+1, je-js+1);
        case 2
            n=0; m=1; js=1; je=lete0+1; is=lxis1+1; ie=lxie1+1; jsr=js; jer=je;
            xleft = xq; yleft = yq; xright = xqr; yright = yqr;
            xxx = zeros(ie-is+1, je-js+1); yyy = zeros(ie-is+1, je-js+1);
        case 3
            n=0; m=2; js=1; je=lete0+1; is=lxis2+1; ie=lxit+1; jsr=js; jer=je;
            xleft = xqr; yleft = yqr; xright = xqr; yright = yqr;
            xxx = zeros(ie-is+1, je-js+1); yyy = zeros(ie-is+1, je-js+1);
        case 4
            n=2; m=0; js=lets1+1; je=lett+1; is=1; ie=lxie0+1; jsr=js; jer=je;
            xleft = xq; yleft = yq; xright = xq; yright = yq;
            xxx = zeros(ie-is+1, je-js+1); yyy = zeros(ie-is+1, je-js+1);
        case 5
            n=2; m=1; js=lets1+1; je=lett+1; is=lxis1+1; ie=lxie1+1; jsr=lets2+1; jer=lettr+1;
            xleft = xq; yleft = yq; xright = xqr; yright = yqr;
            xxx = zeros(ie-is+1, je-js+1); yyy = zeros(ie-is+1, je-js+1);
        case 6
            n=2; m=2; js=lets2+1; je=lettr+1; is=lxis2+1; ie=lxit+1; jsr=js; jer=je;
            xleft = xqr; yleft = yqr; xright = xqr; yright = yqr;
            xxx = zeros(ie-is+1, je-js+1); yyy = zeros(ie-is+1, je-js+1);
        case 7
            n=1; m=2; js=lets1+1; je=lete1+1; is=lxis2+1; ie=lxit+1; jsr=js; jer=je;
            xleft = xqr; yleft = yqr; xright = xqr; yright = yqr;
            xxx = zeros(ie-is+1, je-js+1); yyy = zeros(ie-is+1, je-js+1);
    end
    xtop = xp;
    ytop = yp;
    xbot = xp;
    ybot = yp;
    for j=js:je
        xx(1,j-js+1) = xleft(j,[m+1]);
        yy(1,j-js+1) = yleft(j,[m+1]);
    end
    for j=jsr:jer
        xx(ie-is+1,j-jsr+1) = xright(j,[m+1+1]);
        yy(ie-is+1,j-jsr+1) = yright(j,[m+1+1]);
    end
    for i=is:ie
        xx(i-is+1,1) = xbot(i,[n+1]);
        yy(i-is+1,1) = ybot(i,[n+1]);
    end
    for i=is:ie
        xx(i-is+1,je-js+1) = xtop(i,[n+1+1]);
        yy(i-is+1,je-js+1) = ytop(i,[n+1+1]);
    end
    xco(1,:) = [xx(1,1), xx(ie-is+1,1), xx(ie-is+1,je-js+1), xx(1,je-js+1)];
    yco(1,:) = [yy(1,1), yy(ie-is+1,1), yy(ie-is+1,je-js+1), yy(1,je-js+1)];
    
    xlefts = xleft(js+1:je-1,:);
    ylefts = yleft(js+1:je-1,:);
    xrights = xright(jsr+1:jer-1,:);
    yrights = yright(jsr+1:jer-1,:);
    for j=1:(je-1-(js+1)+1)
        for i=is+1:ie-1
            xcij(:,1) = [xbot(i,n+1), xbot(i,n+1), xtop(i,n+1+1), xtop(i,n+1+1)] - xco(1,:);
            xcij(:,2) = [xlefts(j,m+1), xrights(j,m+1+1), xrights(j,m+1+1), xlefts(j,m+1)] - xco(1,:);
            ycij(:,1) = [ybot(i,n+1), ybot(i,n+1), ytop(i,n+1+1), ytop(i,n+1+1)] - yco(1,:);
            ycij(:,2) = [ylefts(j,m+1), yrights(j,m+1+1), yrights(j,m+1+1), ylefts(j,m+1)] - yco(1,:);
            jj=js+j;
            xt(1:4,1) = [(i-is)*(jj-js), (ie-i)*(jj-js), (ie-i)*(je-jj), (i-is)*(je-jj)].^2;
            xt(1,2) = xt(2,1) * xt(3,1) * xt(4,1);
            xt(2,2) = xt(3,1) * xt(4,1) * xt(1,1);
            xt(3,2) = xt(4,1) * xt(1,1) * xt(2,1);
            xt(4,2) = xt(1,1) * xt(2,1) * xt(3,1);
            res = 1/sum(xt(1:4,2));
            cha(1:4) = xcij(:,1) + xcij(:,2) + xco(1,:)';
            dha(1:4) = ycij(:,1) + ycij(:,2) + yco(1,:)';
            ii=i-is+1;
            xx(ii,j+1) = res * sum(cha(1:4)*xt(1:4,2));
            yy(ii,j+1) = res * sum(dha(1:4)*xt(1:4,2));
        end
    end
    xxx = xx(1:ie-is+1, 1:je-js+1); yyy = yy(1:ie-is+1, 1:je-js+1);
    if k==kplot
        figure(4); hold on; title('Full mesh grid');
        for i=1:ie-is+1
            plot(xxx(i,1:je-js+1), yyy(i,1:je-js+1), 'k');
        end
        for j=1:je-js+1
            plot(xxx(1:ie-is+1,j), yyy(1:ie-is+1,j), 'k');
        end
    end
    if gridsave==1
        filename_output1 = strcat('output/', filename_output);
        filename_output1 = strcat(filename_output1, int2str(zz-1));
        filename_output1 = strcat(filename_output1, '.dat');
        fileID = fopen(filename_output1,'at');
        for j=1:(je-js+1)
            for i=1:(ie-is+1)
                fprintf(fileID,'%.15e %.15e\n',xxx(i,j), yyy(i,j));
            end
        end
        fclose(fileID);
    end
    if gridtec==1
        xxxx{zz}(:,:,k) = xxx(:,:); yyyy{zz}(:,:,k) = yyy(:,:);
    end
end
end

%% Save Tecplot file
if gridtec==1
    filename_tecplot = input('Please enter the grid output filename:', 's');
    filename_tecplot = strcat('output/', filename_tecplot);
    fileTEC = fopen(filename_tecplot,'w');
    lho=2; mq=3;
    gridmino = zeros(3,1);
    gridmaxo = zeros(3,1);
    gridmindel = zeros(3,mbk);
    gridmaxdel = zeros(3,mbk);
    for zz=1:mbk
        gridmindel(1,zz) = min(min(min(xxxx{zz})));
        gridmindel(2,zz) = min(min(min(yyyy{zz})));
        gridmaxdel(1,zz) = max(max(max(xxxx{zz})));
        gridmaxdel(2,zz) = max(max(max(yyyy{zz})));
    end
    gridmino(1) = min(gridmindel(1,:));
    gridmino(2) = min(gridmindel(2,:));
    gridmino(3) = min(zs);
    gridmaxo(1) = max(gridmaxdel(1,:));
    gridmaxo(2) = max(gridmaxdel(2,:));
    gridmaxo(3) = max(zs);
    for zz=1:mbk
        lho=techead_write(fileTEC, zz-1, mbk-1, mq, lho, lximb, letmb, lzemb, 'grid', ["x","y","z"], 0.0, gridmino, gridmaxo);
        switch zz
            case 1
                js=1; je=lete0+1; is=1; ie=lxie0+1;
            case 2
                js=1; je=lete0+1; is=lxis1+1; ie=lxie1+1;
            case 3
                js=1; je=lete0+1; is=lxis2+1; ie=lxit+1;
            case 4
                js=lets1+1; je=lett+1; is=1; ie=lxie0+1;
            case 5
                js=lets1+1; je=lett+1; is=lxis1+1; ie=lxie1+1;
            case 6
                js=lets2+1; je=lettr+1; is=lxis2+1; ie=lxit+1;
            case 7
                js=lets1+1; je=lete1+1; is=lxis2+1; ie=lxit+1;
        end
        for mm=1:mq
           for k=1:lze0+1
                for j=1:je-js+1
                    for i=1:ie-is+1
                        switch mm
                            case 1
                                fwrite(fileTEC, xxxx{zz}(i,j,k), 'single');
                            case 2
                                fwrite(fileTEC, yyyy{zz}(i,j,k), 'single');
                            case 3
                                fwrite(fileTEC, zs(k), 'single');
                        end
                        lho = lho + 1;
                    end
                end
            end
        end
    end
    fclose(fileTEC);
end