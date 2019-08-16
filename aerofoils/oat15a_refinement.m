clear all
clc
close all

lxis1 = 1;
lxie1 = 300;
lxit = lxie1;
lxi1 = lxie1;
xy = dlmread('oat15a.dat');
lnaca = length(xy)/2;
xnaca(:,1) = xy(1:lnaca,1);
xnaca(:,2) = xy(lnaca+1:end,1);
ynaca(:,1) = xy(1:lnaca,2);
ynaca(:,2) = xy(lnaca+1:end,2);
figure(1)
plot(xnaca(:,1),ynaca(:,1), '-ok'); hold on; grid on;
plot(xnaca(:,2),ynaca(:,2), '-ok');

% ll=5;
% nll=10;
% deg=2;
% xxnaca(:,1) = linspace(xnaca(1,1), xnaca(ll,1), nll);
% p = polyfit(xnaca(1:ll,1),ynaca(1:ll,1), deg);
% yynaca(:,1) = polyval(p,xxnaca(:,1));
% figure(1)
% plot(xxnaca(:,1),yynaca(:,1), 'or'); hold on; grid on;

ll=5;
nll=31;
deg=4;
yyle(:) = linspace(ynaca(ll,2), ynaca(ll,1), nll);
yy = [ynaca(ll:-1:1,2), ynaca(1:ll,1)];
xx = [xnaca(ll:-1:1,2), xnaca(1:ll,1)];
p = polyfit(yy,xx, deg);
xxle(:) = polyval(p,yyle(:));
% figure(1)
% plot(xxle(:),yyle(:), 'og'); hold on; grid on;

llte=4;
nllte=5;
degte=3;
xxte(:,1) = linspace(xnaca(end-llte+1,1), xnaca(end,1), nllte);
xx=xnaca(end-llte+1:end,1); yy=ynaca(end-llte+1:end,1);
p=polyfit(xx,yy,degte);
yyte(:,1) = polyval(p,xxte(:,1));
xxte(:,2) = linspace(xnaca(end-llte+1,2), xnaca(end,2), nllte);
xx=xnaca(end-llte+1:end,2); yy=ynaca(end-llte+1:end,2);
p=polyfit(xx,yy,degte);
yyte(:,2) = polyval(p,xxte(:,2));
% figure(1)
% plot(xxte(:,1),yyte(:,1), 'or'); hold on; grid on;
% plot(xxte(:,2),yyte(:,2), 'or');

% [AAA, kk] = min(xxle);
% kkk=min(kk-1, nll-kk);
% flag=0; iii=0;
% while flag==0
%     iii = iii+1;
%     if xxle(kk-kkk)<xnaca(iii,2)
%         flag=1;
%     end
% end
% lnaca = 80;
% a = xnaca(ll,1);
% b = xnaca(end-llte,1);
% for k=1:lnaca
%     xxnaca(k,1) = (a+b)/2 + (b-a)/2 * cos((k-1)/lnaca*pi);
% end
% % xxnaca(:,1) = linspace(a, b, lnaca);
% a = xxle(kk-kkk);
% b = xnaca(end-llte,2);
% for k=1:lnaca
%     xxnaca(k,2) = (a+b)/2 + (b-a)/2 * cos((k-1)/lnaca*pi);
% end
% % xxnaca(:,2) = linspace(a, b, lnaca);
% yynaca(:,1) = spline(xnaca(ll+1:end-llte-1,1), ynaca(ll+1:end-llte-1,1), xxnaca(:,1));
% bbbx = [xxle(kk-kkk); xnaca(iii:end-llte-1,2)];
% bbby = [yyle(kk-kkk); ynaca(iii:end-llte-1,2)];
% yynaca(:,2) = spline(bbbx, bbby, xxnaca(:,2));
% % figure(1)
% % plot(xxnaca(:,1),yynaca(:,1), 'ob'); hold on; grid on;
% % plot(xxnaca(:,2),yynaca(:,2), 'ob');
xxnaca(:,1) = xnaca(ll+1:end-llte,1);
xxnaca(:,2) = xnaca(ll+1:end-llte,2);
yynaca(:,1) = ynaca(ll+1:end-llte,1);
yynaca(:,2) = ynaca(ll+1:end-llte,2);

lnaca = length(xnaca)-llte - (ll+1) +1;
[AAA, kk] = min(xxle);
kkk=min(kk-1, nll-kk);
xnaca = zeros(lnaca+kkk+1+nllte,2);
ynaca = zeros(lnaca+kkk+1+nllte,2);
xnaca(1,1) = xxle(kk);
xnaca(1,2) = xxle(kk);
ynaca(1,1) = yyle(kk);
ynaca(1,2) = yyle(kk);
xnaca(2:kkk+1,1) = xxle(kk+1:kkk+kk);
xnaca(2:kkk+1,2) = xxle(kk-1:-1:kk-kkk);
ynaca(2:kkk+1,1) = yyle(kk+1:kkk+kk);
ynaca(2:kkk+1,2) = yyle(kk-1:-1:kk-kkk);
xnaca(end-nllte+1:end,1) = xxte(:,1);
ynaca(end-nllte+1:end,1) = yyte(:,1);
xnaca(end-nllte+1:end,2) = xxte(:,2);
ynaca(end-nllte+1:end,2) = yyte(:,2);
xnaca(kkk+2:end-nllte,:) = xxnaca(1:end,:);
ynaca(kkk+2:end-nllte,:) = yynaca(1:end,:);

% figure(1)
% plot(xnaca(:,1),ynaca(:,1), '-or'); hold on; grid on;
% plot(xnaca(:,2),ynaca(:,2), '-ob');

xend(:,1) = [xnaca(1:kkk+1,1); xnaca(kkk+2:end-nllte,1); xnaca(end-nllte+1:end,1)];
xend(:,2) = [xnaca(1:kkk+1,2); xnaca(kkk+2:end-nllte,2); xnaca(end-nllte+1:end,2)];
yend(:,1) = [ynaca(1:kkk+1,1); ynaca(kkk+2:end-nllte,1); ynaca(end-nllte+1:end,1)];
yend(:,2) = [ynaca(1:kkk+1,2); ynaca(kkk+2:end-nllte,2); ynaca(end-nllte+1:end,2)];

figure(1)
plot(xend(:,1),yend(:,1), '-or'); hold on; grid on;
plot(xend(:,2),yend(:,2), '-ob')

fileID = fopen('oat15a_ref.dat','w');
for i=1:length(xend(:,1))
    fprintf(fileID,'%e %e\n', xend(i,1), yend(i,1));
end
for i=1:length(xend(:,1))
    fprintf(fileID,'%e %e\n', xend(i,2), yend(i,2));
end
fclose(fileID);