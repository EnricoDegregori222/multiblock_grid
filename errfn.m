function [err] = errfn(i,ii,n,shg,xp,yp)

err=sqrt((xp(i,n+1)-xp(ii,n+1))^2+(yp(i,n+1)-yp(ii,n+1))^2)/shg-1;

end