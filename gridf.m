function [x,xxi] = gridf(x,pxi,xo,xn,dxo,dxn,lxi,mxin,ip)
 xxi=pxi;
 free=1e6;
 dxoo=dxo;
 dxnn=dxn;
 xi=mxin;
 sml=1e-8;
 if(int64((dxo+sml)/free)==1)
    dxoo=2*(xn-xo)/xi-dxnn;
 end
 if(int64((dxn+sml)/free)==1)
    dxnn=2*(xn-xo)/xi-dxoo;
 end
    aa=6.0*(xn-xo)-3.0*xi*(dxoo+dxnn);
    bb=15.0*(xo-xn)+xi*(8.0*dxoo+7.0*dxnn);
    cc=10.0*(xn-xo)-xi*(6.0*dxoo+4.0*dxnn);
    dd=xi*dxoo;
    fctr=1/xi;
 for i=0:mxin
    ii=i+ip;
    xi=i*fctr;
    x(ii)=aa*xi^5+bb*xi^4+cc*xi^3+dd*xi+xo;
    xxi(ii)=fctr*(5*aa*xi^4+4*bb*xi^3+3*cc*xi^2+dd);
 end

end