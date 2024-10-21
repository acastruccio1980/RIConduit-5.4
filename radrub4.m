
global name_parameters

eval(name_parameters);


Tc=T-273.15;

Tf=1050 - (sio2tot + h2o - 50)*12;

Tres=Tc;
H=-H;
deltT=Tc - Tf;

visc=((1-xi/phimax)^(-2.5))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,h2o/2,f2o,Tc);

dTdx=Tres/H;

if strcmpi(composition,'silicic')
c=1300; %1500 basalts, 1300 rhyolites
k=0.6e-6; %0.3e6 basalts, 0.6e6 rhyolites
L=2e5; %5e5 basalts, 2e5 rhyolites

else
    c=1500; %1500 basalts, 1300 rhyolites
    k=0.3e-6; %0.3e6 basalts, 0.6e6 rhyolites
    L=5e5; %5e5 basalts, 2e5 rhyolites
end

beta=0.15;
G=1e10;
v=0.25;
T=1e6;

f1=L/(L + (4/pi)*c*deltT);

lo=deltT/dTdx;


c1=1/0.27;
c2=sqrt((c*deltT/(L+c*deltT))/(1+sqrt((L+c*deltT)/(c*deltT))));
c3=sqrt(3*visc*k)*(G^2);


deltP=(2*c*dTdx*sqrt(3*visc*k)*(G^2)/(L*beta*sqrt(pi)))^(2/5);


deltP2=(c1*c2*c3/lo)^(2/5);

pcal=25.653e6;

lcalc=c1*c2*c3/(pcal^5/2);

w=(deltP-T)*2*(1-v^2)*H/G;

w2=(deltP2-T)*2*(1-v^2)*H/G;
