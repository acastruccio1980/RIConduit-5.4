%1D conduit model to calculate eruptive parameters
%based in equations of Slezin (2003) and Kozono and Koyaguchi (2009)
%numerical procedure: differential system of equations solved numerically with matlab solver
%ode15s
%shooting method using bisection technique.
% Border conditions cited by Yoshida and Koyaguchi (1999) and Kozono and Koyaguchi (2009).
%Viscosity model: Melt: Giordano et al. (2008). 
%Crystal content: Einstein-Roscoe
%Melt density is taken as constant and calculated at the beginning with
%chemical composition, water content, temperature and pressure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Angelo Castruccio, Universidad de Chile, 2024%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results=RIconduitex5_4crangesr
tic;
global R T Tc rho_m co C1 beta phicrit H pfinal overP Pi vinicial cg wr rho_ti xi xmax  fgi Nd q wex;
global sio2 tio2 al2o3 feo mno mgo cao na2o k2o p2o5 f2o h2o;
global sol zsol count rbub visctot limphi1 limphi2 errtol dvdzf dvdr drdrbub lsup linf limperh limperl xfinal;
global name_parameters odesolver;
global deltPused  h2oused Tused xused;

eval(name_parameters);

overP=deltPused;

h2o=h2oused;
T=Tused;
Tc=T-273.15;
xi=xused;
co=h2o/100;

Pi=rcrust*g*H*(-1) + overP; %prssure at the inlet of conduit (Pa)

exi=(1-xi)*(co-C1*Pi^beta)/(1-C1*Pi^beta); %initial exsolved water
dis=C1*Pi^beta; %initial dissolved water

if dis>co
    dis=co;
end

rho_m=density(sio2,tio2,al2o3,feo,mgo,cao,na2o,k2o,dis*100,Tc,Pi/1e6); %density of melt. Constant through the conduit

if exi<=0
    rho_ti=rho_m;
    fgi=0;
    phini=0;
   
else
    fgi=(1-xi)*(co - C1*Pi^beta)/(1-C1*Pi^beta);
    rho_ti=Pi*rho_m/(rho_m*R*T*(fgi)+Pi*(1 - fgi));
    phini=1/(1 + (Pi/(fgi*R*T))*(1-fgi)/rho_m);
    
end

vinicial=1; %initial guess of velocity at the conduit inlet (m/s)

vmin=0.1; % lower limit for velocity in the shooting method (m/s)
vmax=40; % upper limit

count=1;
vsound=0.99*sqrt(R*T);

osolv=str2func(odesolver);

% Here shooting method starts

while count<50 % From experience, if the method doesnt find a solution by the 30th iteration there is no solution, so I put 50 as a safe bet
    
fprintf('Count = %d \n',count);
fprintf('Initial velocity = %d \n',vinicial);

q=vinicial*rho_ti;
velc=sqrt(15e9/rho_m);
visc=fvrel(model,xi,xi,ar1,ar2,xmax,(vinicial/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,dis*100,f2o,Tc);
dpdzcalc = (-rho_ti*(9.81 + cg*visc*vinicial/((wr^2)*rho_ti)))/(1-(vinicial^2)/(velc^2));


dpdt=-dpdzcalc*vinicial;
coefNd=log10(dpdt);
Nd=10^(1.5*coefNd + 5);

if C1*Pi^beta<co
    Nd=1e8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhol=rho_m;
phisel=0.2;


drhogdp=1/(R*T);

Pmax=(co/C1)^(1/beta);

Pmin=pfinal;

Pcalc=(Pmax+Pmin)/2;

count2=1;

while count2<50
    
    ncalc=(1-xi)*(co - C1*Pcalc^beta)/(1-C1*Pcalc^beta);
    rhogcalc=Pcalc/(R*T);
    phicalc=1/((1/ncalc-1)*rhogcalc/rhol + 1);
    
    
    if phicalc>(phisel - 0.002) && phicalc<(phisel + 0.002)
        
        break
        
    else
        if phicalc<0.198
            Pmax=Pcalc;
            Pcalc=(Pmax+Pmin)/2;
        else
            Pmin=Pcalc;
            Pcalc=(Pmax+Pmin)/2;
        end
    end
    count2=count2+1;
end

dfgdpcalc=(C1*beta*Pcalc^(beta-1))*(co-1)/((1-C1*Pcalc^beta)^2);

dphidpcalc=-(-dfgdpcalc*rhogcalc/(rhol*ncalc^2) + (1/ncalc - 1)*drhogdp/rhol)/(((1/ncalc - 1)*rhogcalc/rhol + 1)^2);

velcc=sqrt(15e9/rhol);
viscalc=fvrel(model,xi,xi,ar1,ar2,xmax,(vinicial/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,C1*(Pcalc^beta)*100,f2o,Tc);
rho_ticalc=Pcalc*rhol/(rhol*R*T*(ncalc)+Pcalc*(1 - ncalc));

dpdzcalc2 = (-rho_ticalc*(9.81 + cg*viscalc*vinicial/((wr^2)*rho_ticalc)))/(1-(vinicial^2)/(velcc^2));

dvdz=vinicial*(-dfgdpcalc*dpdzcalc2*(1-phicalc) + (1-ncalc)*dphidpcalc*dpdzcalc2)/((1-phicalc)^2);

dvdz2=vinicial/wr;

dvdz3=(dvdz + dvdz2)/2;

dvdz=dvdz3;

rbcalc=((phicalc/((4/3)*3.14159*Nd*(1-phicalc))))^(1/3);

Ca=abs(dvdz*viscalc*rbcalc/0.3);

phicrit=((lsup - linf)/2)*erf(log10(Ca)) + (lsup + linf)/2;
limphi1=((limperl - limperh)/2)*erf(log10(Ca)) + (limperl + limperh)/2;
limphi2=limphi1 + 0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exi<=0

Pcrit=(co/C1)^(1/beta);

if Pcrit<1.013e5
    Pcrit=1.013e5;
end

deltP=Pi-Pcrit;
deltH=-deltP/dpdzcalc;

Hi=H+deltH;
hspacing=500;

zetash=zeros(hspacing,1);
solaux=zeros(hspacing,6);

for p=1:1:hspacing
    zetash(p,1)=H+(p-1)*deltH/hspacing;
    solaux(p,1)=Pi+(p-1)*deltH*dpdzcalc/hspacing;
    solaux(p,2)=0;
    solaux(p,3)=Nd;
    solaux(p,4)=xi;
    solaux(p,5)=vinicial;
    solaux(p,6)=vinicial;
    
end
    
zsoladi=zetash;
sol1=solaux;


zsol=[zsoladi; Hi];
sol=[sol1; Pcrit 0 Nd xi vinicial vinicial];


if Pcrit>1.013e5

if Hi<0 && Pcrit>1.013e5


zsoladi=[zsoladi; Hi];
sol1=[sol1; Pcrit 0 Nd xi vinicial vinicial];

% Calculo punto adicional por problemas numericos
epsd=1; %delta de distancia vertical
if Hi>(-epsd)
    epsd=-Hi/10;
end

had=Hi+epsd;
Pad=Pcrit + dpdzcalc*epsd;
fgad=(1-xi)*(co - C1*Pad^beta)/(1-C1*Pad^beta);
phiad=1/(1 + (Pad/(fgad*R*T))*(1-fgad)/rho_m);
rho_gad=Pad/(R*T);
viadm= q*(1-fgad)/((1-phiad)*rho_m);
viadg= q*fgad/(phiad*rho_gad);
% fin calculo



  M=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0 ;0 0 0 0 0 0;0 0 0 0 0 0 ];   
  opts = odeset('Mass',M,'Events', @myEvent,'RelTol',errtol,'AbsTol',errtol);
 [zsol,sol]=osolv(@momentumeq,[had 0],[Pad phiad Nd xi viadm viadg],opts);

zsol=[zsoladi;zsol];
sol=[sol1;sol];
zsol=real(zsol);
sol=real(sol);

else
    zsol=[zsoladi;Hi];
    sol=[sol1; Pi+dpdzcalc*deltH 0 Nd xi vinicial vinicial];
    zsol=real(zsol);
    sol=real(sol);
end


numb=size(sol,1);

Hi1=zsol(numb);
Pi1=sol(numb,1);
phini1=sol(numb,2);

if (Hi1<0 && Pi1>pfinal) && phini1>=limphi1
    
    um1=sol(numb,5);
    ug1=sol(numb,6);
    nd1=sol(numb,3);
    x1=sol(numb,4);
    
    M=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0 ;0 0 0 0 0 0;0 0 0 0 0 0 ];   
    opts = odeset('Mass',M,'Events', @myEvent2,'RelTol',errtol,'AbsTol',errtol);
    [zsol2,sol2]=osolv(@momentumeq2,[Hi1 0],[Pi1 phini1 nd1 x1 um1 ug1],opts);

    zsol=[zsol;zsol2];
    sol=[sol;sol2];
  
end

numb=size(sol,1);
Hi2=zsol(numb);
Pi1=sol(numb,1);
phini2=sol(numb,2);

if (Hi2<0 && Pi1>pfinal) && phini2>=limphi2
    
    
    um1=sol(numb,5);
    ug1=sol(numb,6);
    nd1=sol(numb,3);
    x1=sol(numb,4);
    
        M=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0 ;0 0 0 0 0 0;0 0 0 0 0 0 ];   
    opts = odeset('Mass',M,'Events', @myEvent3,'RelTol',errtol,'AbsTol',errtol);
    [zsol2,sol2]=osolv(@momentumeq3,[Hi2 0],[Pi1 phini2 nd1 x1 um1 ug1],opts);

    zsol=[zsol;zsol2];
    sol=[sol;sol2];
  
end

numb=size(sol,1);
Hi3=zsol(numb);
Pi1=sol(numb,1);
ndfinal=sol(numb,3);
phini3=sol(numb,2);
xfinal=sol(numb,4);

if (Hi3<0 && Pi1>pfinal) && phini3>=phicrit
    
    
    opts = odeset('RelTol',errtol,'AbsTol',errtol);
    [zsol3,sol3]=osolv(@momentumeq4,[Hi3 0],[Pi1 phini3],opts);
    
    numb2=size(sol3,1);
    um=zeros(numb2,1);
    ug=zeros(numb2,1);
    nd3=zeros(numb2,1);
    xf3=zeros(numb2,1);
    
    
    n=zeros(numb2,1);
    rho_g=zeros(numb2,1);
    
    for j=1:1:numb2
     
     test=co - C1*sol3(j,1)^beta;

    if test<=0
        n(j,1)=0;
    else
        n(j,1)=(co - C1*sol3(j,1)^beta)/(1-C1*sol3(j,1)^beta);
    end
    
    rho_g(j,1)=sol3(j,1)/(R*T);
    um(j,1)=(1-n(j,1))*q/(rho_m*(1-sol3(j,2)));
    ug(j,1)=n(j,1)*q/(rho_g(j,1)*sol3(j,2));
    nd3(j,1)=ndfinal;
    xf3(j,1)=xfinal;
    
    
    end

    sol3=[sol3 nd3 xf3 um ug];
    zsol=[zsol;zsol3];
    sol=[sol;sol3];
  
end

end


else 
         M=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0 ;0 0 0 0 0 0;0 0 0 0 0 0 ];   
    opts = odeset('Mass',M,'Events', @myEvent,'RelTol',errtol,'AbsTol',errtol);
    [zsol,sol]=osolv(@momentumeq,[H 0],[Pi phini Nd xi vinicial vinicial],opts);
    zsol=real(zsol);
    sol=real(sol);
    
    numb=size(sol,1);
    Hi1=zsol(numb);
    Pi1=sol(numb,1);
    phini1=sol(numb,2);

    if (Hi1<0 && Pi1>pfinal) && phini1>=limphi1
    
    
    um1=sol(numb,5);
    ug1=sol(numb,6);
    Nd1=sol(numb,3);
    x1=sol(numb,4);
  
        M=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0 ;0 0 0 0 0 0;0 0 0 0 0 0 ];   
    opts = odeset('Mass',M,'Events',@myEvent2,'RelTol',errtol,'AbsTol',errtol);
    [zsol2,sol2]=osolv(@momentumeq2,[Hi1 0],[Pi1 phini1 Nd1 x1 um1 ug1],opts);

    zsol=[zsol;zsol2];
    sol=[sol;sol2];
  
    end
    
    numb=size(sol,1);
    Hi2=zsol(numb);
    Pi1=sol(numb,1);
    phini2=sol(numb,2);
 
if (Hi2<0 && Pi1>pfinal) && phini2>=limphi2
    
   
    um1=sol(numb,5);
    ug1=sol(numb,6);
    Nd1=sol(numb,3);
    x1=sol(numb,4);
    
        M=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0 ;0 0 0 0 0 0;0 0 0 0 0 0 ];   
    opts = odeset('Mass',M,'Events',@myEvent3,'RelTol',errtol,'AbsTol',errtol);
    [zsol2,sol2]=osolv(@momentumeq3,[Hi2 0],[Pi1 phini2 Nd1 x1 um1 ug1],opts);

    zsol=[zsol;zsol2];
    sol=[sol;sol2];
  
end

numb=size(sol,1);
Hi3=zsol(numb);
Pi1=sol(numb,1);
ndfinal=sol(numb,3);
phini3=sol(numb,2);
xfinal=sol(numb,4);

if (Hi3<0 && Pi1>pfinal) && phini3>=phicrit
    
 
    opts = odeset('RelTol',errtol,'AbsTol',errtol);
    [zsol3,sol3]=osolv(@momentumeq4,[Hi3 0],[Pi1 phini3],opts);
    
    numb2=size(sol3,1);
    um=zeros(numb2,1);
    ug=zeros(numb2,1);
    nd3=zeros(numb2,1);
    xf3=zeros(numb2,1);

    
    n=zeros(numb2,1);
    rho_g=zeros(numb2,1);
    
    for j=1:1:numb2
     
     test=(1-xi)*co - (1-xfinal)*C1*sol3(j,1)^beta;

    if test<=0
        n(j,1)=0;
    else
        n(j,1)=((1-xi)*co - (1-xfinal)*C1*sol3(j,1)^beta)/(1-C1*sol3(j,1)^beta);
    end
    
    rho_g(j,1)=sol3(j,1)/(R*T);
    um(j,1)=(1-n(j,1))*q/(rho_m*(1-sol3(j,2)));
    ug(j,1)=n(j,1)*q/(rho_g(j,1)*sol3(j,2));
    nd3(j,1)=ndfinal;
    xf3(j,1)=xfinal;
    
    
    
    end

    sol3=[sol3 nd3 xf3 um ug];
    zsol=[zsol;zsol3];
    sol=[sol;sol3];
  
end
end

numb=size(sol,1);
ugexit=sol(numb,6);

pexit=sol(numb,1);

zexit=zsol(numb);

phiexit=sol(numb,2);


if zexit>=-5 && ((pexit>=pfinal && (ugexit>0.95*vsound && ugexit<=1.05*vsound)) || ((pexit>(pfinal - 0.05e5) && pexit<(pfinal + 0.05e5))&& phiexit>=phicrit))
    break
end

if pexit<pfinal
    vmax=vinicial;
    vinicial = vmin + (vmax-vmin)/2;
    
else
    
if zexit<-5 || ugexit>1.02*vsound
    vmax=vinicial;
    vinicial = vmin + (vmax-vmin)/2;
    
else

if ugexit<=0.98*vsound 
    vmin=vinicial;
    vinicial = vmin + (vmax-vmin)/2;
end
end

end



count=count+1;
wex=wr;

if vinicial<0.2
end




end
% End shooting method


results=[zsol,sol];

numb=size(sol,1);


 n=zeros(1,numb);
 visctot=zeros(1,numb);
 captot=zeros(1,numb);
 rbtot=zeros(1,numb);
 rbub=zeros(1,numb);
 dvdzf=zeros(1,numb);
 drdrbub=zeros(1,numb);
 dvdr=zeros(1,numb);
 
  for j=1:1:numb
     
     test=(1-xi)*co - (1-sol(j,4))*C1*sol(j,1)^beta;
     
    if test<=0
        n(1,j)=0;
    else
        n(1,j)=((1-xi)*co - (1-sol(j,4))*C1*sol(j,1)^beta)/(1-C1*sol(j,1)^beta);
        
    end
   
    
    if test>0

        viscl=fvrel(model,sol(j,4),xi,ar1,ar2,xmax,(sol(j,5)/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,(C1*sol(j,1)^beta)*100,f2o,Tc);
        
    else
    
        viscl=fvrel(model,sol(j,4),xi,ar1,ar2,xmax,(sol(j,5)/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,h2o,f2o,Tc);
    
    end

rb=((sol(j,2)/((4/3)*3.14159*sol(j,3)*(1-sol(j,2)))))^(1/3);



c1=-0.2895*sol(j,2) + 0.8132;
c2=sol(j,2);
nca=rb*viscl*(sol(j,5)/wr)/3;

phicritbub=phicrit+0.05;
AA=(1-(sol(j,2)/phicritbub))^(-phicritbub);
BB=(1-(sol(j,2)/phicritbub))^(5*phicritbub/3);


viscrel=0.5*(AA-BB)*(1-erf(real(c1*log(nca)+c2)))+BB;

visc=viscrel*viscl;

    rbtot(1,j)=viscrel;
    captot(1,j)=nca;
    visctot(1,j)=visc;
    if sol(j,2)>phicrit
        visctot(1,j)=NaN;
    end
    rbub(1,j)=rb;

    if j==1

        dvdzf(1,j)=0;
        drdrbub(1,j)=0;
        dvdr(1,j)=0;

    else

        dvdzf(1,j)=(sol(j,4)-sol(j-1,4))/(zsol(j,1)-zsol(j-1,1));
        drdrbub(1,j)=(rbub(1,j)-rbub(1,j-1))*sol(j,4)/((zsol(j,1)-zsol(j-1,1))*rbub(1,j));
        dvdr(1,j)=sol(j,4)/wr;
    end

    
    
  end
 


toc;
end

function dy=momentumeq(t,y)

global R T Tc rho_m co C1 beta  cg wr  xi q vinicial tcar xmax Fc F1 F2 phicrit fragcrit;
global sio2 tio2 al2o3 feo mno mgo cao na2o k2o p2o5 f2o h2o model ar1 ar2;

test=(1-xi)*co - (1-y(4))*C1*y(1)^beta;

%calculation of fraction gas exsolved
    if test<=0
        fg=0;  
    else
        fg=(co*(1-xi) - C1*(1-(y(4)))*y(1)^beta)/((1-C1*y(1)^beta));  
    end

%calculation of volume fraction of bubbles

rho_g=y(1)/(R*T);

rb=((y(2)/((4/3)*3.14159*y(3)*(1-y(2)))))^(1/3);



if test>0   
    viscl=fvrel(model,y(4),xi,ar1,ar2,xmax,(y(5)/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,(C1*y(1)^beta)*100,f2o,Tc);
else    
    viscl=fvrel(model,y(4),xi,ar1,ar2,xmax,(y(5)/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,h2o,f2o,Tc);   
end

nca=rb*viscl*(y(5)/wr)/3;

phicritbub=phicrit+0.05;
AA=(1-(y(2)/phicritbub))^(-phicritbub);
BB=(1-(y(2)/phicritbub))^(5*phicritbub/3);

c1=-0.2895*y(2) + 0.8132;
c2=y(2);

viscrel=0.5*(AA-BB)*(1-erf(real(c1*log(nca)+c2)))+BB;

visc=viscrel*viscl;

%%%
f2cx=(co*(1-xi)-((1-y(4))*C1*y(1)^beta))/(co*(1-xi) - (1-xmax)*C1*(1.01e5)^beta);
if f2cx<0
    f2cx=0;
end

xteo=xi + (xmax-xi)*f2cx;

f3cx=1-y(4)/xteo;

if f3cx<0
    f3cx=0;
end

dxdp=(xmax-xi)*f2cx*f3cx/(tcar*y(5));
if dxdp<0
    dxdp=0;
end

dfgdp=(-(-dxdp*C1*y(1)^beta + (1-y(4))*C1*beta*y(1)^(beta-1))+(co*(1-xi)-(1-y(4))*C1*y(1)^beta)*C1*beta*y(1)^(beta-1))/((1-C1*y(1)^beta)^2);


%%%


Fmw=cg*visc*y(5)/(wr^2);

Fmg=3*visc*(y(6)-y(5))*y(2)*(1-y(2))/(rb^2);


Fgw=0;

aco=rho_g*(y(6)^2);
bco=y(2) - (y(6)^2)*y(2)/(R*T) + dfgdp*q*y(6);
cco=rho_m*(y(5)^2);
dco=dfgdp*q*y(5) - (1-y(2));
eco=-rho_m*(1-y(2))*9.8 + Fmg - Fmw;
fco=rho_g*y(2)*9.8 + Fmg + Fgw;

dphidz= (-eco*bco + dco*fco)/(aco*dco - bco*cco);
dpdz= (-eco*aco + fco*cco)/(aco*dco - bco*cco);

dvdz=vinicial*(-dfgdp*dpdz*(1-y(2)) + (1-fg)*dphidz)/((1-y(2))^2);

fragcrit=dvdz*visc/(0.01*1e10);


lim=0.5;


if rb<lim*wr && y(2)<phicrit
    
dNdt=-(y(3)^(2/3))*((1/(1-y(2)))^(1/3))*((1/9)*(rho_m-rho_g)*9.81/visc)*((3*y(2)/(4*pi))^(2/3))*(F1*F1 - F2*F2)*(1/(1-((F1+F2)*((3*y(2)*(pi/(6*phicrit))/(4*pi))^(1/3)))))*Fc*(1-y(2)/phicrit)*(wr-rb)/wr;

else
    
   dNdt=0; 
   
end

f2=((1-xi)*co-(1-y(4))*(C1*y(1)^beta))/((1-xi)*co - (1-xmax)*C1*(1.01e5)^beta);
if f2<0
    f2=0;
end

xteo=xi + (xmax-xi)*f2;

f3=1-y(4)/xteo;

if f3<0
    f3=0;
end

dxdz=(xmax-xi)*f2*f3/(tcar);
if dxdz<0
    dxdz=0;
end

dy=zeros(6,1);

dy(1)=dpdz;
dy(2)= dphidz;
dy(3)=dNdt/y(5);
dy(4)=dxdz;
dy(5)= y(5) - q*(1-fg)/((1-y(2))*rho_m);
dy(6)= y(6) - q*fg/(y(2)*rho_g);

end

function dy=momentumeq2(t,y)

global R T  Tc rho_m co C1 beta  cg wr xi q vinicial tcar xmax Fc F1 F2 phicrit fragcrit limphi1 limphi2;
global sio2 tio2 al2o3 feo mno mgo cao na2o k2o p2o5 f2o h2o model ar1 ar2;

test=(1-xi)*co - (1-y(4))*C1*y(1)^beta;

%calculation of fraction gas exsolved
    if test<=0
        fg=0;  
    else
        fg=(co*(1-xi) - C1*(1-(y(4)))*y(1)^beta)/((1-C1*y(1)^beta));  
    end

%calculation of volume fraction of bubbles

rho_g=y(1)/(R*T);

rb=((y(2)/((4/3)*3.14159*y(3)*(1-y(2)))))^(1/3);



if test>0   
    viscl=fvrel(model,y(4),xi,ar1,ar2,xmax,(y(5)/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,(C1*y(1)^beta)*100,f2o,Tc);
else    
    viscl=fvrel(model,y(4),xi,ar1,ar2,xmax,(y(5)/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,h2o,f2o,Tc);   
end

nca=rb*viscl*(y(5)/wr)/3;

phicritbub=phicrit+0.05;
AA=(1-(y(2)/phicritbub))^(-phicritbub);
BB=(1-(y(2)/phicritbub))^(5*phicritbub/3);

c1=-0.2895*y(2) + 0.8132;
c2=y(2);

viscrel=0.5*(AA-BB)*(1-erf(real(c1*log(nca)+c2)))+BB;

visc=viscrel*viscl;

%%%
f2cx=(co*(1-xi)-((1-y(4))*C1*y(1)^beta))/(co*(1-xi) - (1-xmax)*C1*(1.01e5)^beta);
if f2cx<0
    f2cx=0;
end

xteo=xi + (xmax-xi)*f2cx;

f3cx=1-y(4)/xteo;

if f3cx<0
    f3cx=0;
end

dxdp=(xmax-xi)*f2cx*f3cx/(tcar*y(5));
if dxdp<0
    dxdp=0;
end

dfgdp=(-(-dxdp*C1*y(1)^beta + (1-y(4))*C1*beta*y(1)^(beta-1))+(co*(1-xi)-(1-y(4))*C1*y(1)^beta)*C1*beta*y(1)^(beta-1))/((1-C1*y(1)^beta)^2);


tt=(y(2)-limphi1)/(limphi2-limphi1);

Fmw=cg*visc*y(5)/(wr^2);

Re=2*rb*rho_g*(y(6)-y(5))/(1e-5);

if Re>2200
    Fmg=(((0.33/(4*rb))*rho_g*abs(y(6)-y(5)))^tt)*(((3)*visc*(1/(rb^2)))^(1-tt))*(y(6)-y(5))*y(2)*(1-y(2));
else
    %kper=(rb^2)*(y(2))/8;
    kper=0.125*(rb^2)*((y(2)-limphi1 + 0.05)^2.1);
    Fmg=((1e-5/kper)^tt)*(((3)*visc*(1/(rb^2)))^(1-tt))*(y(6)-y(5))*y(2)*(1-y(2));
end



Fgw=0;

aco=rho_g*(y(6)^2);
bco=y(2) - (y(6)^2)*y(2)/(R*T) + dfgdp*q*y(6);
cco=rho_m*(y(5)^2);
dco=dfgdp*q*y(5) - (1-y(2));
eco=-rho_m*(1-y(2))*9.8 + Fmg - Fmw;
fco=rho_g*y(2)*9.8 + Fmg + Fgw;

dphidz= (-eco*bco + dco*fco)/(aco*dco - bco*cco);
dpdz= (-eco*aco + fco*cco)/(aco*dco - bco*cco);

dvdz=vinicial*(-dfgdp*dpdz*(1-y(2)) + (1-fg)*dphidz)/((1-y(2))^2);
fragcrit=dvdz*visc/(0.01*1e10);

Ca=abs(dvdz*viscl*rb/0.1);

lim=0.5;


if rb<lim*wr && y(2)<phicrit
    
dNdt=-(y(3)^(2/3))*((1/(1-y(2)))^(1/3))*((1/9)*(rho_m-rho_g)*9.81/visc)*((3*y(2)/(4*pi))^(2/3))*(F1*F1 - F2*F2)*(1/(1-((F1+F2)*((3*y(2)*(pi/(6*phicrit))/(4*pi))^(1/3)))))*Fc*(1-y(2)/phicrit)*(wr-rb)/wr;

else
    
   dNdt=0; 
   
end

f2=(co-(C1*y(1)^beta))/(co - C1*(1.01e5)^beta);
if f2<0
    f2=0;
end

xteo=xi + (xmax-xi)*f2;

f3=1-y(4)/xteo;

if f3<0
    f3=0;
end

dxdz=(xmax-xi)*f2*f3/(tcar);
if dxdz<0
    dxdz=0;
end

dy=zeros(6,1);

dy(1)=dpdz;
dy(2)= dphidz;
dy(3)=dNdt/y(5);
dy(4)=dxdz;
dy(5)= y(5) - q*(1-fg)/((1-y(2))*rho_m);
dy(6)= y(6) - q*fg/(y(2)*rho_g);

end

function dy=momentumeq3(t,y)

global R T Tc  rho_m  co C1 beta  cg wr xi q vinicial tcar xmax F1 F2 Fc phicrit fragcrit limphi1;
global sio2 tio2 al2o3 feo mno mgo cao na2o k2o p2o5 f2o h2o model ar1 ar2;

test=(1-xi)*co - (1-y(4))*C1*y(1)^beta;

%calculation of fraction gas exsolved
    if test<=0
        fg=0;  
    else
        fg=(co*(1-xi) - C1*(1-(y(4)))*y(1)^beta)/((1-C1*y(1)^beta));  
    end

%calculation of volume fraction of bubbles

rho_g=y(1)/(R*T);

rb=((y(2)/((4/3)*3.14159*y(3)*(1-y(2)))))^(1/3);



if test>0   
    viscl=fvrel(model,y(4),xi,ar1,ar2,xmax,(y(5)/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,(C1*y(1)^beta)*100,f2o,Tc);
else    
    viscl=fvrel(model,y(4),xi,ar1,ar2,xmax,(y(5)/wr))*viscosity(sio2,tio2,al2o3,feo,mno,mgo,cao,na2o,k2o,p2o5,h2o,f2o,Tc);   
end

nca=rb*viscl*(y(5)/wr)/3;

phicritbub=phicrit+0.05;
AA=(1-(y(2)/phicritbub))^(-phicritbub);
BB=(1-(y(2)/phicritbub))^(5*phicritbub/3);

c1=-0.2895*y(2) + 0.8132;
c2=y(2);

viscrel=0.5*(AA-BB)*(1-erf(real(c1*log(nca)+c2)))+BB;

visc=viscrel*viscl;

%%%
f2cx=(co*(1-xi)-((1-y(4))*C1*y(1)^beta))/(co*(1-xi) - (1-xmax)*C1*(1.01e5)^beta);
if f2cx<0
    f2cx=0;
end

xteo=xi + (xmax-xi)*f2cx;

f3cx=1-y(4)/xteo;

if f3cx<0
    f3cx=0;
end

dxdp=(xmax-xi)*f2cx*f3cx/(tcar*y(5));
if dxdp<0
    dxdp=0;
end

dfgdp=(-(-dxdp*C1*y(1)^beta + (1-y(4))*C1*beta*y(1)^(beta-1))+(co*(1-xi)-(1-y(4))*C1*y(1)^beta)*C1*beta*y(1)^(beta-1))/((1-C1*y(1)^beta)^2);


%%%
Fmw=cg*visc*y(5)/(wr^2);

Re=2*rb*rho_g*(y(6)-y(5))/(1e-5);

if Re>2200
    Fmg=(0.33/(4*rb))*rho_g*(y(6)-y(5))*(y(6)-y(5))*y(2)*(1-y(2));
else
    %kper=(rb^2)*(y(2))/8;
    kper=0.125*(rb^2)*((y(2)-limphi1 + 0.05)^2.1);
    Fmg=(1e-5/kper)*(y(6)-y(5))*y(2)*(1-y(2));
end

Fgw=0;

aco=rho_g*(y(6)^2);
bco=y(2) - (y(6)^2)*y(2)/(R*T) + dfgdp*q*y(6);
cco=rho_m*(y(5)^2);
dco=dfgdp*q*y(5) - (1-y(2));
eco=-rho_m*(1-y(2))*9.8 + Fmg - Fmw;
fco=rho_g*y(2)*9.8 + Fmg + Fgw;

dphidz= (-eco*bco + dco*fco)/(aco*dco - bco*cco);
dpdz= (-eco*aco + fco*cco)/(aco*dco - bco*cco);

dvdz=vinicial*(-dfgdp*dpdz*(1-y(2)) + (1-fg)*dphidz)/((1-y(2))^2);
fragcrit=dvdz*visc/(0.01*1e10);

Ca=abs(dvdz*viscl*rb/0.1);

lim=0.5;


if rb<lim*wr && y(2)<0.5
    
dNdt=-(y(3)^(2/3))*((1/(1-y(2)))^(1/3))*((1/9)*(rho_m-rho_g)*9.81/visc)*((3*y(2)/(4*pi))^(2/3))*(F1*F1 - F2*F2)*(1/(1-((F1+F2)*((3*y(2)*(pi/(6*phicrit))/(4*pi))^(1/3)))))*Fc*(1-y(2)/phicrit)*(wr-rb)/wr;

else
    
   dNdt=0; 
   
end

f2=(co-(C1*y(1)^beta))/(co - C1*(1.01e5)^beta);
if f2<0
    f2=0;
end

xteo=xi + (xmax-xi)*f2;

f3=1-y(4)/xteo;

if f3<0
    f3=0;
end

dxdz=(xmax-xi)*f2*f3/(tcar*y(5));
if dxdz<0
    dxdz=0;
end

dy=zeros(6,1);

dy(1)=dpdz;
dy(2)= dphidz;
dy(3)= dNdt/y(5);
dy(4)=dxdz;
dy(5)= y(5) - q*(1-fg)/((1-y(2))*rho_m);
dy(6)= y(6) - q*fg/(y(2)*rho_g);

end

function dy=momentumeq4(t,y)

global R T rho_m  co C1 beta wr Nd q phicrit xfinal xi xmax tcar;

test=(1-xi)*co - (1-xfinal)*C1*y(1)^beta;

%calculation of fraction gas exsolved
    
    if test<=0
        fg=0;
              
    else
         fg=(co*(1-xi) - C1*(1-xfinal)*y(1)^beta)/((1-C1*y(1)^beta));  
        
        
    end

%calculation of volume fraction of bubbles
rho_mf=rho_m;
rho_g=y(1)/(R*T);

ra=1e-3;
cd=0.8;

if test<=0
    dfgdp=0;
else
    


%%%

f2cx=(co*(1-xi)-(C1*(1-xfinal)*y(1)^beta))/(co*(1-xi) - (1-xmax)*C1*(1.01e5)^beta);
if f2cx<0
    f2cx=0;
end

xteo=xi + (xmax-xi)*f2cx;

f3cx=1-xfinal/xteo;

if f3cx<0
    f3cx=0;
end

dxdp=(xmax-xi)*f2cx*f3cx/(tcar);
if dxdp<0
    dxdp=0;
end

dfgdp=(-(-dxdp*C1*y(1)^beta + (1-xfinal)*C1*beta*y(1)^(beta-1))+(co*(1-xi)-(1-xfinal)*C1*y(1)^beta)*C1*beta*y(1)^(beta-1))/((1-C1*y(1)^beta)^2);


%%%
end

um=(1-fg)*q/(rho_mf*(1-y(2)));

ug=fg*q/(rho_g*y(2));

%calculation of density and viscosity
rb=((y(2)/((4/3)*3.14159*Nd*(1-y(2)))))^(1/3);

Fgw=0.01*rho_g*abs(ug)*ug/(4*wr);

if y(2)<phicrit + 0.05
    tt=(y(2)-phicrit)/(0.05);
    Fmg= ((0.33/(4*rb))^(1-tt))*((3*cd/(8*ra))^tt)*rho_g*abs(ug-um)*(ug-um)*y(2)*(1-y(2));
else
    Fmg=3*cd*rho_g*abs(ug-um)*(ug-um)*y(2)*(1-y(2))/(8*ra);
end

aco=rho_g*(ug^2);
bco=y(2) - (ug^2)*y(2)/(R*T) + dfgdp*q*ug;
cco=rho_mf*(um^2);
dco=dfgdp*q*um - (1-y(2));
eco=-rho_mf*(1-y(2))*9.8 + Fmg;
fco=rho_g*y(2)*9.8 + Fmg + Fgw;

dphidz= (-eco*bco + dco*fco)/(aco*dco - bco*cco);
dpdz= (-eco*aco + fco*cco)/(aco*dco - bco*cco);

dy=zeros(2,1);

dy(1)=dpdz;
dy(2)= dphidz;

end


function [value, isterminal, direction] = myEvent(~, sol)
global limphi1
value      = (sol(2) - limphi1);
isterminal = 1;   % Stop the integration
direction  = 0;
end

function [value, isterminal, direction] = myEvent2(~, sol)
global limphi2
value      = (sol(2) - limphi2);
isterminal = 1;   % Stop the integration
direction  = 0;
end

function [value, isterminal, direction] = myEvent3(~, sol)
global phicrit 
value = sol(2) - phicrit;

isterminal = 1;   % Stop the integration
direction  = 0;
end


