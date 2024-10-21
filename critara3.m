global name_parameters
global sol zsol vinicial rho_m wr overP countar phicrit rho_ti;
global vmin vmax xi h2o co T Tc rcrust wex h2oused Tused xused wtrans;

figure(1);
%pos = get(gcf, 'Position');
set(gca, 'FontSize', 14, 'LineWidth', 2);

eval(name_parameters);

wr=wtrans;
countar=1;

rmax=50;
rmin=1;



while countar<40
 
    fprintf('Count = %d \n',countar);
fprintf('conduit radius = %d \n',wr);

RIconduitex5_4c



Qex=vinicial*rho_m*3.1415*wr*wr;
vex=vinicial;
zsolex=zsol;
soluex=sol;

numbz=size(zsol,1);
c=1e6;
fang=32;
vstressg=rcrust*9.81;
hstressg=0.7*vstressg;


q=(tan(pi/4 + fang*pi/360))^2;
    cosang=cos(fang*pi/180);
    sinang=sin(fang*pi/180);

cMC1=zeros(numbz);
cMC2=zeros(numbz);

for i=1:1:numbz

    sigmar=soluex(i,1);
    sigmate=-2*(hstressg*zsolex(i)) - sigmar;
    sigmaz=-(vstressg*zsolex(i));

    Pp=-(vstressg*zsolex(i))*0.2;
    A=-(3*hstressg*zsolex(i) - hstressg*zsolex(i));
    B=-(vstressg*zsolex(i));
    C=2*c*cosang/(1-sinang) - Pp*(q-1);
    K=2*c*cosang + sinang*(B - 2*Pp);
    G=K + sinang*A;
    H=(A^2)*(4*(sinang^2) - 3) + (B^2 - A*B)*(4*(sinang^2) - 12);

    if sigmaz>sigmate && sigmate>sigmar

        Pc1=(B-C)/q;
        Pc2= (1/(6-2*(sinang^2)))*((3*A + 2*sinang*K)-sqrt(H + 12*(K^2 +sinang*A*K)));

    else

    if sigmate>sigmaz && sigmaz>sigmar
        Pc1=(A-C)/(1+q);
        Pc2=A/2 - (1/6)*sqrt(12*((2*c*cosang + sinang*(A-2*Pp))^2) - 3*((A - 2*B)^2));

    else

    if sigmate>sigmar && sigmar>sigmaz
        Pc1=A-C-q*B;
        Pc2=(1/(6-2*(sinang^2)))*((3*A - 2*sinang*G) - sqrt(H + 12*(G^2  - sinang*A*G)));
    else

        Pc1=0;
        Pc2=0;
    end
    end
    end
    
dif1=soluex(i,1)-Pc1;
dif2=soluex(i,1)-Pc2;

    cMC1(i)=Pc1;
    cMC2(i)=Pc2;
radius=wr;
if i==1
    maxc1=dif1;
    maxc2=dif2;
else
    
    if dif1<maxc1
    maxc1=dif1;
    end

if dif2<maxc2
    maxc2=dif2;
end
    

end

end



maxc1
maxc2

wex=wr;

if (maxc1<=1e5 && maxc1>=-1e5) || (maxc1>=-1e5 && countar==1)

    break

        
else
    if maxc1>1e5
        rmax=wr;
        wr=(rmax + rmin)/2;
        countar=countar + 1
    else
        rmin=wr;
        wr=(rmax + rmin)/2;
        countar=countar + 1
    end
end


    
end

if wr<wtrans 
    wr=wtrans;
    
end



    





