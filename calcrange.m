clear varialbes 
global Tused xused h2oused deltPused wused;
global name_parameters count vinicial wr dl geometry wtrans;

name_parameters='caulle2011c';

eval(name_parameters);

maxind=3;
factor=1;

Temp=zeros(maxind,1);
cryst=zeros(maxind,1);
water=zeros(maxind,1);


resultsex=zeros(maxind,maxind,maxind);
resultsef=zeros(maxind,maxind,maxind);

for i=1:1:maxind
    
    Temp(i,1)=T + factor*(i-2)*25;
    cryst(i,1)=xi + factor*(i-2)*0.05;
    water(i,1)=h2o + factor*(i-2)*0.5;
    
end


for ind1=1:1:maxind
    Tused=Temp(ind1,1);
    for ind2=1:1:maxind
        xused=cryst(ind2,1);
        for ind3=1:1:maxind
            h2oused=water(ind3,1);
            
    radrub4range
    
    
    if isreal(deltP2) && (~isnan(deltP2) && deltT>10)
        deltPused=deltP2;
    else
        deltPused=deltP;
    end



    if isreal(w2) && (~isnan(deltP2) && deltT>10)
        wused=w2;
    else
        wused=w;
    end

wtrans=wused;

 RIconduitex5_4crange

if count>=49
    solex=false;
else
    solex=true;
end

if solex && strcmp('cylinder',geometry)
  critara2 
  Qex=vinicial*3.1415*wr*wr;

else
    if solex
    Qex=vinicial*wr*dl;
    else
        Qex=0;
    end
end


RIconduitef5_4crange

if count>=49
    solef=false;
else
    solef=true;
end


if solef
    
if strcmpi(geometry,'dyke')
    Qef=vinicial*wr*dl;
else
    Qef=vinicial*3.1415*wr*wr;
end

else
    Qef=0;
    
end

resultsex(ind1,ind2,ind3)=Qex;
resultsef(ind1,ind2,ind3)=Qef;
            
            
        end
    end
end