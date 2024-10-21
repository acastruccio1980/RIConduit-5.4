%Magma density function, based on the program of Tom D Pering - University of Sheffield

%Based on Bottinga & Weill (1970).
%Densities of liquid silicate systems calculated from partial molar volumes of oxide components.
%American Journal of Science 269, pp. 169-182

function [density_gas]=density(Wsio2,Wtio2,Wal2o3,Wfeo,Wmgo,Wcao,Wna2o,Wk2o,Wh2o,mTc,mPp) 

 totalElements=(Wsio2+Wtio2+Wal2o3+Wfeo+Wmgo+Wcao+Wna2o+Wk2o); 
 
Wsio2=(100-Wh2o)*Wsio2/totalElements;
Wtio2=(100-Wh2o)*Wtio2/totalElements;
Wal2o3=(100-Wh2o)*Wal2o3/totalElements;
Wfeo=(100-Wh2o)*Wfeo/totalElements;
Wmgo=(100-Wh2o)*Wmgo/totalElements;
Wcao=(100-Wh2o)*Wcao/totalElements;
Wna2o=(100-Wh2o)*Wna2o/totalElements;
Wk2o=(100-Wh2o)*Wk2o/totalElements;
    
Msio2=60.085; 
Mtio2=79.899; 
Mal2o3=101.961; 
Mfeo=71.846; 
Mmgo=40.304; 
Mcao=56.079; 
Mna2o=61.979; 
Mk2o=94.203; 
Mh2o=18.01528; 

%Calculations
    %Step 1 Mol Prop
a=Wsio2/Msio2; 
b=Wtio2/Mtio2; 
c=Wal2o3/Mal2o3; 
d=Wfeo/Mfeo; 
e=Wmgo/Mmgo; 
f=Wcao/Mcao;
g=Wna2o/Mna2o; 
h=Wk2o/Mk2o; 
j=Wh2o/Mh2o; 

% Molecular volume at 1400C
MVsio2=26.9;  
MVtio2=23.16;  
MVal2o3=37.11; 
MVfeo=13.65; 
MVmgo=11.45; 
MVcao=16.57; 
MVna2o=28.78; 
MVk2o=45.84;
MVh2o=17; 

%DV/DT

dvdtsio2=0;  
dvdttio2=7.24;  
dvdtal2o3=2.62; 
dvdtfeo=2.92; 
dvdtmgo=2.62; 
dvdtcao=2.92; 
dvdtna2o=7.41; 
dvdtk2o=11.91;
dvdth2o=9.46;

%DV/DP

dvdpsio2=-1.89;  
dvdptio2=-2.31;  
dvdpal2o3=-2.26; 
dvdpfeo=-0.45; 
dvdpmgo=0.27; 
dvdpcao=0.34; 
dvdpna2o=-2.4; 
dvdpk2o=-6.75;
dvdph2o=-3.15;

a4=MVsio2+dvdtsio2*0.001*(mTc-1400)+dvdpsio2*0.001*(mPp-0.1);
b4=MVtio2+dvdttio2*0.001*(mTc-1400)+dvdptio2*0.001*(mPp-0.1);
c4=MVal2o3+dvdtal2o3*0.001*(mTc-1400)+dvdpal2o3*0.001*(mPp-0.1);
d4=MVfeo+dvdtfeo*0.001*(mTc-1400)+dvdpfeo*0.001*(mPp-0.1);
e4=MVmgo+dvdtmgo*0.001*(mTc-1400)+dvdpmgo*0.001*(mPp-0.1);
f4=MVcao+dvdtcao*0.001*(mTc-1400)+dvdpcao*0.001*(mPp-0.1);
g4=MVna2o+dvdtna2o*0.001*(mTc-1400)+dvdpna2o*0.001*(mPp-0.1);
h4=MVk2o+dvdtk2o*0.001*(mTc-1400)+dvdpk2o*0.001*(mPp-0.1);
j4=MVh2o+dvdth2o*0.001*(mTc-1400)+dvdph2o*0.001*(mPp-0.1);

a5=a*a4;
b5=b*b4;
c5=c*c4;
d5=d*d4;
e5=e*e4;
f5=f*f4;
g5=g*g4;
h5=h*h4;
j5=j*j4;

suma=a5+b5+c5+d5+e5+f5+g5+h5+j5;

density_gas=100000/suma;

clear Wsio2 Wtio2 Wal2o3 Wfeo Wmgo Wcao Wna2o Wk2o Wco2 Wh2o Wso2
clear Msio2 Mtio2 Mal2o3 Mfeo Mmgo Mcao Mna2o Mk2o Mco2 Mh2o Mso2
clear MVsio2 MVtio2 MVal2o3 MVfeo MVmgo MVcao MVna2o MVk2o MVco2 MVh2o MVso2
clear a b c d e f g h  j   a4 b4 c4 d4 e4 f4 g4 h4  j4  
clear a5 b5 c5 d5 e5 f5 g5 h5 i5 j5 k5 suma
clear dvdtsio2 dvdttio2 dvdtal2o3 dvdtfeo dvdtmgo dvdtcao dvdtna2o dvdth2o
clear dvdpsio2 dvdptio2 dvdpal2o3 dvdpfeo dvdpmgo dvdpcao dvdpna2o dvdph2o

end