function [relvisc]=fvrel(model,cx,cxi,ar1,ar2,cmax,sr)

% NOTE THAT NOT ALL MODELS USE ALL THE INPUT PARAMETERS

% model: ER52=Einstein Roscoe 1952, vona11=Vona et al. 2011, costa09=Costa
% et al. 2009, costa09v11=Costa et al., 2009 but with parameters of Vona et
% al., 2011. EM2s=effective medium with phenocrysts and microlites

%cx= total cristal content

%cxi=phenocryst content

%ar1=aspect ratio of phenocrysts (width/length)

%ar1=aspect ratio of microlites (width/length)

%cmax=maximum packing fraction

%sr=strain rate (1/s)

if strcmpi(model,'ER52')
    relvisc=(1-cx/cmax)^-2.5;
end

if strcmpi(model,'vona11')

phimax1=0.656*exp(-(log10(ar1)^2)/(2*(1.08^2)));
phimax2=0.656*exp(-(log10(ar2)^2)/(2*(1.08^2)));
phimax=(cxi*phimax1 + (cx-cxi)*phimax2)/cx;
relvisc=(1-cx/phimax)^(-2*(1-0.06*log10(sr)));

end

if strcmpi(model,'EM2s')
    
phimax1=0.656*exp(-(log10(ar1)^2)/(2*(1.08^2)));
phimax2=0.656*exp(-(log10(ar2)^2)/(2*(1.08^2)));
relvisc1=(1-cxi/phimax1)^(-2.5);
relvisc2=(1-(cx-cxi)/((1-cxi)*phimax2))^(-2.5);
relvisc=relvisc1*relvisc2;
end

if strcmpi(model,'costa09')
    
phimax=0.066499*tanh(0.913424*log10(sr) + 3.850623) + 0.591806;
delta=-6.301095*tanh(0.818496*log10(sr) + 2.86) + 7.462405;
alpha=-0.000378*tanh(1.148101*log10(sr) + 3.92) + 0.999572;
gamma=3.987815*tanh(0.8908*log10(sr) + 3.24) + 5.099645;

f1= 1 + (cx/phimax)^delta;
f2=1 + (cx/phimax)^gamma;
f3=alpha*erf(sqrt(pi)*cx*f2/(2*alpha*phimax));
relvisc=f1/((1-f3)^(2.5*phimax));
    
end


if strcmpi(model,'costa09v11')
    
phimax=0.274;
delta=13 - 0.84;
alpha=1-0.0327;
gamma=0.84;

f1= 1 + (cx/phimax)^delta;
f2=1 + (cx/phimax)^gamma;
f3=alpha*erf(sqrt(pi)*cx*f2/(2*alpha*phimax));
relvisc=f1/((1-f3)^(2.5*phimax));
   
end


end



