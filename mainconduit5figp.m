clear variables

global name_parameters

name_parameters='villarrica2015f';

global sol zsol vinicial wr dl count rho_ti;
global rbub rbub2 visctot visctot2 xtot1 xtot;
global wused wex deltPused Tused h2oused xused wtrans

figure(1);
set(gca, 'FontSize', 14, 'LineWidth', 2);


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

Tused=Temp(2,1);
xused=cryst(2,1);
h2oused=water(2,1);

%Tused=900 + 273.15;
%xused=0.25;
%h2oused=4;

radrub4

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

%deltPused=deltP;
%wused=w;
wtrans=wused;

RIconduitex5_4crange

if count>=39
    solex=false;
else
    solex=true;
end



if solex && strcmp('cylinder',geometry)
  critara2
end

if strcmpi(geometry,'dyke')
    Qex=vinicial*wr*dl;
else
    Qex=vinicial*3.1415*wr*wr;
end


vex=vinicial;
zsolex=zsol;
soluex=sol;
rbubex=rbub;
viscex=visctot;
xex=xtot1;


if isreal(w2) && (~isnan(deltP2) && deltT>10)
        wr=w2;
    else
        wr=w;
    end

RIconduitef5_4csmrange

if count>=49
    solef=false;
else
    solef=true;
end

if strcmpi(geometry,'dyke')
    Qef=vinicial*wr*dl;
else
    Qef=vinicial*3.1415*wr*wr;
end


vef=vinicial;
zsolef=zsol;
soluef=sol;
rbubef=rbub2;
viscef=visctot2;
xef=xtot;
  


if solef

subplot(1,4,1)
semilogx(soluef(:,1),zsolef,'linewidth',2)
xlabel('Pressure (Pa)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,4,2)
plot(soluef(:,2),zsolef,'linewidth',2)
xlabel('gas volume fraction','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,4,3)
semilogx(soluef(:,5),zsolef,soluef(:,6),zsolef,'--','linewidth',2)
xlabel('velocity (m/s)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
dim = [.1328 .2963 .1411 .0322];
str2 = ['Eff mass flow rate = ' num2str(Qef*rho_ti,'%.3g ') ' kg/s'];
str22 = ['(' num2str(Qef,'%.3g ') 'm^3/s)'];
str33 = {str2, str22};
 
annotation('textbox',dim,'String',str33,'FitBoxToText','on','fontsize',12);
hold on
  
subplot(1,4,4)
semilogx(viscef,zsolef,'linewidth',2)

dim = [.1391 .2026 .1036 .0322];
str2 = [' Conduit radius = ' num2str(wused,'%.3g ') ' m'];
%annotation('textbox',dim,'String',str2,'FitBoxToText','on','fontsize',12);
dim = [.1381 .16 .1359 .0322];
str19 = [' Inlet Overpressure = ' num2str(deltPused,'%.3g ') ' Pa'];
str20= {str2,str19};
dim10 = [.1381 .56 .1359 .0322];
str10 = [' Melt temperature = ' num2str(Tused-273.15,'%.3g ') ' °C'];
str11 = [' Initial crystal content = ' num2str(xused*100,'%.3g ') ' %'];
str12 = [' Water content = ' num2str(h2oused,'%.3g ') ' %'];
str13 = {str10, str11, str12};
annotation('textbox',dim,'String',str20,'FitBoxToText','on','Backgroundcolor','w','fontsize',12);
annotation('textbox',dim10,'String',str13,'FitBoxToText','on','Backgroundcolor','w','fontsize',12);

xlabel('viscosity (Pa.s)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
dim5 = [.4571 .9357 .1164 .0477];
str5=namevolc;
annotation('textbox',dim5,'String',str5,'FitBoxToText','on','fontsize',24,'fontweight','bold','Edgecolor','none','LineWidth',2, 'FontName', 'Arial');
hold on

end

if solex
    
subplot(1,4,1)
semilogx(soluex(:,1),zsolex,'linewidth',2)
xlabel('Pressure (Pa)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,4,2)
plot(soluex(:,2),zsolex,'linewidth',2)
xlabel('gas volume fraction','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
dim5 = [.4571 .9357 .1164 .0477];
str5=namevolc;
annotation('textbox',dim5,'String',str5,'FitBoxToText','on','fontsize',24,'fontweight','bold','Edgecolor','none','LineWidth',2, 'FontName', 'Arial');
hold on

subplot(1,4,3)
semilogx(soluex(:,5),zsolex,soluex(:,6),zsolex,'--','linewidth',2)
xlabel('velocity (m/s)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
dim = [.1308 .2504 .1474 .0322];
str2 = [' Exp mass flow rate = ' num2str(Qex*rho_ti,'%.3g ') ' kg/s'];
annotation('textbox',dim,'String',str2,'FitBoxToText','on','fontsize',12);
hold on

subplot(1,4,4)
semilogx(viscex,zsolex,'linewidth',2)
dim = [.1391 .2026 .1036 .0322];
str2 = [' Conduit radius = ' num2str(wex,'%.3g ') ' m'];
%annotation('textbox',dim,'String',str2,'FitBoxToText','on','fontsize',12);
dim = [.1381 .16 .1359 .0322];
str19 = [' Inlet Overpressure = ' num2str(deltPused,'%.3g ') ' Pa'];
str20= {str2,str19};
dim10 = [.1381 .56 .1359 .0322];
str10 = [' Melt temperature = ' num2str(Tused-273.15,'%.3g ') ' °C'];
str11 = [' Initial crystal content = ' num2str(xused*100,'%.3g ') ' %'];
str12 = [' Water content = ' num2str(h2oused,'%.3g ') ' %'];
str13 = {str10, str11, str12};
annotation('textbox',dim,'String',str20,'FitBoxToText','on','Backgroundcolor','w','fontsize',12);
annotation('textbox',dim10,'String',str13,'FitBoxToText','on','Backgroundcolor','w','fontsize',12);
xlabel('viscosity (Pa.s)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

end






