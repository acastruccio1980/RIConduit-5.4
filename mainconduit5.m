clear variables

global name_parameters

name_parameters='calbuco2015d';

global sol zsol vinicial wr overP dl count rho_ti;
global rbub rbub2 visctot visctot2 xtot1 xtot;

figure(1);
set(gca, 'FontSize', 14, 'LineWidth', 2);

eval(name_parameters);

wr=radius1;
overP=overP1;

RIconduitex5_4c

if count>=49
    solex=false;
else
    solex=true;
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




RIconduitef5_4c

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

subplot(1,6,1)
semilogx(soluef(:,1),zsolef,'linewidth',2)
xlabel('Pressure (Pa)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,6,2)
plot(soluef(:,2),zsolef,'linewidth',2)
xlabel('gas volume fraction','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,6,3)
semilogx(soluef(:,5),zsolef,soluef(:,6),zsolef,'--','linewidth',2)
xlabel('velocity (m/s)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
dim = [.65 .5 .3 .3];
str2 = [' mass flow rate = ' num2str(Qef*rho_ti,'%.3g ') ' kg/s'];
annotation('textbox',dim,'String',str2,'FitBoxToText','on');
hold on
   
subplot(1,6,4)
semilogx(soluef(:,3),zsolef,'linewidth',2)
xlabel('Nd','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,6,5)
plot(soluef(:,4),zsolef,'linewidth',2)
xlabel('crystal content','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,6,6)
semilogx(viscef,zsolef,'linewidth',2)
xlabel('viscosity (Pa.s)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

end

if solex
    
subplot(1,6,1)
semilogx(soluex(:,1),zsolex,'linewidth',2)
xlabel('Pressure (Pa)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,6,2)
plot(soluex(:,2),zsolex,'linewidth',2)
xlabel('gas volume fraction','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,6,3)
semilogx(soluex(:,5),zsolex,soluex(:,6),zsolex,'--','linewidth',2)
xlabel('velocity (m/s)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
dim = [.65 .4 .3 .3];
str2 = [' mass flow rate = ' num2str(Qex*rho_ti,'%.3g ') ' kg/s'];
annotation('textbox',dim,'String',str2,'FitBoxToText','on');
hold on

subplot(1,6,4)
semilogx(soluex(:,3),zsolex,'linewidth',2)
xlabel('Nd','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on
    

subplot(1,6,5)
plot(soluex(:,4),zsolex,'linewidth',2)
xlabel('crystal content','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

subplot(1,6,6)
semilogx(viscex,zsolex,'linewidth',2)
xlabel('viscosity (Pa.s)','fontweight','bold','fontsize',14)
ylabel('Depth (m)','fontweight','bold','fontsize',14)
hold on

end






