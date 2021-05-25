% Juan Camilo Gutierrez 
% 201710009014

clc
clear
close all

% Este programa calcula el input para el programa IBEM del profesor Juan
% Diego Jaramillo, en el caso de una viga prismatica en voladizo con una
% carga en el extremo izquierdo y empotrada en el extremo derecho

%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
% Presión manometrica aplicada N/mm^2
P = 250;

% Radio de tanque mm
R = 5000; 

% Espesor de tanque cilindrico mm
t = 500;

%Esfuerzo principal
Sigma1 = P*R/t;
sprintf('Esfuerzo principal analitico es %d', Sigma1)

% Numero de nodos por lado ingresado
n = 20; 

E = 72400;
Nu = 0.33;
G = 28000.00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WCS = eye(3);

Xin = zeros(size(t));
Yin = zeros(size(t));
Xout = zeros(size(t));
Yout = zeros(size(t));

Xnorm = zeros(size(t));
Ynorm = zeros(size(t));
Pxnorm = zeros(size(t));
Pynorm = zeros(size(t));

cont = 1;
for i= 0:2*pi/n:2*pi
    i2 = 2*pi - i;
    Xin(cont) = R*cos(i); 
    Yin(cont) = R*sin(i);  
    Xout(cont) = (R+t)*cos(i2);
    Yout(cont) = (R+t)*sin(i2); 
    
    Xnorm(cont) = cos(i+2*pi/n/2); %deriv cos = -sin -> sin (derivo para tangente, y alterno cambiando uno de signo)
    Ynorm(cont) = sin(i+2*pi/n/2); %deriv sin = cos -> cos
    Pxnorm(cont) = ( R*cos(i) + R*cos(i+2*pi/n) ) /2;
    Pynorm(cont) = ( R*sin(i) + R*sin(i+2*pi/n) ) /2;
    
    cont = cont + 1;   
end

figure(1)
subplot(2,2,1)
hold on
o = quiver(Pxnorm(1:n),Pynorm(1:n),Xnorm(1:n),Ynorm(1:n),'r','AutoScaleFactor',0.5,'LineWidth',1.2);
plot(Xin,Yin, 'b','LineWidth',2)
plot(Xout,Yout, 'k','LineWidth',1.2)
plt_axes_str_2D(WCS, 800, 'k', 'r','x','y','')
title('Descripcion de elementos')
axis equal

InF = splitCoord(Xin,Yin);
OutF = splitCoord(Xout,Yout);

CoordenadasF = [InF OutF]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(size(CoordenadasF,1),1);
Y = zeros(size(CoordenadasF,1),1);
Y(1:n,1) = Ynorm(1:n)*P;   % Lugar de aplicacion carga
X(1:n,1) = Xnorm(1:n)*P; 
% X =Xnorm*P;
% Y =Ynorm*P;
CoordenadasF = [CoordenadasF X Y] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = 250;                    %%%%%%%%% DIST A PARTIR DE RADIO INT
cont2=1;
ptsobs = 20;               %%%%%%% NUM PTS OBS
deltaobs = 2*pi/(ptsobs-1);                      
rangeobs = 0:deltaobs:2*pi;
xobs = zeros(size(rangeobs));
yobs = zeros(size(rangeobs));
for i= rangeobs
    xobs(cont2) = (R+tt)*cos(i);
    yobs(cont2) = (R+tt)*sin(i);
    cont2 = cont2 +1;
end

Observacion = [xobs' yobs'];

% figure(2)
subplot(2,2,2)
hold on
plot(CoordenadasF(:,1),CoordenadasF(:,3),'b.')
% plot(CoordenadasF(:,2),CoordenadasF(:,4),'k.')
plot(Observacion(:,1),Observacion(:,2),'mo')
plt_axes_str_2D(WCS, 800, 'k', 'r','x','y','')
title('Puntos de Observacion')
axis equal


%%% IMPORT EXCEL 
Xob = xlsread('IBEM2Doutput.xlsx','','A6:A86');
Yob = xlsread('IBEM2Doutput.xlsx','','B6:B86');
DespX = xlsread('IBEM2Doutput.xlsx','','C6:C86');
DespY = xlsread('IBEM2Doutput.xlsx','','D6:D86');
SigmaXX = xlsread('IBEM2Doutput.xlsx','','E6:E86');
TaoXY = xlsread('IBEM2Doutput.xlsx','','F6:F86');
TaoYX = xlsread('IBEM2Doutput.xlsx','','G6:G86');
SigmaYY = xlsread('IBEM2Doutput.xlsx','','H6:H86');
Resul = [Xob, Yob, DespX, DespY, SigmaXX, TaoXY, TaoYX, SigmaYY];


% figure(3)
subplot(2,2,3)
hold on
plot(Xob,Yob,'b','LineWidth',1.2)
plot(Xob+DespX,Yob+DespY,'r')
quiver(Xob,Yob,DespX,DespY,'b','AutoScaleFactor',0.1)
plt_axes_str_2D(WCS, 800, 'k', 'r','x','y','')
title('Desplazamientos resultantes')
axis equal

error = (abs(max(SigmaXX)-Sigma1)/Sigma1)*100;
sprintf('porcentaje de error de esfuerzo es %d', error)

despSolid = 184.513;
errordesp = (abs(max(DespY)-despSolid)/despSolid)*100;
sprintf('porcentaje de error de esfuerzo es %d', errordesp)


% figure(4)
subplot(2,2,4)
hold on
plot(Xob,Yob,'b','LineWidth',1.2)
plot3(Xob,Yob,SigmaXX,'-ro')
plot3(Xob,Yob,SigmaYY,'-bo')
plot3(Xob,Yob,TaoXY,'-go')
plt_axes_str_2D(WCS, 800, 'k', 'r','x','y','')
title('SigmaXX and SigmaYY')
legend('','XX','YY','TXY')
axis equal

% subplot(1,2,2)
% hold on
% plot(Xob,Yob,'b','LineWidth',1.2)
% plot3(Xob,Yob,SigmaYY,'ro')
% plt_axes_str_2D(WCS, 800, 'k', 'r','x','y','')
% title('SigmaXX')
% axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Empty = strings(size(CoordenadasF,1),1); 

Index = transpose(1:size(CoordenadasF,1)); 
CoordsyCargas = sprintfc('%0.4f',CoordenadasF);
CoordsyCargas(:,7) = num2cell(Index);
CoordsyCargas(:,8) = num2cell(Empty);
Nele1 = 0;
Nele2 = size(CoordenadasF,1);
Nele3 = 0;
Nele4 = 0; 
Nele5 = 0; 
Nele6 = size(Observacion,1);

TitleCell = {'TANQUE', '', '', '', '', '', '', ''};
SpacesCell = {'1', '', '', '', '', '', '', ''};
ELEMENTOSCell = {'ELEMENTOS', '', '', '', '', '', '', ''};
OBSERVACIONCell = {'OBSERV.', '', '', '', '', '', '', ''};
CellMaterial = [sprintfc('%0.2f',G), string(Nu), string(Nele1), string(Nele2), string(Nele3), string(Nele4), string(Nele5), string(Nele6)];

Empty2 = strings(size(Observacion,1),6); 
ObservacionCell = sprintfc('%0.4f',Observacion);
ObservacionCell(:,3:8) = num2cell(Empty2);

TFinal = table([TitleCell;SpacesCell;CellMaterial;ELEMENTOSCell;CoordsyCargas;OBSERVACIONCell;ObservacionCell]);
delete('IBEM2D.xlsx')
filename = 'IBEM2D.xlsx';
writetable(TFinal,filename,'Sheet',1,'Range','A1','WriteVariableNames',false)
