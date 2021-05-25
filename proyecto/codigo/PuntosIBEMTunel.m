% Juan Camilo Gutierrez 
% 201710009014

clc
clear
close all

% Este programa calcula el input para el programa IBEM del profesor Juan
% Diego Jaramillo, en el caso de una viga prismatica en voladizo con una
% carga en el extremo izquierdo y empotrada en el extremo derecho

%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
% Presión aplicada N/mm^2
P = 100;

% Radio de tunel mm  (tunel de oriente 12,5m)
R = 15000; 

% Espesor de puente cilindrico mm
t = 3000;

% Numero de nodos por lado ingresado
n = 10; 

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
for i= 0:pi/n:pi
    i2 = pi - i;
    Xin(cont) = R*cos(i); 
    Yin(cont) = R*sin(i);  
    Xout(cont) = (R+t)*cos(i2);
    Yout(cont) = (R+t)*sin(i2); 
    
    Xnorm(cont) = -cos(i+pi/n/2); %deriv cos = -sin -> sin (derivo para tangente, y alterno cambiando uno de signo)
    Ynorm(cont) = -sin(i+pi/n/2); %deriv sin = cos -> cos
    Pxnorm(cont) = ( (R+t)*cos(i) + (R+t)*cos(i+pi/n) ) /2;
    Pynorm(cont) = ( (R+t)*sin(i) + (R+t)*sin(i+pi/n) ) /2;
    
    cont = cont + 1;   
end

Xbleft = -R:-t/n:-R-t;
Xbright = R+t:-t/n:R;
Ybleft = zeros(size(Xbleft));
Ybright = zeros(size(Xbright));

figure(1)
subplot(2,1,1)
hold on
quiver(Pxnorm(1:n),Pynorm(1:n),Xnorm(1:n),Ynorm(1:n),'r','AutoScaleFactor',0.2,'LineWidth',1.2);
plot(Xbright,Ybright, 'k','LineWidth',2)
plot(Xbleft,Ybleft, 'k','LineWidth',2)
plot(Xin,Yin, 'k','LineWidth',2)
plot(Xout,Yout, 'b','LineWidth',1.2)
plt_axes_str_2D(WCS, 5000, 'k', 'r','x','y','')
title('Descripcion de elementos')
axis equal

BleftF=splitCoord(Xbleft,Ybleft);
BrightF=splitCoord(Xbright,Ybright);
InF = splitCoord(Xin,Yin);
OutF = splitCoord(Xout,Yout);

CoordenadasF = [BleftF BrightF InF OutF]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(size(CoordenadasF,1),1);
Y = zeros(size(CoordenadasF,1),1);
Y(3*n+1:4*n,1) = Ynorm(1:n)*P;   % Lugar de aplicacion carga
X(3*n+1:4*n,1) = Xnorm(1:n)*P; 
CoordenadasF = [CoordenadasF X Y] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tt = 2500;                   %%%%%%%%% DIST A PARTIR DE RADIO INT
cont2=1;
ptsobs = 81;                 %%%%%%% NUM PTS OBS
deltaobs = pi/(ptsobs-1);                       
rangeobs = 0:deltaobs:pi;
xobs = zeros(size(rangeobs));
yobs = zeros(size(rangeobs));
for i= rangeobs
    xobs(cont2) = (R+tt)*cos(i);
    yobs(cont2) = (R+tt)*sin(i);
    cont2 = cont2 +1;
end

Observacion = [xobs' yobs'];

figure(1)
subplot(2,1,2)
hold on
plot(CoordenadasF(:,1),CoordenadasF(:,3),'b.')
% plot(CoordenadasF(:,2),CoordenadasF(:,4),'k.')
plot(Observacion(:,1),Observacion(:,2),'mo')
plt_axes_str_2D(WCS, 5000, 'k', 'r','x','y','')
title('Puntos de Observacion')
axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Empty = strings(size(CoordenadasF,1),1); 

Index = transpose(1:size(CoordenadasF,1)); 
CoordsyCargas = sprintfc('%0.3f',CoordenadasF);
CoordsyCargas(:,7) = num2cell(Index);
CoordsyCargas(:,8) = num2cell(Empty);
Nele1 = 2*n;
Nele2 = 2*n;
Nele3 = 0;
Nele4 = 0; 
Nele5 = 0; 
Nele6 = size(Observacion,1);

TitleCell = {'TUNEL', '', '', '', '', '', '', ''};
SpacesCell = {'1', '', '', '', '', '', '', ''};
ELEMENTOSCell = {'ELEMENTOS', '', '', '', '', '', '', ''};
OBSERVACIONCell = {'OBSERV.', '', '', '', '', '', '', ''};
CellMaterial = [sprintfc('%0.2f',G), string(Nu), string(Nele1), string(Nele2), string(Nele3), string(Nele4), string(Nele5), string(Nele6)];

Empty2 = strings(size(Observacion,1),6); 
ObservacionCell = sprintfc('%0.3f',Observacion);
ObservacionCell(:,3:8) = num2cell(Empty2);

TFinal = table([TitleCell;SpacesCell;CellMaterial;ELEMENTOSCell;CoordsyCargas;OBSERVACIONCell;ObservacionCell]);
delete('IBEM2D.xlsx')
filename = 'IBEM2D.xlsx';
writetable(TFinal,filename,'Sheet',1,'Range','A1','WriteVariableNames',false)
