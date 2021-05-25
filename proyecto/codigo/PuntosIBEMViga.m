% Juan Camilo Gutierrez 
% 201710009014

clc
clear
close all

% Este programa calcula el input para el programa IBEM del profesor Juan
% Diego Jaramillo, en el caso de una viga prismatica en voladizo con una
% carga en el extremo izquierdo y empotrada en el extremo derecho

%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%
% Carga aplicada en extremo izquierdo de viga
P = 100;

% Longitud m de la viga
L = 200; 

% Altura m de la viga
h = 20;

% Base de la viga
b = 1;

%Punto de evalucación de desplazamiento vertical
X = 0;

% Numero de nodos por lado ingresado
n = 1200; 

E = 72400;
Nu = 0.33;
G = 28000.00;
I = (b*h^3)/(12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_A = -P/(1*h); 

vmax = -(P*L^3)/(3*E*I);
disp('Desplazamiento Y maximo es igual a ')
disp(vmax)

vx = @(x) (P/(6*E*I))*(-x^3+3*L^2*x-2*L^3);
sprintf('Desplazamiento en X igual a %d', X)
disp(vx(X))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


WCS = eye(3);

verts = [L 1/2*h 1
         L -1/2*h 1
         0 -1/2*h 1
         0 1/2*h 1]';
loops = 4;

% Amplio el numero de puntos de muestra y armo las matrices de puntos
Resampled_P = resample_polyg(verts, n);
CoordenadasF = zeros(4,size(Resampled_P,2)-1);

figure(1)
subplot(2,1,1)
hold on
draw_polyg_faces_2D(Resampled_P, loops*n, 'b')
plt_axes_str_2D(WCS, 20, 'k', 'r','x','y','o')
axis([-L*0.3 L -h*1.8 h*1.8])
title('Descripcion de elementos')
axis equal

for i = 1:size(Resampled_P,2)
    if i ~= size(Resampled_P,2)
        Xi=Resampled_P(1,i);
        Xf=Resampled_P(1,i+1);
        Yi=Resampled_P(2,i);
        Yf=Resampled_P(2,i+1);
    else
        Xi=Resampled_P(1,i);
        Xf=Resampled_P(1,1);
        Yi=Resampled_P(2,i);
        Yf=Resampled_P(2,1);        
    end
    coord= [Xi Xf Yi Yf]';
    CoordenadasF(:,i) = coord; 
end

X = zeros(size(CoordenadasF,2),1);
Y = zeros(size(CoordenadasF,2),1);
Y(2*n+1:3*n,1) = ones(n,1)*F_A;   % Lugar de aplicacion carga

CoordenadasF = [CoordenadasF' X Y] ;

xobs = 0:4:L;
yobs = h/2:-h/(size(xobs,2)-1):-h/2;
Observacion = [xobs' ones(size(xobs,2),1)*yobs(2)
               xobs' ones(size(xobs,2),1)*yobs(size(yobs,2)-1)
               xobs' ones(size(xobs,2),1)*yobs(ceil(size(yobs,2)/2))]; % longitud superior viga
disp(size(Observacion,1))

figure(1)
subplot(2,1,2)
hold on
plot(CoordenadasF(:,1),CoordenadasF(:,3),'b.')
% plot(CoordenadasF(:,2),CoordenadasF(:,4),'k.')
plot(Observacion(:,1),Observacion(:,2),'mo')
plt_axes_str_2D(WCS, 20, 'k', 'r','x','y','o')
axis([-L*0.3 L -h*1.8 h*1.8])
axis equal
title('Puntos de Observacion')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Empty = strings(size(CoordenadasF,1),1); 
% CoordenadasFstr = [CoordenadasF Empty] ;

Index = transpose(1:size(CoordenadasF,1)); 
CoordsyCargas = sprintfc('%0.4f',CoordenadasF);
CoordsyCargas(:,7) = num2cell(Index);
CoordsyCargas(:,8) = num2cell(Empty);
Nele1 = 1*n;
Nele2 = 3*n;
Nele3 = 0;
Nele4 = 0; 
Nele5 = 0; 
Nele6 = size(Observacion,1);

TitleCell = {'VIGA EN VOLADIZO', '', '', '', '', '', '', ''};
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

file = 'IBEM2doutput.csv';