function []=plt_axes_str_2D(cords, L, A, B,label_x,label_y,label_o)

% Px=cords(1:3,1)';
% Py=cords(1:3,2)';
% Pz=cords(1:3,3)';
% S=cords(1:3,4);

Px=cords(1:2,1)'; 
Py=cords(1:2,2)';
S=cords(1:2,3);

cg=S;

label_x = cat(2,'  \bf\it',label_x);
label_y = cat(2,'  \bf\it',label_y);
label_o = cat(2,'  \bf\it',label_o);

plot([cg(1) cg(1)+L*Px(1)],...   %Ploting X axis
      [cg(2) cg(2)+L*Px(2)], A, 'LineWidth', 1.5)
  hold on
plot([cg(1) cg(1)+L*Py(1)],...   %Ploting Y axis
      [cg(2) cg(2)+L*Py(2)], B, 'LineWidth', 1.5)
  hold on
%--------------------------------------------------------------------------
% Modificate Matrix Householder transformation
text(cg(1)+L*Px(1),cg(2)+L*Px(2),label_x)
text(cg(1)+L*Py(1),cg(2)+L*Py(2),label_y)
text(cg(1),cg(2),label_o)
