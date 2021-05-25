function CoordenadasF=splitCoord(X,Y)

CoordenadasF = zeros(4,size(X,2)-1);

for i = 1:size(X,2)-1
    Xi=X(i);
    Xf=X(i+1);
    Yi=Y(i);
    Yf=Y(i+1);
    coord= [Xi Xf Yi Yf]';
    CoordenadasF(:,i) = coord;
end

end