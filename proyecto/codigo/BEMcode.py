from IPython import get_ipython;   
get_ipython().magic('reset -sf')

import numpy as np
import math as m

# pinches sumas

def Gij(ziX,ziY,Xi,Xf,Yi,Yf,G,Nu):  # (xc,yc,Xi,Xf,Yi,Yf,G,Nu,G11,G12,G21,G22)   xc, yc = coordenadas de x

    xX = (Xi+Xf)/2 # coord Zi(x) nodo
    xY = (Yi+Yf)/2 # coord Zi(y) nodo
    rx = (Xf-Xi)
    ry = (Yf-Yi)    
    ele= m.sqrt(rx**2 + ry**2) #Tamaño nodo l
    nx=-ry/ele # normal del elemento
    ny=rx/ele 

    C1 = 1/(8*m.pi*G*(1-Nu))
    C2 = (3-4*Nu)

    if (m.sqrt((xX-ziX)**2+(xY-ziY)**2) < ele/10): # Si el vector X es igual a Zi
        G11 = C1*ele*(ny**2+C2*(1-m.log(ele/2)))
        G22 = C1*ele*(nx**2+C2*(1-m.log(ele/2)))
        G12 = -C1*(ele*nx*ny)
        G21 = G12

    else: # Si el vector Zi está en otra posicion a X
        nx=rx/ele #direccion del elemento
        ny=ry/ele 
        # Cuadratura de Gauss
        Pt1 = -m.sqrt(3/5)
        Pt2 = 0
        Pt3 = m.sqrt(3/5)

        W1 = 5/9
        W2 = 8/9
        W3 = 5/9

        ## coordenadas puntos de Gauss
        x1=xX+nx*Pt1*ele/2
        y1=xY+ny*Pt1*ele/2

        x2=xX+nx*Pt2*ele/2
        y2=xY+ny*Pt2*ele/2

        x3=xX+nx*Pt3*ele/2
        y3=xY+ny*Pt3*ele/2

        ## Calculo radio, a partir de puntos de cuadratura de Gauss y Zi  r = x - zi
        rx1=x1-ziX
        ry1=y1-ziY
        r1=m.sqrt(rx1**2+ry1**2)

        rx2=x2-ziX
        ry2=y2-ziY
        r2=m.sqrt(rx2**2+ry2**2)

        rx3=x3-ziX
        ry3=y3-ziY
        r3=m.sqrt(rx3**2+ry3**2)

        G11p1=C1*(C2*m.log(1/r1)+(rx1*rx1)/r1**2)
        G22p1=C1*(C2*m.log(1/r1)+(ry1*ry1)/r1**2)
        G12p1=C1*((rx1*ry1)/r1**2)
        G21p1=G12p1       

        G11p2=C1*(C2*m.log(1/r2)+(rx2*rx2)/r2**2)
        G22p2=C1*(C2*m.log(1/r2)+(ry2*ry2)/r2**2)
        G12p2=C1*((rx2*ry2)/r2**2)
        G21p2=G12p2

        G11p3=C1*(C2*m.log(1/r3)+(rx3*rx3)/r3**2)
        G22p3=C1*(C2*m.log(1/r3)+(ry3*ry3)/r3**2)
        G12p3=C1*(rx3*ry3/r3/r3)
        G21p3=G12p3

        G11=(G11p1*W1+G11p2*W2+G11p3*W3)*ele/2
        G22=(G22p1*W1+G22p2*W2+G22p3*W3)*ele/2
        G12=(G12p1*W1+G12p2*W2+G12p3*W3)*ele/2
        G21=(G21p1*W1+G21p2*W2+G21p3*W3)*ele/2       

    return G11,G22,G12,G21


def Hij(ziX,ziY,Xi,Xf,Yi,Yf,n1,n2,G,Nu):    # (xc,yc,Xi,Xf,Yi,Yf,n1,n2,G,Nu,H11,H12,H21,H22,Iob)

    xX = (Xi+Xf)/2 # coord Zi(x) nodo
    xY = (Yi+Yf)/2 # coord Zi(y) nodo
    rx = (Xf-Xi)
    ry = (Yf-Yi)    
    ele= m.sqrt(rx**2 + ry**2) #Tamaño nodo l
    nx=rx/ele # direccion del elemento en x
    ny=ry/ele # direccion del elemento en y
    
    C1 = -1/(4*m.pi*(1-Nu))
    C2 = 1-2*Nu
    
    if (m.sqrt((xX-ziX)**2+(xY-ziY)**2) <= ele/10): # Si el vector X es igual a Zi
        H11 = -1/2
        H22 = -1/2
        H12 = 0
        H21 = 0
        
    else:
        # Cuadratura de Gauss
        Pt1 = -m.sqrt(3/5)
        Pt2 = 0
        Pt3 = m.sqrt(3/5)

        W1 = 5/9
        W2 = 8/9
        W3 = 5/9

        ## coordenadas puntos de Gauss en xi, debe ser en x
        x1=xX+nx*Pt1*ele/2
        y1=xY+ny*Pt1*ele/2

        x2=xX+nx*Pt2*ele/2
        y2=xY+ny*Pt2*ele/2

        x3=xX+nx*Pt3*ele/2
        y3=xY+ny*Pt3*ele/2
        
        ## Calculo radios, a partir de puntos de cuadratura de Gauss y Zi  r = x - zi
        rx1=x1-ziX
        ry1=y1-ziY
        r1=m.sqrt(rx1**2+ry1**2)

        rx2=x2-ziX
        ry2=y2-ziY
        r2=m.sqrt(rx2**2+ry2**2)

        rx3=x3-ziX
        ry3=y3-ziY
        r3=m.sqrt(rx3**2+ry3**2)

        rn1=rx1*n1+ry1*n2
        rn2=rx2*n1+ry2*n2
        rn3=rx3*n1+ry3*n2

        H11p1=(C1/(r1**2))*(rn1*(C2+2*rx1*rx1/r1**2))
        H22p1=(C1/(r1**2))*(rn1*(C2+2*ry1*ry1/r1**2))
        H12p1=(C1/(r1**2))*(rn1*(2*rx1*ry1/r1**2)+C2*(rx1*n2-ry1*n1))
        H21p1=(C1/(r1**2))*(rn1*(2*ry1*rx1/r1**2)+C2*(ry1*n1-rx1*n2))
        
        H11p2=(C1/(r2**2))*(rn2*(C2+2*rx2*rx2/r2**2))
        H22p2=(C1/(r2**2))*(rn2*(C2+2*ry2*ry2/r2**2))
        H12p2=(C1/(r2**2))*(rn2*(2*rx2*ry2/r2**2)+C2*(rx2*n2-ry2*n1))
        H21p2=(C1/(r2**2))*(rn2*(2*ry2*rx2/r2**2)+C2*(ry2*n1-rx2*n2))
        
        H11p3=(C1/(r3**2))*(rn3*(C2+2*rx3*rx3/r3**2))
        H22p3=(C1/(r3**2))*(rn3*(C2+2*ry3*ry3/r3**2))
        H12p3=(C1/(r3**2))*(rn3*(2*rx3*ry3/r3**2)+C2*(rx3*n2-ry3*n1))
        H21p3=(C1/(r3**2))*(rn3*(2*ry3*rx3/r3**2)+C2*(ry3*n1-rx3*n2))


        H11=(H11p1*W1+H11p2*W2+H11p3*W3)*ele/2
        H22=(H22p1*W1+H22p2*W2+H22p3*W3)*ele/2
        H12=(H12p1*W1+H12p2*W2+H12p3*W3)*ele/2
        H21=(H21p1*W1+H21p2*W2+H21p3*W3)*ele/2

    return H11,H22,H12,H21

archivo = open('tank3.csv','r')
texto = archivo.read()
lineas = texto.split('\n')

title = lineas[0].split(';')[0]
G = float(lineas[2].split(';')[0])
Nu = float(lineas[2].split(';')[1])
Nele1 = int(lineas[2].split(';')[2])
Nele2 = int(lineas[2].split(';')[3])
Nele3 = int(lineas[2].split(';')[4])
Nele4 = int(lineas[2].split(';')[5])
PtsObs = int(lineas[2].split(';')[7])

NeleT = Nele1+Nele2

Elems = np.zeros([NeleT,7])
for i in range(NeleT):
    for j in range(7):
        Elems[i][j] = lineas[i+4].split(';')[j]
        
Obs = np.zeros([PtsObs,2])
for i in range(PtsObs):
    for j in range(2):
        Obs[i][j] = lineas[i+NeleT+5].split(';')[j]

eleT = len(Elems)
K = np.zeros([eleT*2,eleT*2])
V = np.zeros([eleT*2,1])

zaux = 0
for z in range(0,eleT): ########### SISTEMA DE ECUACIONES
    ziX = (Elems[z][0]+Elems[z][1])/2
    ziY = (Elems[z][2]+Elems[z][3])/2
    xaux = 0
    for x in range(0,eleT):
        ele = m.sqrt((Elems[x][1]-Elems[x][0])**2 + (Elems[x][3]-Elems[x][2])**2)

        t1 = Elems[x][4]
        t2 = Elems[x][5]
        
        nx = (Elems[x][2]-Elems[x][3])/ele
        ny = (Elems[x][1]-Elems[x][0])/ele
        
        Hij_Zi1_X1 = Hij(ziX,ziY,Elems[x][0],Elems[x][1],Elems[x][2],Elems[x][3],nx,ny,G,Nu)
        H11 = Hij_Zi1_X1[0]
        H22 = Hij_Zi1_X1[1]
        H12 = Hij_Zi1_X1[2]
        H21 = Hij_Zi1_X1[3]

        Gij_Zi1_X1 = Gij(ziX,ziY,Elems[x][0],Elems[x][1],Elems[x][2],Elems[x][3],G,Nu)
        G11 = Gij_Zi1_X1[0]
        G22 = Gij_Zi1_X1[1]
        G12 = Gij_Zi1_X1[2]
        G21 = Gij_Zi1_X1[3]
        
        V[zaux] += G11*t1 + G21*t2 
        V[zaux+1] +=  G12*t1 + G22*t2    
        
        if x == z:
            K[zaux][xaux] = H11+1
            K[zaux][xaux+1] = H21
            K[zaux+1][xaux] = H12
            K[zaux+1][xaux+1] = H22+1

        else:
            K[zaux][xaux] = H11
            K[zaux][xaux+1] = H21
            K[zaux+1][xaux] = H12
            K[zaux+1][xaux+1] = H22
       
        xaux += 2
    zaux += 2
u = np.linalg.solve(K,V)       

############# DESPLAZAMIENTO PTS DE OBSERVACION
ObsF = np.zeros([len(Obs)*2,1])
zaux = 0
for zobs in range(len(Obs)):
    obsX = Obs[zobs][0]
    obsY = Obs[zobs][1]
    xaux = 0
    for x in range(0,eleT):   
        ele = m.sqrt((Elems[x][1]-Elems[x][0])**2 + (Elems[x][3]-Elems[x][2])**2)
        
        t1 = Elems[x][4]
        t2 = Elems[x][5]
        
        u1 = u[xaux][0]
        u2 = u[xaux+1][0]

        nx = (Elems[x][2]-Elems[x][3])/ele
        ny = (Elems[x][1]-Elems[x][0])/ele
        
        Hij_Zi1_X1 = Hij(obsX,obsY,Elems[x][0],Elems[x][1],Elems[x][2],Elems[x][3],nx,ny,G,Nu)
        H11 = Hij_Zi1_X1[0]
        H22 = Hij_Zi1_X1[1]
        H12 = Hij_Zi1_X1[2]
        H21 = Hij_Zi1_X1[3]
    
        Gij_Zi1_X1 = Gij(obsX,obsY,Elems[x][0],Elems[x][1],Elems[x][2],Elems[x][3],G,Nu)
        G11 = Gij_Zi1_X1[0]
        G22 = Gij_Zi1_X1[1]
        G12 = Gij_Zi1_X1[2]
        G21 = Gij_Zi1_X1[3]
        
        ObsF[zaux] += G11*t1 + G21*t2 -H11*u1 -H21*u2
        ObsF[zaux+1] +=  G12*t1 + G22*t2 -H12*u1 -H22*u2
        
        xaux += 2
    zaux += 2



alpha = 1
beta = 2
phi = 4

def Gkij(ziX,ziY,Xi,Xf,Yi,Yf,G,Nu):
    xX = (Xi+Xf)/2 # coord Zi(x) nodo
    xY = (Yi+Yf)/2 # coord Zi(y) nodo
    rx = (Xf-Xi)
    ry = (Yf-Yi)    
    ele= m.sqrt(rx**2 + ry**2) #Tamaño nodo l
    
    nx=rx/ele #direccion del elemento
    ny=ry/ele 
    # Cuadratura de Gauss
    Pt1 = -m.sqrt(3/5)
    Pt2 = 0
    Pt3 = m.sqrt(3/5)

    W1 = 5/9
    W2 = 8/9
    W3 = 5/9

    ## coordenadas puntos de Gauss
    x1=xX+nx*Pt1*ele/2
    y1=xY+ny*Pt1*ele/2

    x2=xX+nx*Pt2*ele/2
    y2=xY+ny*Pt2*ele/2

    x3=xX+nx*Pt3*ele/2
    y3=xY+ny*Pt3*ele/2

    ## Calculo radio, a partir de puntos de cuadratura de Gauss y Zi  r = x - zi
    rx1=x1-ziX
    ry1=y1-ziY
    r1=m.sqrt(rx1**2+ry1**2)
    rx_1 = -rx1/r1
    ry_1 = -ry1/r1

    rx2=x2-ziX
    ry2=y2-ziY
    r2=m.sqrt(rx2**2+ry2**2)
    rx_2 = -rx2/r2
    ry_2 = -ry2/r2

    rx3=x3-ziX
    ry3=y3-ziY
    r3=m.sqrt(rx3**2+ry3**2)    
    rx_3 = -rx3/r3
    ry_3 = -ry3/r3
    
    G111p1 = 1/(4*m.pi*alpha*(1-Nu))*(1/r1**alpha)*((1-2*Nu)*(rx_1 + rx_1 - rx_1) + beta*rx_1*rx_1*rx_1)
    G112p1 = 1/(4*m.pi*alpha*(1-Nu))*(1/r1**alpha)*((1-2*Nu)*(ry_1 + 0 - 0) + beta*rx_1*rx_1*ry_1)
    
    G121p1 = 1/(4*m.pi*alpha*(1-Nu))*(1/r1**alpha)*((1-2*Nu)*(0 + rx_1 - 0) + beta*rx_1*ry_1*rx_1)
    G122p1 = 1/(4*m.pi*alpha*(1-Nu))*(1/r1**alpha)*((1-2*Nu)*(0 + 0 - rx_1) + beta*rx_1*ry_1*ry_1)
    
    G211p1 = 1/(4*m.pi*alpha*(1-Nu))*(1/r1**alpha)*((1-2*Nu)*(0 + 0 - ry_1) + beta*ry_1*rx_1*rx_1)
    G212p1 = 1/(4*m.pi*alpha*(1-Nu))*(1/r1**alpha)*((1-2*Nu)*(0 + ry_1 - 0) + beta*ry_1*rx_1*ry_1)
    
    G221p1 = 1/(4*m.pi*alpha*(1-Nu))*(1/r1**alpha)*((1-2*Nu)*(rx_1 + 0 - 0) + beta*ry_1*ry_1*rx_1)
    G222p1 = 1/(4*m.pi*alpha*(1-Nu))*(1/r1**alpha)*((1-2*Nu)*(ry_1 + ry_1 - ry_1) + beta*ry_1*ry_1*ry_1)
    #######################  
  
    G111p2 = 1/(4*m.pi*alpha*(1-Nu))*(1/r2**alpha)*((1-2*Nu)*(rx_2 + rx_2 - rx_2) + beta*rx_2*rx_2*rx_2)
    G112p2 = 1/(4*m.pi*alpha*(1-Nu))*(1/r2**alpha)*((1-2*Nu)*(ry_2 + 0 - 0) + beta*rx_2*rx_2*ry_2)
    
    G121p2 = 1/(4*m.pi*alpha*(1-Nu))*(1/r2**alpha)*((1-2*Nu)*(0 + rx_2 - 0) + beta*rx_2*ry_2*rx_2)
    G122p2 = 1/(4*m.pi*alpha*(1-Nu))*(1/r2**alpha)*((1-2*Nu)*(0 + 0 - rx_2) + beta*rx_2*ry_2*ry_2)
    
    G211p2 = 1/(4*m.pi*alpha*(1-Nu))*(1/r2**alpha)*((1-2*Nu)*(0 + 0 - ry_2) + beta*ry_2*rx_2*rx_2)
    G212p2 = 1/(4*m.pi*alpha*(1-Nu))*(1/r2**alpha)*((1-2*Nu)*(0 + ry_2 - 0) + beta*ry_2*rx_2*ry_2)
    
    G221p2 = 1/(4*m.pi*alpha*(1-Nu))*(1/r2**alpha)*((1-2*Nu)*(rx_2 + 0 - 0) + beta*ry_2*ry_2*rx_2)
    G222p2 = 1/(4*m.pi*alpha*(1-Nu))*(1/r2**alpha)*((1-2*Nu)*(ry_2 + ry_2 - ry_2) + beta*ry_2*ry_2*ry_2) 
    ######################
    
    G111p3 = 1/(4*m.pi*alpha*(1-Nu))*(1/r3**alpha)*((1-2*Nu)*(rx_3 + rx_3 - rx_3) + beta*rx_3*rx_3*rx_3)
    G112p3 = 1/(4*m.pi*alpha*(1-Nu))*(1/r3**alpha)*((1-2*Nu)*(ry_3 + 0 - 0) + beta*rx_3*rx_3*ry_3)
    
    G121p3 = 1/(4*m.pi*alpha*(1-Nu))*(1/r3**alpha)*((1-2*Nu)*(0 + rx_3 - 0) + beta*rx_3*ry_3*rx_3)
    G122p3 = 1/(4*m.pi*alpha*(1-Nu))*(1/r3**alpha)*((1-2*Nu)*(0 + 0 - rx_3) + beta*rx_3*ry_3*ry_3)
    
    G211p3 = 1/(4*m.pi*alpha*(1-Nu))*(1/r3**alpha)*((1-2*Nu)*(0 + 0 - ry_3) + beta*ry_3*rx_3*rx_3)
    G212p3 = 1/(4*m.pi*alpha*(1-Nu))*(1/r3**alpha)*((1-2*Nu)*(0 + ry_3 - 0) + beta*ry_3*rx_3*ry_3)
    
    G221p3 = 1/(4*m.pi*alpha*(1-Nu))*(1/r3**alpha)*((1-2*Nu)*(rx_3 + 0 - 0) + beta*ry_3*ry_3*rx_3)
    G222p3 = 1/(4*m.pi*alpha*(1-Nu))*(1/r3**alpha)*((1-2*Nu)*(ry_3 + ry_3 - ry_3) + beta*ry_3*ry_3*ry_3)
    
    G111 = (G111p1*W1 + G111p2*W2 + G111p3*W3) * ele/2
    G112 = (G112p1*W1 + G112p2*W2 + G112p3*W3) * ele/2
    
    G121 = (G121p1*W1 + G121p2*W2 + G121p3*W3) * ele/2
    G122 = (G122p1*W1 + G122p2*W2 + G122p3*W3) * ele/2
    
    G211 = (G211p1*W1 + G211p2*W2 + G211p3*W3) * ele/2
    G212 = (G212p1*W1 + G212p2*W2 + G212p3*W3) * ele/2
    
    G221 = (G221p1*W1 + G221p2*W2 + G221p3*W3) * ele/2
    G222 = (G222p1*W1 + G222p2*W2 + G222p3*W3) * ele/2    
    
    return G111,G112,G121,G122,G211,G212,G221,G222

def Hkij(ziX,ziY,Xi,Xf,Yi,Yf,n1,n2,G,Nu):    # (xc,yc,Xi,Xf,Yi,Yf,n1,n2,G,Nu,H11,H12,H21,H22,Iob)

    xX = (Xi+Xf)/2 # coord Zi(x) nodo
    xY = (Yi+Yf)/2 # coord Zi(y) nodo
    rx = (Xf-Xi)
    ry = (Yf-Yi)    
    ele= m.sqrt(rx**2 + ry**2) #Tamaño nodo l
    nx=rx/ele # direccion del elemento en x
    ny=ry/ele # direccion del elemento en y
    
    C1 = G/(2*m.pi*alpha*(1-Nu))
    C2 = 1-2*Nu
    C3 = 1-4*Nu
    
    # Cuadratura de Gauss
    Pt1 = -m.sqrt(3/5)
    Pt2 = 0
    Pt3 = m.sqrt(3/5)

    W1 = 5/9
    W2 = 8/9
    W3 = 5/9

    ## coordenadas puntos de Gauss en xi, debe ser en x
    x1=xX+nx*Pt1*ele/2
    y1=xY+ny*Pt1*ele/2

    x2=xX+nx*Pt2*ele/2
    y2=xY+ny*Pt2*ele/2

    x3=xX+nx*Pt3*ele/2
    y3=xY+ny*Pt3*ele/2
    
    ## Calculo radios, a partir de puntos de cuadratura de Gauss y Zi  r = x - zi
    rx1=x1-ziX
    ry1=y1-ziY
    r1=m.sqrt(rx1**2+ry1**2)
    drdn1 = (rx1*nx+ry1*ny)/r1
    rx_1 = -rx1/r1
    ry_1 = -ry1/r1
    
    rx2=x2-ziX
    ry2=y2-ziY
    r2=m.sqrt(rx2**2+ry2**2)
    drdn2 = (rx2*nx+ry2*ny)/r2
    rx_2 = -rx2/r2
    ry_2 = -ry2/r2

    rx3=x3-ziX
    ry3=y3-ziY
    r3=m.sqrt(rx3**2+ry3**2)
    drdn3 = (rx3*nx+ry3*ny)/r3
    rx_3 = -rx3/r3
    ry_3 = -ry3/r3

    rn1=rx1*n1+ry1*n2
    rn2=rx2*n1+ry2*n2
    rn3=rx3*n1+ry3*n2

    H111p1= C1*(1/r1**beta)*( beta* drdn1 (C2 * rx_1 + Nu(rx_1 + rx_1) - phi*rx_1*rx_1*rx_1 ) + beta*Nu*(nx*rx_1*rx_1 + nx*rx_1*rx_1) + C2*(beta*nx*rx_1*rx_1 + nx+nx) - C3*nx )
    H112p1= C1*(1/r1**beta)*( beta* drdn1 (C2 * 0 + Nu(rx_1 + 0) - phi*rx_1*ry_1*rx_1 ) + beta*Nu*(nx*rx_1*rx_1 + ny*rx_1*rx_1) + C2*(beta*nx*rx_1*ry_1 + ny+0) - C3*0 )
    
    H121p1= C1*(1/r1**beta)*( beta* drdn1 (C2 * 0 + Nu(0 + ry_1) - phi*ry_1*rx_1*rx_1 ) + beta*Nu*(ny*rx_1*rx_1 + nx*ry_1*rx_1) + C2*(beta*nx*ry_1*rx_1 + 0+ny) - C3*0 )
    H122p1= C1*(1/r1**beta)*( beta* drdn1 (C2 * rx_1 + Nu(0 + 0) - phi*ry_1*ry_1*rx_1 ) + beta*Nu*(ny*ry_1*rx_1 + ny*ry_1*rx_1) + C2*(beta*nx*ry_1*ry_1 + 0+0) - C3*nx )
    
    H211p1= C1*(1/r1**beta)*( beta* drdn1 (C2 * ry_1 + Nu(0 + 0) - phi*rx_1*rx_1*ry_1 ) + beta*Nu*(nx*rx_1*ry_1 + nx*rx_1*ry_1) + C2*(beta*ny*rx_1*rx_1 + 0+0) - C3*ny )
    H212p1= C1*(1/r1**beta)*( beta* drdn1 (C2 * 0 + Nu(0 + rx_1) - phi*rx_1*ry_1*ry_1 ) + beta*Nu*(nx*ry_1*ry_1 + ny*rx_1*ry_1) + C2*(beta*ny*rx_1*ry_1 + 0+nx) - C3*0 )
    
    H221p1= C1*(1/r1**beta)*( beta* drdn1 (C2 * 0 + Nu(ry_1 + 0) - phi*ry_1*rx_1*ry_1 ) + beta*Nu*(ny*ry_1*ry_1 + nx*ry_1*ry_1) + C2*(beta*ny*ry_1*rx_1 + nx+0) - C3*0 )
    H222p1= C1*(1/r1**beta)*( beta* drdn1 (C2 * ry_1 + Nu(ry_1 + ry_1) - phi*ry_1*ry_1*ry_1 ) + beta*Nu*(ny*ry_1*ry_1 + ny*ry_1*ry_1) + C2*(beta*ny*ry_1*ry_1 + ny+ny) - C3*ny )

###############
    H111p2= C1*(1/r2**beta)*( beta* drdn2 (C2 * rx_2 + Nu(rx_2 + rx_2) - phi*rx_2*rx_2*rx_2 ) + beta*Nu*(nx*rx_2*rx_2 + nx*rx_2*rx_2) + C2*(beta*nx*rx_2*rx_2 + nx+nx) - C3*nx )
    H112p2= C1*(1/r2**beta)*( beta* drdn2 (C2 * 0 + Nu(rx_2 + 0) - phi*rx_2*ry_2*rx_2 ) + beta*Nu*(nx*rx_2*rx_2 + ny*rx_2*rx_2) + C2*(beta*nx*rx_2*ry_2 + ny+0) - C3*0 )
    
    H121p2= C1*(1/r2**beta)*( beta* drdn2 (C2 * 0 + Nu(0 + ry_2) - phi*ry_2*rx_2*rx_2 ) + beta*Nu*(ny*rx_2*rx_2 + nx*ry_2*rx_2) + C2*(beta*nx*ry_2*rx_2 + 0+ny) - C3*0 )
    H122p2= C1*(1/r2**beta)*( beta* drdn2 (C2 * rx_2 + Nu(0 + 0) - phi*ry_2*ry_2*rx_2 ) + beta*Nu*(ny*ry_2*rx_2 + ny*ry_2*rx_2) + C2*(beta*nx*ry_2*ry_2 + 0+0) - C3*nx )
    
    H211p2= C1*(1/r2**beta)*( beta* drdn2 (C2 * ry_2 + Nu(0 + 0) - phi*rx_2*rx_2*ry_2 ) + beta*Nu*(nx*rx_2*ry_2 + nx*rx_2*ry_2) + C2*(beta*ny*rx_2*rx_2 + 0+0) - C3*ny )
    H212p2= C1*(1/r2**beta)*( beta* drdn2 (C2 * 0 + Nu(0 + rx_2) - phi*rx_2*ry_2*ry_2 ) + beta*Nu*(nx*ry_2*ry_2 + ny*rx_2*ry_2) + C2*(beta*ny*rx_2*ry_2 + 0+nx) - C3*0 )
    
    H221p2= C1*(1/r2**beta)*( beta* drdn2 (C2 * 0 + Nu(ry_2 + 0) - phi*ry_2*rx_2*ry_2 ) + beta*Nu*(ny*ry_2*ry_2 + nx*ry_2*ry_2) + C2*(beta*ny*ry_2*rx_2 + nx+0) - C3*0 )
    H222p2= C1*(1/r2**beta)*( beta* drdn2 (C2 * ry_2 + Nu(ry_2 + ry_2) - phi*ry_2*ry_2*ry_2 ) + beta*Nu*(ny*ry_2*ry_2 + ny*ry_2*ry_2) + C2*(beta*ny*ry_2*ry_2 + ny+ny) - C3*ny )

##############
    H111p3= C1*(1/r3**beta)*( beta* drdn3 (C2 * rx_3 + Nu(rx_3 + rx_3) - phi*rx_3*rx_3*rx_3 ) + beta*Nu*(nx*rx_3*rx_3 + nx*rx_3*rx_3) + C2*(beta*nx*rx_3*rx_3 + nx+nx) - C3*nx )
    H112p3= C1*(1/r3**beta)*( beta* drdn3 (C2 * 0 + Nu(rx_3 + 0) - phi*rx_3*ry_3*rx_3 ) + beta*Nu*(nx*rx_3*rx_3 + ny*rx_3*rx_3) + C2*(beta*nx*rx_3*ry_3 + ny+0) - C3*0 )
    
    H121p3= C1*(1/r3**beta)*( beta* drdn3 (C2 * 0 + Nu(0 + ry_3) - phi*ry_3*rx_3*rx_3 ) + beta*Nu*(ny*rx_3*rx_3 + nx*ry_3*rx_3) + C2*(beta*nx*ry_3*rx_3 + 0+ny) - C3*0 )
    H122p3= C1*(1/r3**beta)*( beta* drdn3 (C2 * rx_3 + Nu(0 + 0) - phi*ry_3*ry_3*rx_3 ) + beta*Nu*(ny*ry_3*rx_3 + ny*ry_3*rx_3) + C2*(beta*nx*ry_3*ry_3 + 0+0) - C3*nx )
    
    H211p3= C1*(1/r3**beta)*( beta* drdn3 (C2 * ry_3 + Nu(0 + 0) - phi*rx_3*rx_3*ry_3 ) + beta*Nu*(nx*rx_3*ry_3 + nx*rx_3*ry_3) + C2*(beta*ny*rx_3*rx_3 + 0+0) - C3*ny )
    H212p3= C1*(1/r3**beta)*( beta* drdn3 (C2 * 0 + Nu(0 + rx_3) - phi*rx_3*ry_3*ry_3 ) + beta*Nu*(nx*ry_3*ry_3 + ny*rx_3*ry_3) + C2*(beta*ny*rx_3*ry_1 + 0+nx) - C3*0 )
    
    H221p3= C1*(1/r3**beta)*( beta* drdn3 (C2 * 0 + Nu(ry_3 + 0) - phi*ry_3*rx_3*ry_3 ) + beta*Nu*(ny*ry_3*ry_3 + nx*ry_3*ry_3) + C2*(beta*ny*ry_3*rx_3 + nx+0) - C3*0 )
    H222p3= C1*(1/r3**beta)*( beta* drdn3 (C2 * ry_3 + Nu(ry_3 + ry_3) - phi*ry_3*ry_3*ry_3 ) + beta*Nu*(ny*ry_3*ry_3 + ny*ry_3*ry_3) + C2*(beta*ny*ry_3*ry_3 + ny+ny) - C3*ny )
    
    H111=(H111p1*W1+H111p2*W2+H111p3*W3)*ele/2
    H112=(H112p1*W1+H112p2*W2+H112p3*W3)*ele/2

    H121=(H121p1*W1+H121p2*W2+H121p3*W3)*ele/2
    H122=(H122p1*W1+H122p2*W2+H122p3*W3)*ele/2
    
    H211=(H211p1*W1+H211p2*W2+H211p3*W3)*ele/2
    H212=(H212p1*W1+H212p2*W2+H212p3*W3)*ele/2

    H221=(H221p1*W1+H221p2*W2+H221p3*W3)*ele/2
    H222=(H222p1*W1+H222p2*W2+H222p3*W3)*ele/2
    
    return H111,H112,H121,H122,H211,H212,H221,H222   
    
# G[i,j] +=  1/(8*pi*g*(1-mu))*((3-4*mu)*(i==j)*np.log(1/r) + (R[i,0]*R[j,0])/r**2)*w
# H[i,j] += -1/(4*pi*(1-mu)*r**2)*(R.T.dot(N)*((1-2*mu)*(i==j) + (2*R[i,0]*R[j,0])/r**2) + (1-2*mu)*(R[i,0]*N[j,0] - R[j,0]*N[i,0]))*w
                        
# for k in (0,1):
#     C1 = 1/(4*pi*(1-mu)*r)
#     ri = (-R[i,0]/r)
#     rj = (-R[j,0]/r)
#     rk = (-R[k,0]/r)
    
#     kirj = (k == i) * rj
#     kjri = (k == j) * ri
#     ijrk = (i == j) * rk
    
#     Gkij[k,i,j] +=  C1*((1-2*mu)*(kirj + kjri - ijrk) + 2*(ri)*(rj)*(rk))*w
    
#     C2 = g/(2*pi*(1-mu)*r**2)
#     drdn = -(np.dot(R.T,N))/r
#     A1 = (1 - 2*mu) * ijrk
#     A2 = mu * (kirj + kjri)
#     A3 = 4*ri*rj*rk
    
#     B1 = 2*drdn*(A1 + A2 - A3)
#     B2 = 2*mu*(N[i,0]*rj*rk + N[j,0]*ri*rk)
#     B3 = (1-2*mu)*(2*N[k,0]*ri*rj + N[j,0]*(i == k) + N[i,0]*(j == k))
#     B4 = (1-4*mu)*N[k,0]*(i == j)
    
#     Fkij[k,i,j] +=  C2*(B1 + B2 + B3 - B4)*w
                            

despTobs = np.zeros(int(len(ObsF)/2))
j=0
for i in range(0,int(len(ObsF)/2)):
    despTobs[i] = m.sqrt(ObsF[j]**2+ObsF[j+1]**2)    
    j+=2


########################################## 
despTfront = np.zeros(int(len(u)/2))
j=0
for i in range(0,int(len(u)/2)):
    despTfront[i] = m.sqrt(u[j]**2+u[j+1]**2)    
    j+=2

##########################################
import matplotlib.pyplot as plt

xcoords = np.zeros(eleT)
ycoords = np.zeros(eleT)
for i in range(eleT):
    xcoords[i] = Elems[i][1]
    ycoords[i] = Elems[i][3]

xcoordsobs = np.zeros(len(Obs))
ycoordsobs = np.zeros(len(Obs))
for i in (range(len(Obs))):
    xcoordsobs[i] = Obs[i][0]
    ycoordsobs[i] = Obs[i][1]

xc = np.zeros(int(len(ObsF)/2))
yc = np.zeros(int(len(ObsF)/2))
aux = 0
for i in range(int(len(ObsF)/2)):
    xc[i] = ObsF[aux]+xcoordsobs[i]
    yc[i] = ObsF[aux+1]+ycoordsobs[i]
    aux += 2

plt.plot(xcoords,ycoords, 'o', color='blue')  
plt.plot(xcoordsobs,ycoordsobs, 'o', color='red')  
plt.plot(xc,yc, 'o', color='orange')  
plt.show()  

