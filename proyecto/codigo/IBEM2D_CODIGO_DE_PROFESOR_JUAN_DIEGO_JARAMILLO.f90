	PROGRAM IBEM2D

!C    INCLUDE 'link_fnl_shared_hpc.h'
!C    INCLUDE 'link_fnl_static_hpc.h'

	USE NUMERICAL_LIBRARIES

    
	CHARACTER*80 Titulo

	REAL*8 G,Nu,K,U,Xi,Xf,Yi,Yf,DesX,DesY,XSP,TraX,TraY,xc,yc,nx,ny,elle
	REAL*8 G11,G12,G21,G22,H11,H12,H21,H22
	REAL*8 Sxx,Sxy,Syx,Syy,DesX1,DesY1,alfa,elleM

	ALLOCATABLE G(:),Nu(:),Nele(:,:),K(:,:),U(:),XSP(:,:),ISP(:,:),NC(:),elleM(:)

      OPEN(5,FILE='IBEM 2D input.prn',STATUS='UNKNOWN')
      OPEN(6,FILE='IBEM 2D  output.sal',STATUS='UNKNOWN')


!C		************************     LEE ALGUNOS DATOS DE ENTRADA		*********************************
!C
!C				Nspace:		Número de espacios del problema.  Si hay más de un espacio, y por supuesto fronteras compartidas,
!C							es necesario que se tenga al menos un desplazamiento prescrito en X y uno en Y en al menos uno de
!C							los espacios del problema.

!C				G(I):		Módulos de cortante de los espacios I
!C				Nu(I):		Módulos de Poisson de los espacios I
!C				Nele(I,1):	Número de elementos de cada espacio I con desplazamientos X y Y prescritos
!C				Nele(I,2):	Número de elementos de cada espacio I con tracciones X y Y prescritos
!C				Nele(I,3):	Número de elementos de cada espacio I con desplazamientos X prescritos y tracciones Y prescritas
!C				Nele(I,4):	Número de elementos de cada espacio I con tracciones X prescritas y desplazamientos Y prescritos 
!C				Nele(I,5):	Número de elementos de cada espacio I con frontera compartida
!C				Nele(I,6):	Número de puntos de observación del espacio I
!C
!C
!C		**************************************************************************************************

	READ(5,'(A80)') Titulo
	WRITE(6,'(A80)') Titulo
	READ(5,'(I10)')Nspace
	WRITE(6,'(I10)')Nspace

	ALLOCATE (G(Nspace),Nu(Nspace),Nele(Nspace,6),NC(Nspace),elleM(Nspace))

	NeleT=0
	NelFc=0
	NCT=0
	DO I=1,Nspace
		READ(5,'(2E10.3,6I10)')G(I),Nu(I),Nele(I,1),Nele(I,2),Nele(I,3),Nele(I,4),Nele(I,5),Nele(I,6)
		WRITE(6,'(2E10.3,6I10)')G(I),Nu(I),Nele(I,1),Nele(I,2),Nele(I,3),Nele(I,4),Nele(I,5),Nele(I,6)
		NeleT=NeleT+Nele(I,1)+Nele(I,2)+Nele(I,3)+Nele(I,4)+Nele(I,5)
		NelFc=NelFc+Nele(I,5)
		NC(I)=0
		
		IF (Nele(I,1).GT. 0)THEN
		    NC(I)=NC(I)+2
		ELSE IF (Nele(I,5) .GT. 0) THEN
		    NC(I)=NC(I)+2
		ELSE IF (Nele(I,3) .GT. 0)THEN
		    NC(I)=NC(I)+1
		    IF(Nele(I,4) .GT. 0)THEN
		        NC(I)=NC(I)+1
		    END IF
		ELSE IF (Nele(I,4) .GT. 0) THEN
		    NC(I)=NC(I)+1
		END IF
		NCT=NCT+NC(I)
	END DO
		
	NH=2*(NeleT+NelFc)

	ALLOCATE (XSP(NeleT,6),ISP(NeleT,3),K(NH+NCT,NH+NCT),U(NH+NCT))

!C		************************     LEE DATOS DE ENTRADA DE CADA UNO DE LOS ELEMENTOS		*********************************
!C
!C				El orden en el que se deben ingresar los elementos es el siguiente:
!C
!C				elementos de cada espacio I con desplazamientos X y Y prescritos
!C				elementos de cada espacio I con tracciones X y Y prescritos
!C				elementos de cada espacio I con desplazamientos X prescritos y tracciones Y prescritas
!C				elementos de cada espacio I con tracciones X prescritas y desplazamientos Y prescritas
!C				elementos de cada espacio I con frontera compartida
!C
!C				Nspace:		Número de espacios del problema
!C
!C				Xi, XSP(J,1):			Coordenada X inicial
!C				Xf, XSP(J,2):			Coordenada X final
!C				Yi, XSP(J,3):			Coordenada Y inicial
!C				Yf, XSP(J,4):			Coordenada Y final
!C
!C              El sistema de referencia x,y debe ser derecho, es decir, que el eje z debe apuntar hacia el observador.
!C				Los elementos de las fronteras externas de todos los espacios se definen en el sentido de las manecillas del reloj
!C              Los elementos de las fronteras internas (agujeros) de los espacios se definen en sentido contrario de las manecillas del reloj
!C
!C				DesX, XSP(J,5):		Desplazamiento prescrito en X
!C				DesY, XSP(J,6):		Desplazamiento prescrito en Y
!C				TraX, XSP(J,5):		Tracción prescrita en X
!C				TraY, XSP(J,6):		Tracción prescrita en Y
!C				IeqF, ISP(J,1):		Número de ecuación de la densidad de fuerza en el elemento
!C				IeqD, ISP(J,2):		Número de ecuación del desplazamiento desconocido en el elemento (solo existe en fronteras compartidas)
!C							El número de la ecuación por un lado de la frontera compartida debe ser igual al número de
!C							la ecuación por el otro lado de la frontera.
!C				IeqT, ISP(J,3):		Número de ecuación de la tracción desconocida en el elemento (solo existe en fronteras compartidas)
!C							El número de la ecuación por un lado de la frontera compartida debe ser igual y de signo contrario al	!C							número de la ecuación por el otro lado de la frontera
!C
!C		********************************************************************************************************************
	READ(5,'(A80)') Titulo
	WRITE(6,'(A80)') Titulo
	J=0
	DO L=1,Nspace
	    elleM(L)=0.0D00
		DO I=1,Nele(L,1)
			READ(5,'(6F10.3,2I10)')Xi,Xf,Yi,Yf,DesX,DesY,IeqF
			J=J+1
			XSP(J,1)=Xi
			XSP(J,2)=Xf
			XSP(J,3)=Yi
			XSP(J,4)=Yf
			elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
			IF(elle .GT. elleM(L)) THEN
			    elleM(L)=elle
			END IF
			XSP(J,5)=DesX
			XSP(J,6)=DesY
			ISP(J,1)=IeqF
		END DO

		DO I=1,Nele(L,2)
			READ(5,'(6F10.3,2I10)')Xi,Xf,Yi,Yf,TraX,TraY,IeqF
			J=J+1
			XSP(J,1)=Xi
			XSP(J,2)=Xf
			XSP(J,3)=Yi
			XSP(J,4)=Yf
			XSP(J,5)=TraX
			XSP(J,6)=TraY
			ISP(J,1)=IeqF
		END DO

		DO I=1,Nele(L,3)
			READ(5,'(6F10.3,2I10)')Xi,Xf,Yi,Yf,DesX,TraY,IeqF
			J=J+1
			XSP(J,1)=Xi
			XSP(J,2)=Xf
			XSP(J,3)=Yi
			XSP(J,4)=Yf
			elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
			IF(elle .GT. elleM(L)) THEN
			    elleM(L)=elle
			END IF
			XSP(J,5)=DesX
			XSP(J,6)=TraY
			ISP(J,1)=IeqF
		END DO

		DO I=1,Nele(L,4)
			READ(5,'(6F10.3,2I10)')Xi,Xf,Yi,Yf,TraX,DesY,IeqF
			J=J+1
			XSP(J,1)=Xi
			XSP(J,2)=Xf
			XSP(J,3)=Yi
			XSP(J,4)=Yf
			elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
			IF(elle .GT. elleM(L)) THEN
			    elleM(L)=elle
			END IF
			XSP(J,5)=TraX
			XSP(J,6)=DesY
			ISP(J,1)=IeqF
		END DO

		DO I=1,Nele(L,5)
			READ(5,'(4F10.3,3I10)')Xi,Xf,Yi,Yf,IeqD,IeqT,IeqF
			J=J+1
			XSP(J,1)=Xi
			XSP(J,2)=Xf
			XSP(J,3)=Yi
			XSP(J,4)=Yf
			elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
			IF(elle .GT. elleM(L)) THEN
			    elleM(L)=elle
			END IF
			ISP(J,2)=IeqD
			ISP(J,3)=IeqT
			ISP(J,1)=IeqF
		END DO
	END DO	

!C		**************		Se inicia el ensamble de la matriz de ecuaciones y del vector de cargas	******************

	DO I=1,NH+NCT
		U(I)=0.0D00
		DO J=1,NH+NCT
			K(I,J)=0.0D00
		END DO
	END DO

	JI=0
	JKI=0
	JIE=0
	NCS=0
	DO L=1,Nspace
		Nel=Nele(L,1)+Nele(L,2)+Nele(L,3)+Nele(L,4)+Nele(L,5)
	    alfa=5.0D00*G(L)/elleM(L)
		
!C		*************	Se ensamblan las ecuaciones del espacio que garantizan que la suma de fuerzas ************
!C		*************	equivalentes sobre la frontera son iguales a cero							  ************

        IF(NC(L) .GT. 0) THEN
            IF(NC(L) .GT. 1) THEN
		        DO J=1,Nel
			        JK=JKI+J
			        Xi=XSP(JK,1)
			        Xf=XSP(JK,2)
			        Yi=XSP(JK,3)
			        Yf=XSP(JK,4)
			        elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)

			        K(NH+NCS+1,2*ISP(JK,1)-1)=elle/elleM(L)
			        K(NH+NCS+2,2*ISP(JK,1))=elle/elleM(L)
		        END DO
		    ELSE IF(Nele(L,3) .GT. 0)THEN
		        DO J=1,Nel
			        JK=JKI+J
			        Xi=XSP(JK,1)
			        Xf=XSP(JK,2)
			        Yi=XSP(JK,3)
			        Yf=XSP(JK,4)
			        elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)

			        K(NH+NCS+1,2*ISP(JK,1)-1)=elle/elleM(L)
		        END DO
		    ELSE IF(Nele(L,4) .GT. 0)THEN
		        DO J=1,Nel
			        JK=JKI+J
			        Xi=XSP(JK,1)
			        Xf=XSP(JK,2)
			        Yi=XSP(JK,3)
			        Yf=XSP(JK,4)
			        elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)

			        K(NH+NCS+1,2*ISP(JK,1))=elle/elleM(L)
		        END DO
		    END IF
		END IF	

!C		**************		Se ensamblan las fronteras con ambos desplazamientos prescritos	******************

		DO I=1,Nele(L,1)

			JI=JI+1
			JIE=JIE+1
			
			xc=(XSP(JI,1)+XSP(JI,2))/2
			yc=(XSP(JI,3)+XSP(JI,4))/2

			U(2*JIE-1)=XSP(JI,5)*alfa
			U(2*JIE)=XSP(JI,6)*alfa

			K(2*JIE-1,NH+NCS+1)=1.0D00
			K(2*JIE,NH+NCS+2)=1.0D00			
	
			DO J=1,Nel

				JK=JKI+J
				Xi=XSP(JK,1)
				Xf=XSP(JK,2)
				Yi=XSP(JK,3)
				Yf=XSP(JK,4)
				
				CALL IGREENG(xc,yc,Xi,Xf,Yi,Yf,G(L),Nu(L),G11,G12,G21,G22)

				K(2*JIE-1,2*ISP(JK,1)-1)=G11*alfa
				K(2*JIE-1,2*ISP(JK,1))=G12*alfa
				K(2*JIE,2*ISP(JK,1)-1)=G21*alfa
				K(2*JIE,2*ISP(JK,1))=G22*alfa

			END DO
		END DO

!C		**************		Se ensamblan las fronteras con ambas tracciones prescritas	******************

		DO I=1,Nele(L,2)

			JI=JI+1
			JIE=JIE+1
			
			Xi=XSP(JI,1)
			Xf=XSP(JI,2)
			Yi=XSP(JI,3)
			Yf=XSP(JI,4)
			xc=(Xi+Xf)/2
			yc=(Yi+Yf)/2
			elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
			nx=(Yi-Yf)/elle
			ny=(Xf-Xi)/elle

			U(2*JIE-1)=XSP(JI,5)
			U(2*JIE)=XSP(JI,6)

			DO J=1,Nel

				JK=JKI+J
				Xi=XSP(JK,1)
				Xf=XSP(JK,2)
				Yi=XSP(JK,3)
				Yf=XSP(JK,4)
				
			CALL IGREENH(xc,yc,Xi,Xf,Yi,Yf,nx,ny,G(L),Nu(L),H11,H12,H21,H22,0)

				K(2*JIE-1,2*ISP(JK,1)-1)=H11
				K(2*JIE-1,2*ISP(JK,1))=H12
				K(2*JIE,2*ISP(JK,1)-1)=H21
				K(2*JIE,2*ISP(JK,1))=H22

			END DO
		END DO			

!C		**************		Se ensamblan las fronteras con desplazamientos X prescritos y tracciones Y prescritas	******************

		DO I=1,Nele(L,3)

			JI=JI+1
			JIE=JIE+1
			
			Xi=XSP(JI,1)
			Xf=XSP(JI,2)
			Yi=XSP(JI,3)
			Yf=XSP(JI,4)
			xc=(Xi+Xf)/2
			yc=(Yi+Yf)/2
			elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
			nx=(Yi-Yf)/elle
			ny=(Xf-Xi)/elle

			U(2*JIE-1)=XSP(JI,5)*alfa
			U(2*JIE)=XSP(JI,6)

			K(2*JIE-1,NH+NCS+1)=1.0D00

			DO J=1,Nel

				JK=JKI+J
				Xi=XSP(JK,1)
				Xf=XSP(JK,2)
				Yi=XSP(JK,3)
				Yf=XSP(JK,4)
								
				CALL IGREENG(xc,yc,Xi,Xf,Yi,Yf,G(L),Nu(L),G11,G12,G21,G22)

				K(2*JIE-1,2*ISP(JK,1)-1)=G11*alfa
				K(2*JIE-1,2*ISP(JK,1))=G12*alfa

			CALL IGREENH(xc,yc,Xi,Xf,Yi,Yf,nx,ny,G(L),Nu(L),H11,H12,H21,H22,0)

				K(2*JIE,2*ISP(JK,1)-1)=H21
				K(2*JIE,2*ISP(JK,1))=H22

			END DO
		END DO

!C		**************		Se ensamblan las fronteras con tracciones X prescritas y desplazamientos Y prescritos	******************

		DO I=1,Nele(L,4)

			JI=JI+1
			JIE=JIE+1
			
			Xi=XSP(JI,1)
			Xf=XSP(JI,2)
			Yi=XSP(JI,3)
			Yf=XSP(JI,4)
			xc=(Xi+Xf)/2
			yc=(Yi+Yf)/2
			elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
			nx=(Yi-Yf)/elle
			ny=(Xf-Xi)/elle

			U(2*JIE-1)=XSP(JI,5)
			U(2*JIE)=XSP(JI,6)*alfa
            
            IF(NC(L) .GT. 1)THEN
			    K(2*JIE,NH+NCS+2)=1.0D00
			ELSE
			    K(2*JIE,NH+NCS+1)=1.0D00
			END IF

			DO J=1,Nel

				JK=JKI+J
				Xi=XSP(JK,1)
				Xf=XSP(JK,2)
				Yi=XSP(JK,3)
				Yf=XSP(JK,4)
								
				CALL IGREENG(xc,yc,Xi,Xf,Yi,Yf,G(L),Nu(L),G11,G12,G21,G22)

				K(2*JIE,2*ISP(JK,1)-1)=G21*alfa
				K(2*JIE,2*ISP(JK,1))=G22*alfa

			CALL IGREENH(xc,yc,Xi,Xf,Yi,Yf,nx,ny,G(L),Nu(L),H11,H12,H21,H22,0)

				K(2*JIE-1,2*ISP(JK,1)-1)=H11
				K(2*JIE-1,2*ISP(JK,1))=H12

			END DO
		END DO

!C		**************		Se ensamblan las fronteras compartidas	******************

		DO I=1,Nele(L,5)

			JI=JI+1
			JIE=JIE+1
			
			Xi=XSP(JI,1)
			Xf=XSP(JI,2)
			Yi=XSP(JI,3)
			Yf=XSP(JI,4)
			xc=(Xi+Xf)/2
			yc=(Yi+Yf)/2
			elle=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
			nx=(Yi-Yf)/elle
			ny=(Xf-Xi)/elle

			U(2*JIE-1)=0.0D00
			U(2*JIE)=0.0D00

			K(2*JIE-1,2*ISP(JI,2)-1)=-1.0D00*alfa
			K(2*JIE,2*ISP(JI,2))=-1.0D00*alfa

			K(2*JIE-1,NH+NCS+1)=1.0D00
			K(2*JIE,NH+NCS+2)=1.0D00

			DO J=1,Nel

				JK=JKI+J
				Xi=XSP(JK,1)
				Xf=XSP(JK,2)
				Yi=XSP(JK,3)
				Yf=XSP(JK,4)

				CALL IGREENG(xc,yc,Xi,Xf,Yi,Yf,G(L),Nu(L),G11,G12,G21,G22)

				K(2*JIE-1,2*ISP(JK,1)-1)=G11*alfa
				K(2*JIE-1,2*ISP(JK,1))=G12*alfa
				K(2*JIE,2*ISP(JK,1)-1)=G21*alfa
				K(2*JIE,2*ISP(JK,1))=G22*alfa

			END DO

			JIE=JIE+1

			U(2*JIE-1)=0.0D00
			U(2*JIE)=0.0D00

			ISIG=SIGN(1,ISP(JI,3))
			
			K(2*JIE-1,2*ISP(JI,3)*ISIG-1)=-1.0D00*ISIG
			K(2*JIE,2*ISP(JI,3)*ISIG)=-1.0D00*ISIG

			DO J=1,Nel

				JK=JKI+J
				Xi=XSP(JK,1)
				Xf=XSP(JK,2)
				Yi=XSP(JK,3)
				Yf=XSP(JK,4)
				ISIG=SIGN(1,ISP(JK,3))
							
			CALL IGREENH(xc,yc,Xi,Xf,Yi,Yf,nx,ny,G(L),Nu(L),H11,H12,H21,H22,0)

				K(2*JIE-1,2*ISP(JK,1)-1)=H11
				K(2*JIE-1,2*ISP(JK,1))=H12
				K(2*JIE,2*ISP(JK,1)-1)=H21
				K(2*JIE,2*ISP(JK,1))=H22

			END DO
		END DO
		JKI=JKI+Nel
		NCS=NCS+NC(L)	    
	END DO			

	CALL DLSLRG (NH+NCT,K,NH+NCT,U,1,U)

!C		************************     LEE COORDENADAS DE LOS PUNTOS DE OBSERVACIÓN Y CALCULA 	*********************************
!C		************************	 DESPLAZAMIENTOS Y TENSOR DE ESFUERZOS      **************************************************
!C
!C				Las coordenadas de los puntos de observación hay que especificarlas en el mismo orden
!C				en que se especificaron los espacios.
!C
!C              Si no se especificaron desplazamientos en la frontera, los desplazamientos impresos en los puntos de observación 
!C              quedan referidos a un "datum" o referencia arbitraria.
!C
!C		***********************************************************************************************************************

	READ(5,'(A80)') Titulo
	WRITE(6,'(A80)') Titulo

	Jo=0
	NCS=0
	DO L=1,Nspace
		Nel=Nele(L,1)+Nele(L,2)+Nele(L,3)+Nele(L,4)+Nele(L,5)
		alfa=5.0D00*G(L)/elleM(L)
		IF(NC(L) .GT. 0) THEN
		    IF(NC(L) .GT. 1) THEN
		        DesX1=U(NH+NCS+1)/alfa
		        DesY1=U(NH+NCS+2)/alfa
		    ELSE IF(Nele(L,3) .GT. 0) THEN
		            DesX1=U(NH+NCS+1)/alfa
		            DesY1=0.0D00
		        ELSE
		            DesX1=0.0D00
		            DesY1=U(NH+NCS+1)/alfa
		    END IF
		ELSE
		    DesX1=0.0
		    DesY1=0.0
		END IF
		
		DO I=1,Nele(L,6)
		    READ(5,'(2F10.3)')xc,yc
		    
		    Sxx=0.0D00
		    Sxy=0.0D00
		    Syx=0.0D00
		    Syy=0.0D00
		    
		    DesX=DesX1
		    DesY=DesY1
		    
			DO JK=1,Nel
				J=Jo+JK
				Xi=XSP(J,1)
				Xf=XSP(J,2)
				Yi=XSP(J,3)
				Yf=XSP(J,4)

				CALL IGREENG(xc,yc,Xi,Xf,Yi,Yf,G(L),Nu(L),G11,G12,G21,G22)

				DesX=DesX+G11*U(2*ISP(J,1)-1)+G12*U(2*ISP(J,1))
				DesY=DesY+G21*U(2*ISP(J,1)-1)+G22*U(2*ISP(J,1))

				CALL IGREENH(xc,yc,Xi,Xf,Yi,Yf,1.0D00,0.0D00,G(L),Nu(L),H11,H12,H21,H22,1)

				Sxx=Sxx+H11*U(2*ISP(J,1)-1)+H12*U(2*ISP(J,1))
				Sxy=Sxy+H21*U(2*ISP(J,1)-1)+H22*U(2*ISP(J,1))

				CALL IGREENH(xc,yc,Xi,Xf,Yi,Yf,0.0D00,1.0D00,G(L),Nu(L),H11,H12,H21,H22,1)

				Syx=Syx+H11*U(2*ISP(J,1)-1)+H12*U(2*ISP(J,1))
				Syy=Syy+H21*U(2*ISP(J,1)-1)+H22*U(2*ISP(J,1))

			END DO

!C		************************     IMPRIME RESULTADOS DE LOS PUNTOS DE OBSERVACIÓN 	*********************************
!C
!C				Los resultados de los puntos de observación se imprimen en espacios de 10 columnas con
!C				3 decimales únicamente. Si se quiere tener más resolución en los resultados se deben cambiar
!C				las unidades en las que está descrito el problema, de tal manera que los resultados tengan 
!C				parte entera y parte decimal
!C
!C				Los resultados se escriben en el siguiente orden:
!C
!C				Xob, Yob, DespX, DespY, SigmaXX, TaoXY, TaoYX, SigmaYY
!C
!C		***********************************************************************************************************************


			WRITE(6,'(8F10.3)')xc,yc,DesX,DesY,Sxx,Sxy,Syx,Syy

		END DO
		Jo=Jo+Nel
		NCS=NCS+NC(L)	    
	END DO

	CLOSE (5)
	CLOSE (6)

	STOP
	END PROGRAM

!C			**************************		SUBROUTINE IGREENG	***********************************

	SUBROUTINE IGREENG(xc,yc,Xi,Xf,Yi,Yf,G,Nu,G11,G12,G21,G22)

	REAL*8 G,Nu,xc,yc,Xi,Xf,Yi,Yf,G11,G12,G21,G22
	REAL*8 ele,nx,ny,Cte1,Cte2,zxc,zyc,PI,S13,S2,W13,W2
	REAL*8 zx1,zy1,zx2,zy2,zx3,zy3,rx1,ry1,rx2,ry2,rx3,ry3,r1,r2,r3
	REAL*8 G11p1,G12p1,G21p1,G22p1,G11p2,G12p2,G21p2,G22p2
	REAL*8 G11p3,G12p3,G21p3,G22p3
	
	PI=4.0D00*DATAN(1.0D00)

	zxc=(Xi+Xf)/2
	zyc=(Yi+Yf)/2

	ele=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
	nx=(Xf-Xi)/ele
	ny=(Yf-Yi)/ele

	Cte1=1.0/8/PI/G/(1-Nu)
	Cte2=3.0-4.0*Nu

	IF(DSQRT((xc-zxc)**2+(yc-zyc)**2) .LE. ele/10.0) THEN
				
		G11=Cte1*ele*(Cte2*(1.0-DLOG(ele/2.0))+nx*nx)
		G22=Cte1*ele*(Cte2*(1.0-DLOG(ele/2.0))+ny*ny)
		G12=Cte1*ele*nx*ny
		G21=G12

	ELSE

		S13=0.774596669241483
		S2=0.0D00
		W13=0.555555555555556
		W2=0.8888888888888889

		zx1=zxc-nx*S13*ele/2.0
		zy1=zyc-ny*S13*ele/2.0

		zx2=zxc
		zy2=zyc

		zx3=zxc+nx*S13*ele/2.0
		zy3=zyc+ny*S13*ele/2.0

		rx1=xc-zx1
		ry1=yc-zy1
		r1=DSQRT(rx1**2+ry1**2)

		rx2=xc-zx2
		ry2=yc-zy2
		r2=DSQRT(rx2**2+ry2**2)

		rx3=xc-zx3
		ry3=yc-zy3
		r3=DSQRT(rx3**2+ry3**2)

		G11p1=Cte1*(Cte2*DLOG(1.0/r1)+rx1*rx1/r1/r1)
		G22p1=Cte1*(Cte2*DLOG(1.0/r1)+ry1*ry1/r1/r1)
		G12p1=Cte1*(rx1*ry1/r1/r1)
		G21p1=G12p1

		G11p2=Cte1*(Cte2*DLOG(1.0/r2)+rx2*rx2/r2/r2)
		G22p2=Cte1*(Cte2*DLOG(1.0/r2)+ry2*ry2/r2/r2)
		G12p2=Cte1*(rx2*ry2/r2/r2)
		G21p2=G12p2

		G11p3=Cte1*(Cte2*DLOG(1.0/r3)+rx3*rx3/r3/r3)
		G22p3=Cte1*(Cte2*DLOG(1.0/r3)+ry3*ry3/r3/r3)
		G12p3=Cte1*(rx3*ry3/r3/r3)
		G21p3=G12p3

		G11=(G11p1*W13+G11p2*W2+G11p3*W13)*ele/2.0
		G22=(G22p1*W13+G22p2*W2+G22p3*W13)*ele/2.0
		G12=(G12p1*W13+G12p2*W2+G12p3*W13)*ele/2.0
		G21=(G21p1*W13+G21p2*W2+G21p3*W13)*ele/2.0

	END IF

	RETURN
	END

!C			**************************		SUBROUTINE IGREENH	***********************************

	SUBROUTINE IGREENH(xc,yc,Xi,Xf,Yi,Yf,n1,n2,G,Nu,H11,H12,H21,H22,Iob)

	REAL*8 G,Nu,xc,yc,Xi,Xf,Yi,Yf,H11,H12,H21,H22,zxc,zyc
	REAL*8 ele,nx,ny,Cte1,Cte2,PI,S13,S2,W13,W2
	REAL*8 zx1,zy1,zx2,zy2,zx3,zy3,rx1,ry1,rx2,ry2,rx3,ry3,r1,r2,r3
	REAL*8 H11p1,H12p1,H21p1,H22p1,H11p2,H12p2,H21p2,H22p2,n1,n2
	REAL*8 H11p3,H12p3,H21p3,H22p3

	PI=4.0D00*DATAN(1.0D00)

	Cte1=-1.0/4/PI/(1-Nu)
	Cte2=1.0-2.0*Nu

	zxc=(Xi+Xf)/2
	zyc=(Yi+Yf)/2

	ele=DSQRT((Xf-Xi)**2+(Yf-Yi)**2)
	nx=(Xf-Xi)/ele
	ny=(Yf-Yi)/ele

	IF(DSQRT((xc-zxc)**2+(yc-zyc)**2) .LE. ele/10.0) THEN

		IF(Iob .EQ. 1) THEN

			IF(n1 .EQ. 1.0) THEN
				H11=-ny*(nx*nx+ny*ny/2)
				H12=-ny*ny*nx/2
				H21=nx*nx*nx/2
				H22=-ny*ny*ny/2
			ELSE
				H11=nx*nx*nx/2
				H12=-ny*ny*ny/2
				H21=ny*nx*nx/2
				H22=nx*(ny*ny+nx*nx/2)
			END IF

		ELSE

			H11=0.5D00
			H22=0.5D00
			H12=0.0D00
			H21=0.0D00

		END IF

	ELSE

		S13=0.774596669241483
		S2=0.0D00
		W13=0.555555555555556
		W2=0.8888888888888889

		zx1=zxc+nx*S13*ele/2.0
		zy1=zyc+ny*S13*ele/2.0

		zx2=zxc
		zy2=zyc

		zx3=zxc-nx*S13*ele/2.0
		zy3=zyc-ny*S13*ele/2.0

		rx1=xc-zx1
		ry1=yc-zy1
		r1=DSQRT(rx1**2+ry1**2)

		rx2=xc-zx2
		ry2=yc-zy2
		r2=DSQRT(rx2**2+ry2**2)

		rx3=xc-zx3
		ry3=yc-zy3
		r3=DSQRT(rx3**2+ry3**2)

		rn1=rx1*n1+ry1*n2
		rn2=rx2*n1+ry2*n2
		rn3=rx3*n1+ry3*n2

		H11p1=Cte1/r1/r1*(rn1*(Cte2+2*rx1*rx1/r1/r1))
		H22p1=Cte1/r1/r1*(rn1*(Cte2+2*ry1*ry1/r1/r1))
		H12p1=Cte1/r1/r1*(rn1*(2*rx1*ry1/r1/r1)+Cte2*(rx1*n2-ry1*n1))
		H21p1=Cte1/r1/r1*(rn1*(2*ry1*rx1/r1/r1)+Cte2*(ry1*n1-rx1*n2))

		H11p2=Cte1/r2/r2*(rn2*(Cte2+2*rx2*rx2/r2/r2))
		H22p2=Cte1/r2/r2*(rn2*(Cte2+2*ry2*ry2/r2/r2))
		H12p2=Cte1/r2/r2*(rn2*(2*rx2*ry2/r2/r2)+Cte2*(rx2*n2-ry2*n1))
		H21p2=Cte1/r2/r2*(rn2*(2*ry2*rx2/r2/r2)+Cte2*(ry2*n1-rx2*n2))

		H11p3=Cte1/r3/r3*(rn3*(Cte2+2*rx3*rx3/r3/r3))
		H22p3=Cte1/r3/r3*(rn3*(Cte2+2*ry3*ry3/r3/r3))
		H12p3=Cte1/r3/r3*(rn3*(2*rx3*ry3/r3/r3)+Cte2*(rx3*n2-ry3*n1))
		H21p3=Cte1/r3/r3*(rn3*(2*ry3*rx3/r3/r3)+Cte2*(ry3*n1-rx3*n2))

		H11=(H11p1*W13+H11p2*W2+H11p3*W13)*ele/2.0
		H22=(H22p1*W13+H22p2*W2+H22p3*W13)*ele/2.0
		H12=(H12p1*W13+H12p2*W2+H12p3*W13)*ele/2.0
		H21=(H21p1*W13+H21p2*W2+H21p3*W13)*ele/2.0

	END IF

	RETURN
	END
