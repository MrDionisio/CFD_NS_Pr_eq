        Program Pr
        Implicit none
 
        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        INTEGER NI, NJ, NITER, SMAX
        INTEGER I,J
        REAL(8) L,H,U0,MU,Nu,R0,P0
        REAL(8) dx,dy,CFL,EPS, C_f, Re_x
        REAL(8),ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL(8),ALLOCATABLE :: U_n(:,:),V_n(:,:),P_n(:,:),R_n(:,:)
        REAL(8),ALLOCATABLE :: U(:,:),V(:,:),P(:,:),R(:,:)
        write(*,*) 'Read input file'
        open(IO,FILE='Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) EPS
        read(IO,*) SMAX

        read(IO,*) U0
        read(IO,*) MU
        read(IO,*) R0
        read(IO,*) P0
        CLOSE(IO)
   
        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates

!----------------- Node variables -----------------------------
        allocate(U_n(NI,NJ))  ! Velocity U
        allocate(V_n(NI,NJ))  ! Velocity V
        allocate(P_n(NI,NJ))  ! Pressure
        allocate(U(NI,NJ))  
        allocate(V(NI,NJ))  
        allocate(P(NI,NJ))

!----------------- Coordinate of nodes ------------------------
        dx=L/(NI-1)
        dy=H/(NJ-1)

        DO I=1,NI
          DO J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          END DO
        END DO

!----------------- Parameters ------------------------

        NU=MU/R0

        write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
        write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
        write(*,*)'ReL= ', U0*L/NU
      

!----------------- Initial fields -----------------------------

        DO I=1,NI
          DO J=1,NJ
            U_n(I,J)=U0
            V_n(I,J)=1.0e-5
            P_n(I,J)=P0
          ENDDO
        ENDDO

!---------------- Solve Prandtl equations ---------------------
 
        write(*,*) 'Solve Prandtl equations'      
        call  Prandtl(U_n,V_n,NI,NJ, U0, EPS, dx, dy, NU, SMAX)

        write(*,*) 'Output data' 
        Open(IO,FILE='C_f3.txt')
                do i=1,NI
                        Re_x = U0*dx*i/NU
                        C_f=0.664/sqrt(Re_x)
                        write(IO,*) Re_x , C_f
                end do
        Close(IO)


        Open(IO,FILE='x=0.1.txt')
                i=0.1/dx
                print*, i
                do j=1,NJ
                        write(IO,*) U_n(i,j), Y_Node(i,j) 
                end do
        Close(IO)

        Open(IO,FILE='x=0.2.txt')
                i=0.2/dx
                print*, i
                do j=1,NJ
                        write(IO,*) U_n(i,j), Y_Node(i,j) 
                end do
        Close(IO)

        Open(IO,FILE='x=0.5.txt')
                i=0.5/dx
                print*, i
                do j=1,NJ
                        write(IO,*) U_n(i,j), Y_Node(i,j) 
                end do
        Close(IO)




 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data' 
        Open(IO,FILE='Results.plt')
        Call OUTPUT(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
        Close(IO)
     
        END PROGRAM

   SUBROUTINE OUTPUT(IO,NI,NJ,X,Y,U,V,P)
               IMPLICIT NONE
         INTEGER NI,NJ,IO
         REAL(8),DIMENSION(NI,NJ):: X,Y
         REAL(8),DIMENSION(NI,NJ):: U,V,P

         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI,1:NJ)
         Write(IO,'(100E25.16)') V(1:NI,1:NJ)
         Write(IO,'(100E25.16)') P(1:NI,1:NJ)
 	CLOSE(1)
      END SUBROUTINE

   SUBROUTINE Prandtl(U,V,NI,NJ, U0, EPS, dx, dy, NU, SMAX)
      IMPLICIT NONE
	INTEGER NI,NJ, j, i, s, SMAX
        REAL(8) U0, ERR, EPS, dx, dy, NU
        REAL(8),DIMENSION(NI,NJ):: U,V
        REAL(8), DIMENSION (NJ) :: A, B, C, D, US, VS, US1, VS1, U_ni, alpha, beta
        
        A = 0
        B = 1
        C = 0
        D = 0 
        D(NJ) = U0
        do i=2,NI
                ERR=1
                U_ni=U(i-1,:)
                VS=V(i-1,:)
                US=U_ni
                US(1)=0
                US1(NJ)=U0
                VS(1)=0
                s=0
                DO WHILE ((ERR>=EPS).and.(s<SMAX))

                        alpha(1) = - C(1) / B(1)
                        beta(1) = D(1) / B(1)

                        do j = 2, NJ-1, 1
                                A(j) = -(VS(j-1) / 2.0 / dy) - (NU / (dy**2))
                                B(j) = US(j) / dx  + 2.0 * NU / (dy**2)
                                C(j) = (VS(j+1) / 2.0 / dy) - (NU / (dy**2))
                                D(j) = U_ni(j) * U_ni(j) / dx
                                alpha(j) = - C(j) / (A(j) * alpha(j-1) + B(j))
                                beta(j) = (D(j) - A(j) * beta(j-1))/ (A(j) * alpha(j-1) + B(j))
                        end do

                        do j = NJ - 1, 2, -1
                                US1(j) = alpha(j) * US1(j+1) + beta(j)
                        end do

                        do j = 2, NJ, 1
                                VS1(j) = VS1(j - 1) - dy / 2.0 * ((US1(j) - U_ni(j)) / dx + (US1(j-1) - U_ni(j-1)) / dx)
                        end do

                        ERR=MAX(MAXVAL(ABS(US1 - US)), MAXVAL(ABS(VS1 - VS)))
                        s=s+1
                        US=US1
                        VS=VS1
                end do
                print*, ERR, EPS, s
                U(i, :) = US
                V(i, :) = VS
        end do



      END SUBROUTINE
