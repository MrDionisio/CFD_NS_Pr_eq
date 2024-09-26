        Program Pr
        Implicit none
 
        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        INTEGER NI, NJ, NITER, SMAX
        INTEGER I,J
        REAL(8) L,H,U0,MU,NU,R0,P0
        REAL(8) dx,dy,CFL,EPS,dt, C_f
        REAL(8),ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
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
        read(IO,*) CFL
        CLOSE(IO)
   
        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates

!----------------- Node variables -----------------------------
        allocate(U(0:NI,0:NJ))  
        allocate(V(0:NI,0:NJ))  
        allocate(P(0:NI,0:NJ))

!----------------- Coordinate of nodes ------------------------
        dx=L/(NI-1)
        dy=H/(NJ-1)
        dt=CFL/U0*MIN(dx, dy)
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
        write(*,*)'ReL= ', U0*L/NU, 'CFL = ', CFL, "dt = ", dt
      

!----------------- Initial fields -----------------------------

        DO I=0,NI
          DO J=0,NJ
            U(I,J)=U0
            V(I,J)=1.0e-5
            P(I,J)=0
          ENDDO
        ENDDO

!---------------- Solve Prandtl equations ---------------------
 
        write(*,*) 'Solve Novie-Stocks equations'      
        call NOVIE_STOCKS(U, V, P, dx, dy, dt, EPS, NI, NJ, NU, U0, SMAX)


        Open(IO,FILE='C_f.txt')
        do i=1,NI-1
                C_f=MU*(U(i,1)-U(i,0))/dy/R0/MAXVAL(U(i,:))/MAXVAL(U(i,:))*2
                write(IO,*) MAXVAL(U(i,:))*(dx*(i-0.5))/NU, C_f
        end do
        Close(IO)

        Open(IO,FILE='Nx=0.1.txt')
                i=0.1/dx
                print*, i
                do j=0,NJ
                        write(IO,*) U(i,j), (j-0.5)*dy
                end do
        Close(IO)

        Open(IO,FILE='Nx=0.2.txt')
                i=0.2/dx
                print*, i
                do j=0,NJ
                        write(IO,*) U(i,j), (j-0.5)*dy 
                end do
        Close(IO)

        Open(IO,FILE='Nx=0.5.txt')
                i=0.5/dx
                print*, i
                do j=0,NJ
                        write(IO,*) U(i,j), (j-0.5)*dy
                end do
                 
        Close(IO)
 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data' 
        Open(IO,FILE='Results.plt')
        Call OUTPUT(IO,NI,NJ,X_Node,Y_Node,U,V,P)
        Close(IO)
     
        END PROGRAM

        SUBROUTINE OUTPUT(IO,NI,NJ,X,Y,U,V,P)
                IMPLICIT NONE
                INTEGER NI,NJ,IO
                REAL(8),DIMENSION(NI,NJ):: X,Y
                REAL(8),DIMENSION(0:NI,0:NJ):: U,V,P

                Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
                Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
                Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
                Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
                Write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
                Write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
                Write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)
                CLOSE(1)
        END SUBROUTINE

        SUBROUTINE NOVIE_STOCKS(U, V, P, dx, dy, dt, EPS, NI, NJ, NU, U0, SMAX)
                real(8) EPS, dx, dy, dt, NU, U0, ERR_U, ERR_V, ERR_P
                real(8) on, against, Half
                INTEGER NI, NJ, i, j, SMAX, S
                real(8), dimension(0:NI, 0:NJ) :: U, V, P, UN, VN, PN, R_U, R_V, R_P 
                real(8), dimension (3) :: left, right, up, down ! 1 - U; 2 - V; 3 - P
                real(8), dimension (2) :: bar_l, bar_r, bar_u, bar_d ! 1 - U; 2 - V


                ERR_U = 0
                ERR_V = 0
                ERR_P = EPS+1
                S = 0 

                U(: , 0) = - U(:, 1)
                V(: , 0) = - V(:, 1)
                P(: , 0) = P(:, 1)
                U(NI, :) = U(NI-1, :)
                V(NI, :) = V(NI-1, :)
                P(NI, :) = 0.0
                U(0, :) = U0
                V(0, :) = 0.0
                P(0, :) = P(1, :)
                U(:, NJ) = U(:, NJ-1)
                V(:, NJ) = V(:, NJ-1)
                P(:, NJ) = 0.0

                OPEN(1, file="R_U.txt")
                OPEN(2, file="R_V.txt")
                OPEN(3, file="R_P.txt")

                DO WHILE ((MAX(ERR_U, MAX(ERR_V, ERR_P))>EPS).AND.(S<SMAX))

                        U(: , 0) = - U(:, 1)
                        V(: , 0) = - V(:, 1)
                        P(: , 0) = P(:, 1)
                        U(NI, :) = U(NI-1, :)
                        V(NI, :) = V(NI-1, :)
                        P(NI, :) = 0.0
                        U(0, :) = U0
                        V(0, :) = 0.0
                        P(0, :) = P(1, :)
                        U(:, NJ) = U(:, NJ-1)
                        V(:, NJ) = V(:, NJ-1)
                        P(:, NJ) = 0.0

                        do i=1,NI-1
                                do j=1,NJ-1

                                        call bar(U, V, bar_d, bar_l, bar_r, bar_u, NI, NJ, i, j)
                                        
                                        left(1) = Half(bar_l(1), U(i,j), U(i-1,j))!on(U, bar_l(1), NI, NJ, i-1, j)
                                        left(2) = Half(bar_l(1), V(i,j), V(i-1,j))!on(V, bar_l(2), NI, NJ, i-1, j)
                                        left(3) = Half(bar_l(1), P(i-1,j), P(i,j))!against(P, bar_l(1), NI, NJ, i-1, j)

                                        right(1) = Half(bar_r(1), U(i+1,j), U(i,j))!on(U, bar_r(1), NI, NJ, i, j)
                                        right(2) = Half(bar_r(1), V(i+1,j), V(i,j))!on(V, bar_r(2), NI, NJ, i, j)
                                        right(3) = Half(bar_r(1), P(i,j), P(i+1,j))!against(P, bar_r(1), NI, NJ, i, j)

                                        up(1) = Half(bar_u(2), U(i,j+1), U(i,j))!on(U, bar_u(1), NI, NJ, i, j)
                                        up(2) = Half(bar_u(2), V(i,j+1), V(i,j))!on(V, bar_u(2), NI, NJ, i, j)
                                        up(3) = Half(bar_u(2), P(i,j), P(i,j+1))!against(P, bar_d(2), NI, NJ, i, j)

                                        down(1) = Half(bar_d(2), U(i,j), U(i,j-1))!on(U, bar_d(1), NI, NJ, i, j-1)
                                        down(2) = Half(bar_d(2), V(i,j), V(i,j-1))!on(V, bar_d(2), NI, NJ, i, j-1)
                                        down(3) = Half(bar_d(2), P(i,j-1), P(i,j))!against(P, bar_d(2), NI, NJ, i, j-1)
                                        
                                        

                                        if (j==1) then 
                                                R_P(i,j) = - 1/U0**2 * ((right(1)-left(1))/dx &
                                                + (up(2))/dy)
                                        else 
                                                R_P(i,j) = - 1/U0**2 * ((right(1)-left(1))/dx &
                                                + (up(2)-down(2))/dy)
                                        end if

                                        R_U(i,j) = (bar_r(1)*right(1)-bar_l(1)*left(1))/dx &
                                        + (bar_u(2)*up(1)-bar_d(2)*down(1))/dy &
                                        + (right(3)-left(3))/dx &
                                        - (U(i+1, j)-2*U(i,j)+U(i-1,j)) * NU / dx**2 &
                                        - (U(i,j+1)-2*U(i,j)+U(i, j-1)) * NU / dy**2

                                        R_V(i,j) = (bar_r(1)*right(2)-bar_l(1)*left(2))/dx &
                                        + (bar_u(2)*up(2)-bar_d(2)*down(2))/dy &
                                        + (up(3)-down(3))/dy &
                                        - (V(i+1, j)-2*V(i,j)+V(i-1,j)) * NU / dx**2 &
                                        - (V(i,j+1)-2*V(i,j)+V(i, j-1)) * NU / dy**2
                                        


                                end do
                        end do
                       
                        P(:,:) = P(:,:) + dt*R_P(:,:)
                        U(:,:) = U(:,:) - dt*R_U(:,:)
                        V(:,:) = V(:,:) - dt*R_V(:,:)
                      
                

                        ERR_U = MAXVAL(ABS(R_U))
                        ERR_V = MAXVAL(ABS(R_V))
                        ERR_P = MAXVAL(ABS(R_P))
                                                

                        S = S+1

                        call Residuals(1, 2, 3, S, ERR_U, ERR_V, ERR_P)
                END DO

                CLOSE(1)
                CLOSE(2)
                CLOSE(3)

        END SUBROUTINE

        SUBROUTINE Residuals(IO1, IO2, IO3, S, ERR1, ERR2, ERR3)
                integer IO1, IO2, IO3, S
                REAL(8) ERR1, ERR2, ERR3
                Write(IO1, *) S, ERR1
                Write(IO2, *) S, ERR2
                Write(IO3, *) S, ERR3
                Print*, "S", S, ERR1, ERR2, ERR3
        END SUBROUTINE

        SUBROUTINE bar(U, V, bar_d, bar_l, bar_r, bar_u, NI, NJ, i, j)
                integer NI, NJ, i, j
                real(8), dimension(0:NI, 0:NJ) :: U, V
                real(8), dimension(2) ::  bar_d, bar_l, bar_r, bar_u

                bar_l(1) = 0.5*(U(i,j)+U(i-1,j))
                bar_l(2) = 0.5*(V(i,j)+V(i-1,j))
                bar_r(1) = 0.5*(U(i,j)+U(i+1,j))
                bar_r(2) = 0.5*(V(i,j)+U(i+1,j))
                bar_d(1) = 0.5*(U(i,j)+U(i,j-1))
                bar_d(2) = 0.5*(V(i,j)+V(i,j-1))
                bar_u(1) = 0.5*(U(i,j)+V(i,j+1))
                bar_u(2) = 0.5*(V(i,j)+V(i,j+1))
        END SUBROUTINE

        real(8) function on(A, bar, NI, NJ, i, j)
                Implicit none
                integer NI, NJ, i, j
                real(8) bar
                real(8), dimension(0:NI, 0:NJ) :: A


                if (bar>= 0) then
                        on = A(i,j)
                else 
                        on = A(i+1, j)
                end if       
                if (j==0) then
                        on = bar
                end if
        end function 

        real(8) function against(A, bar, NI, NJ, i, j)
                Implicit none
                integer NI, NJ, i, j
                real(8) bar
                real(8), dimension(0:NI, 0:NJ) :: A


                if (bar>= 0) then
                        against = A(i+1,j)
                else 
                        against = A(i, j)
                end if       
                
        end function 

        REAL(8) FUNCTION Half(arg, minus_res, plus_res)
        IMPLICIT NONE
        REAL(8) :: arg, minus_res, plus_res

        IF (arg < 0D0) THEN
                Half = minus_res
        ELSE
                Half = plus_res
        END IF

        END FUNCTION