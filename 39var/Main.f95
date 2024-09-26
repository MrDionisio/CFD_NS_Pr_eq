        Program Pr
        Implicit none
 
        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        INTEGER NI, NJ, NITER, SMAX
        INTEGER I,J
        REAL(8) L,H,U0,MU,R0,P0, NU0, k, gamma, Re_x
        REAL(8) dx,dy,CFL,EPS,dt, C, C_f
        REAL(8),ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL(8),ALLOCATABLE :: U(:,:),V(:,:),P(:,:),R(:,:), NU(:,:), RHO(:,:)


        write(*,*) 'Read input file'
        open(IO,FILE='Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) EPS
        read(IO,*) SMAX

        read(IO,*) U0
        read(IO,*) P0
        read(IO,*) CFL
        read(IO,*) NU0
        read(IO,*) k
        read(IO, *) R0
        read(IO,*) gamma
        CLOSE(IO)
   
        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates

!----------------- Node variables -----------------------------
        allocate(U(0:NI,0:NJ))  
        allocate(V(0:NI,0:NJ))  
        allocate(P(0:NI,0:NJ))
        allocate(NU(0:NI,0:NJ))
        allocate(RHO(0:NI,0:NJ))

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

        

        write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
        write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
        write(*,*)'ReL= ', U0*L/NU0, 'CFL = ', CFL, "dt = ", dt
      

!----------------- Initial fields -----------------------------

        call nu_turbulent(NU, NU0, k, dy, H, NI, NJ)
        RHO = R0
        U=U0
        V=1.0e-5
        P=P0
        C=P0/(R0**gamma)
        print*, NU(1,1), RHO(1,1)

!---------------- Solve Prandtl equations ---------------------
 
        write(*,*) 'Solve Novie-Stocks equations'
        Open(IO,FILE='Results0.plt')
        Call OUTPUT(IO,NI,NJ,X_Node,Y_Node,U,V,P,NU, RHO)
        Close(IO)    

        call NOVIE_STOCKS(U, V, P, dx, dy, dt, EPS, NI, NJ, NU, U0, SMAX, RHO, gamma, c, R0, P0)
        MU = NU0/R0
        Open(IO,FILE='C_f_Re.txt')
        do i=1,NI
                C_f=MU*(U(i,1)-U(i,0))/dy/R0/U0/U0*2
                write(IO,*) U0*(dx*(i-0.5))/NU0, C_f
        end do
        Close(IO)

        Open(IO,FILE='C_f3.txt')
                do i=1,NI
                        Re_x = U0*dx*(i-0.5)/NU0
                        C_f=0.664/sqrt(Re_x)
                        write(IO,*) Re_x , C_f
                end do
        Close(IO)

        Open(IO,FILE='11.txt')
                do j=0,NJ
                        write(IO,*) U(NI-1,j), (j-0.5)*dy
                end do
        Close(IO)
 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data' 
        Open(IO,FILE='Results.plt')
        Call OUTPUT(IO,NI,NJ,X_Node,Y_Node,U,V,P,NU, RHO)
        Close(IO)
        DEALLOCATE(U,V,P,NU,RHO,X_Node,Y_Node)
        
        END PROGRAM

        SUBROUTINE OUTPUT(IO,NI,NJ,X,Y,U,V,P, NU, RHO)
                IMPLICIT NONE
                INTEGER NI,NJ,IO
                REAL(8),DIMENSION(NI,NJ):: X,Y
                REAL(8),DIMENSION(0:NI,0:NJ):: U,V,P,NU, RHO

                Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P", "NU" , "RHO"' 
                Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
                Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
                Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
                Write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
                Write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
                Write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)
                Write(IO,'(100E25.16)') NU(1:NI-1,1:NJ-1)
                Write(IO,'(100E25.16)') RHO(1:NI-1,1:NJ-1)
                CLOSE(1)
        END SUBROUTINE

        SUBROUTINE NOVIE_STOCKS(U, V, P, dx, dy, dt, EPS, NI, NJ, NU, U0, SMAX, RHO, gamma, c, R0, P0)
                Implicit none
                real(8) EPS, dx, dy, dt, U0, ERR_U, ERR_V, ERR_P, ERR_S, gamma, c, R0, P0
                real(8) Half
                INTEGER NI, NJ, i, j, SMAX, S
                real(8), dimension(0:NI, 0:NJ) :: U, V, P, UN, VN, PN, R_U, R_V, R_P, NU, RHO
                real(8), dimension (6) :: left, right, up, down ! 1 - U; 2 - V; 3 - P
                real(8), dimension (2) :: bar_l, bar_r, bar_u, bar_d ! 1 - U; 2 - V

                ERR_U = 0
                ERR_V = 0
                ERR_P = EPS+1
                S = 0 

                R_U=0.0
                R_V=0.0
                R_P=0.0
                

                OPEN(1, file="R_U.txt")
                OPEN(2, file="R_V.txt")
                OPEN(3, file="R_P.txt")

                DO WHILE ((MAX(ERR_U, MAX(ERR_V, ERR_P))>EPS).AND.(S<SMAX))

                        U(: , 0) = - U(:, 1)
                        V(: , 0) = - V(:, 1)
                        P(: , 0) = P(:, 1)
                        

                        U(NI, :) = U(NI-1, :)
                        V(NI, :) = V(NI-1, :)
                        P(NI, :) = P0
                        

                        U(0, :) = U0
                        V(0, :) = 0.0
                        P(0, :) = P(1, :)
                        

                        U(:, NJ) = U(:, NJ-1)
                        V(:, NJ) = - V(:, NJ-1)
                        P(:, NJ) = P(:, NJ-1)

                        
                        
                        do i=1,NI-1
                                do j=1,NJ-1

                                        call bar(U, V, bar_d, bar_l, bar_r, bar_u, NI, NJ, i, j, RHO)
                                        
                                        left(1) = Half(bar_l(1), U(i,j), U(i-1,j))
                                        left(2) = Half(bar_l(1), V(i,j), V(i-1,j))
                                        left(3) = Half(bar_l(1), P(i-1,j), P(i,j))
                                        left(4) = 0.5 * (NU(i,j) + NU(i-1,j))
                                        left(5) = 0.5 * (RHO(i,j) + RHO(i-1,j))
                                        left(6) = left(4)*left(5)

                                        right(1) = Half(bar_r(1), U(i+1,j), U(i,j))
                                        right(2) = Half(bar_r(1), V(i+1,j), V(i,j))
                                        right(3) = Half(bar_r(1), P(i,j), P(i+1,j))
                                        right(4) = 0.5 * (NU(i+1,j) + NU(i,j))
                                        right(5) = 0.5 * (RHO(i+1,j) + RHO(i,j))
                                        right(6) = right(5)*right(4)

                                        up(1) = Half(bar_u(2), U(i,j+1), U(i,j))
                                        up(2) = Half(bar_u(2), V(i,j+1), V(i,j))
                                        up(3) = Half(bar_u(2), P(i,j), P(i,j+1))
                                        up(4) = 0.5 * (NU(i,j+1) + NU(i,j))
                                        up(5) = 0.5 * (RHO(i,j+1) + RHO(i,j))
                                        up(6) = up(5)*up(4)

                                        down(1) = Half(bar_d(2), U(i,j), U(i,j-1))
                                        down(2) = Half(bar_d(2), V(i,j), V(i,j-1))
                                        down(3) = Half(bar_d(2), P(i,j-1), P(i,j))
                                        down(4) = 0.5 * (NU(i,j) + NU(i,j-1))
                                        down(5) = 0.5 * (RHO(i,j) + RHO(i,j-1))
                                        down(6) = down(5)*down(4)

                                        if (j==1) THEN
                                                R_P(i,j) = (bar_r(1)-bar_l(1))/dx+bar_u(2)/dy
                                        else
                                                R_P(i,j) = (bar_r(1)-bar_l(1))/dx+(bar_u(2)-bar_d(2))/dy
                                        end if
                                        
                                        
                                        
                                        R_U(i,j) = -1*(bar_r(1)*right(1)-bar_l(1)*left(1))/dx-(bar_u(2)*up(1)-bar_d(2)*down(1))/dy &
                                        -(right(3)-left(3))/dx + 2*(right(6)*(U(i+1,j)-U(i,j))/dx-left(6)*(U(i,j)-U(i-1,j))/dx)/dx &
                                        + (up(6)*(V(i+1, j+1)-V(i-1,j+1))/2/dx-down(6)*(V(i+1,j-1)-V(i-1,j-1))/2/dx)/2/dy &
                                        + (up(6)*(U(i,j+1)-U(i,j))/dy-down(6)*(U(i,j)-U(i,j-1))/dy)/dy &
                                        - 2/3 * (right(6)*(U(i+1,j)-U(i,j))/dx-left(6)*(U(i,j)-U(i-1,j))/dx)/dx &
                                        - 2/3 * (right(6)*(V(i+1, j+1)-V(i+1,j-1))/2/dy-left(6)*(V(i-1,j+1)-V(i-1,j-1))/2/dy)/2/dx
                                        
                                        R_V(i,j) = -1*(bar_r(1)*right(2)-bar_l(1)*left(2))/dx-(bar_u(2)*up(2)-bar_d(2)*down(2))/dy &
                                        - (up(3)-down(3))/dy+(right(6)*(V(i+1,j)-V(i,j))/dx-left(6)*(V(i,j)-V(i-1,j))/dx)/dx &
                                        + (right(6)*(U(i+1,j+1)-U(i+1,j-1))/2/dy-left(6)*(U(i-1,j+1)-U(i-1,j-1))/2/dy)/2/dx &
                                        + 2 * (up(6)*(V(i,j+1)-V(i,j))/dy-down(6)*(V(i,j)-V(i,j-1))/dy)/dy &
                                        - 2/3 * (up(6)*(U(i+1,j+1)-U(i-1,j+1))/2/dx-down(6)*(U(i+1,j-1)-U(i-1,j-1))/2/dx)/2/dy &
                                        - 2/3 * (up(6)*(V(i,j+1)-V(i,j))/dy-down(6)*(V(i,j)-V(i,j-1))/dy)/dy

                                end do
                        end do
                        
                        P(:,:) = P(:,:) - dt/U0**2*R_P(:,:)
                        U(:,:) = (RHO(:,:)*U(:,:) + dt*R_U(:,:))/(P(:,:)/C)**(1/gamma)
                        V(:,:) = (RHO(:,:)*V(:,:) + dt*R_V(:,:))/(P(:,:)/C)**(1/gamma)
                        

                        RHO(:,:) = (P(:,:)/C)**(1/gamma)
                       
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

        !Вывод невязок
        SUBROUTINE Residuals(IO1, IO2, IO3, S, ERR1, ERR2, ERR3)
                integer IO1, IO2, IO3, S
                REAL(8) ERR1, ERR2, ERR3
                Write(IO1, *) S, ERR1
                Write(IO2, *) S, ERR2
                Write(IO3, *) S, ERR3
                Print*, "S", S, ERR1, ERR2, ERR3
        END SUBROUTINE

        !Вычисление скоростей с шапкой
        SUBROUTINE bar(U, V, bar_d, bar_l, bar_r, bar_u, NI, NJ, i, j, RHO)
                integer NI, NJ, i, j
                real(8), dimension(0:NI, 0:NJ) :: U, V, RHO
                real(8), dimension(2) ::  bar_d, bar_l, bar_r, bar_u

                bar_l(1) = 0.5*(RHO(i,j)*U(i,j)+RHO(i-1,j)*U(i-1,j))
                bar_l(2) = 0.5*(RHO(i,j)*V(i,j)+RHO(i-1,j)*V(i-1,j))
                bar_r(1) = 0.5*(RHO(i,j)*U(i,j)+RHO(i+1,j)*U(i+1,j))
                bar_r(2) = 0.5*(RHO(i,j)*V(i,j)+RHO(i+1,j)*U(i+1,j))
                bar_d(1) = 0.5*(RHO(i,j)*U(i,j)+RHO(i,j-1)*U(i,j-1))
                bar_d(2) = 0.5*(RHO(i,j)*V(i,j)+RHO(i,j-1)*V(i,j-1))
                bar_u(1) = 0.5*(RHO(i,j)*U(i,j)+RHO(i,j+1)*V(i,j+1))
                bar_u(2) = 0.5*(RHO(i,j)*V(i,j)+RHO(i,j+1)*V(i,j+1))
        END SUBROUTINE

        !Поточная и противопоточная схема
        REAL(8) FUNCTION Half(arg, minus_res, plus_res)
        IMPLICIT NONE
        REAL(8) :: arg, minus_res, plus_res

        IF (arg < 0D0) THEN
                Half = minus_res
        ELSE
                Half = plus_res
        END IF

        END FUNCTION
        !Вычисление поля турбулентной вязкости
        SUBROUTINE nu_turbulent (NU, NU0, k , dy, H, NI, NJ)
                Implicit none 
                real(8), dimension(0:NI,0:NJ) ::  NU
                real(8) NU0, k, dy, H
                integer NI, NJ, i, j

                        do j=1, NJ-1
                                NU(:,j) = NU0*(1+k*(H-dy*(j-0.5)))
                        end do

                NU(:, 0) = NU(:, 1)
                NU(:, NJ) = NU(:, NJ-1)
                open(1,file="NU.txt")
                write(1,*) NU(:, :)
                close(1)
        end SUBROUTINE