c----------------------------------------------------------------------
c
c Routines associated with the source functions
c
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c smooth (C2) transition function from x=0 at x=x1, to 1 at x=x2
c e.g. trans=0 everywhere for x1=1.0,x2=1.0, 
c      trans=1 everywhere for x1=0.0,x2=0.0
c----------------------------------------------------------------------
        real*8 function trans(x,x1,x2)
        implicit none
        real*8 x,x1,x2

        real*8 xb

        ! initialize fixed-size variables
        data xb/0.0/

        !--------------------------------------------------------------

        xb=(x2-x)/(x2-x1)

        if (x.ge.x2) then
          trans=1
        else if (x.ge.x1) then
          trans=1-xb**3*(6*xb**2-15*xb+10)
        else
          trans=0
        end if

        return
        end

c----------------------------------------------------------------------
c Evolution of Hb is split into two routines, hb_t_evo and hb_i_evo.
c
c parameters:
c
c gauge : integer controlling which scheme
c
c x1,x2,x3,x4,xi1,xi2: real constants, determine the precise
c       implementation of each gauge
c
c current schemes:
c
c Hb_t:
c
c gauge = 0 : fixed gauge
c gauge 1 and 2 use the prescription of 2011.12970 in Cartesian Kerr-Schild coordinates.
c We call the result of this prescription "target gauge".
c Later: add option for expressions obtained from prescription of 2011.12970 
c in Cartesian Gullstrand-Painleve-like coordinates.
c
c gauge = 1 : set target gauge near AdS4D bdy of perturbed Schwarzschild-AdS4D
c             on slices at fixed x away from AdS5D bdy at x=1, then smoothly connect this
c             with target gauge near AdS5D-bdy of perturbed uniform black string-AdS5D
c gauge = 2 : sets target gauge for perturbed uniform black string-AdS5D outside 
c             of a elliptic region that includes the AdS5D boundary and the AdS4D boundary
c
c Hb_i:
c
c gauge = 0 : fixed gauge
c gauge 1 and 2 use the prescription of 2011.12970 in Cartesian Kerr-Schild coordinates.
c We call the result of this prescription "target gauge".
c Later: add option for expressions obtained from prescription of 2011.12970 
c in Cartesian Gullstrand-Painleve-like coordinates.
c
c gauge = 1 : set target gauge near AdS4D bdy of perturbed Schwarzschild-AdS4D
c             on slices at fixed x away from AdS5D bdy at x=1, then smoothly connect this along x
c             with target gauge near AdS5D-bdy of perturbed uniform black string-AdS5D
c gauge = 2 : sets target gauge for perturbed uniform black string-AdS5D outside 
c             of a elliptic region that includes the AdS5D boundary and the AdS4D boundary
c
c----------------------------------------------------------------------
        subroutine hb_t_evo(res,
     &                      gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                      gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                      gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                      gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                      gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                      gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                      gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                      Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                      Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                      Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                      L,x,y,dt,chr,ex,
     &                      phys_bdy,ghost_width,Nx,Ny,app_dim,
     &                      Hb_t_0,Hb_x_0,Hb_y_0,
     &                      gauge,t_n,x1,x2,x3,x4,y3,y4,xi1,xi2,
     &                      cbulk,unibs_background)

        implicit none

        integer unibs_background
        integer Nx,Ny,app_dim,gauge
        integer phys_bdy(2*app_dim),ghost_width(2*app_dim)
        real*8 res(Nx,Ny),t_n,t_np1
        real*8 Hb_t_0(Nx,Ny),Hb_x_0(Nx,Ny),Hb_y_0(Nx,Ny)
        real*8 Hb_t_np1(Nx,Ny),Hb_x_np1(Nx,Ny),Hb_y_np1(Nx,Ny)
        real*8 Hb_t_n(Nx,Ny),Hb_x_n(Nx,Ny),Hb_y_n(Nx,Ny)
        real*8 Hb_t_nm1(Nx,Ny),Hb_x_nm1(Nx,Ny),Hb_y_nm1(Nx,Ny)
        real*8 chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny),dt
        real*8 x1,x2,x3,x4,y3,y4,xi1,xi2,cbulk
        real*8 gb_tt_np1(Nx,Ny),gb_tx_np1(Nx,Ny),gb_ty_np1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xy_np1(Nx,Ny),gb_yy_np1(Nx,Ny)
        real*8 gb_zz_np1(Nx,Ny)
        real*8 gb_tt_n(Nx,Ny),gb_tx_n(Nx,Ny),gb_ty_n(Nx,Ny)
        real*8 gb_xx_n(Nx,Ny),gb_xy_n(Nx,Ny),gb_yy_n(Nx,Ny)
        real*8 gb_zz_n(Nx,Ny)
        real*8 gb_tt_nm1(Nx,Ny),gb_tx_nm1(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_nm1(Nx,Ny),gb_xy_nm1(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 gb_zz_nm1(Nx,Ny)

        real*8 Hb_t0,Hb_x0,Hb_y0

        real*8 F_t_np1,F_x_np1,F_y_np1
        real*8 F_t_np1_xfixslices
        real*8 F_x_np1_xfixslices
        real*8 F_y_np1_xfixslices

        real*8 x0,y0,rho0,rho1,rho2,rho3,rho4

        integer i,j,k,i1,j1
        real*8 dx,dy
        real*8 f0,f1,f1x,f1y,g0,trans

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j,i1,j1/0,0,0,0/

        data Hb_t0,Hb_x0,Hb_y0/0.0,0.0,0.0/
        data F_t_np1,F_x_np1,F_y_np1/0.0,0.0,0.0/
        data F_t_np1_xfixslices/0.0/
        data F_x_np1_xfixslices/0.0/
        data F_y_np1_xfixslices/0.0/
        data x0,y0/0.0,0.0/
        data f0,f1,f1x,f1y,g0/0.0,0.0,0.0,0.0,0.0/
        data dx,dy/0.0,0.0/

        !--------------------------------------------------------------
 
        dx=x(2)-x(1)
        dy=y(2)-y(1)

        t_np1=t_n+dt

        ! initialize output variables
        do i=1,Nx
          do j=1,Ny
            res(i,j)=0
            if (chr(i,j).eq.ex) then
              Hb_t_np1(i,j)=0
            end if
          end do
        end do

        !--------------------------------------------------------------
        ! gauge 0
        !--------------------------------------------------------------

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 1
        !
        ! set target gauge near AdS4D bdy of perturbed Schwarzschild-AdS4D
        ! on slices at fixed x away from AdS5D bdy at x=1, then smoothly connect this along x
        ! with target gauge near AdS5D-bdy of perturbed uniform black string-AdS5D
        !--------------------------------------------------------------

        if (gauge.eq.1) then
          do i=2,Nx-1
            do j=2,Ny-1
              if (chr(i,j).ne.ex) then
                x0=x(i)
                y0=y(j)

                Hb_t0=Hb_t_0(i,j)

                !define spatial transition functions along y
                f1y=trans(y0,y3,y4)

                if (unibs_background.eq.1) then
                  !stable gauge for perturbations of Schwarzschild-AdS4 
                  !in Cartesian Kerr-Schild coordinates
                  F_t_np1_xfixslices=(3.0d0/2.0d0)*gb_ty_np1(i,j)*f1y

                  !define spatial transition functions along x
                  !smoothly join target gauge near the boundary x=1 with the 
                  !F_t_np1_xfixslices in the bulk
                  f1x=trans(x0,x3,x4)
                  F_t_np1=2*gb_tx_np1(i,j)*f1x
     &                 +2*gb_tx_np1(i,j)*F_t_np1_xfixslices*(1-f1x)
                else if (unibs_background.eq.2) then
                  write(*,*) "Schwarzschild slices in Gullstrand-Painleve coords 
     &                          not yet implemented"
                  stop
                end if

                !--------------------------------------------------------------------------
                ! define time transition function
                ! used to smoothly join the t=0 gauge near-boundary value to the
                ! target, stable gauge near-boundary value 
                ! (xi2 must typically be very small, compared to xi1, 
                ! so that the target gauge is reached quickly)
                !--------------------------------------------------------------------------
                f0=trans(x0,x1,x2)
                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4

                Hb_t_np1(i,j)=F_t_np1+(Hb_t0-F_t_np1)*exp(-g0)

              end if
            end do
          end do
          return
        end if

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 2
        !
        ! sets target gauge for perturbed uniform black string-AdS5D outside 
        ! of a elliptic region that includes the AdS5D boundary and the AdS4D boundary
        !--------------------------------------------------------------

        if (gauge.eq.2) then
          do i=2,Nx-1
            do j=2,Ny-1
              if (chr(i,j).ne.ex) then
                x0=x(i)
                y0=y(j)
                rho0=sqrt(x0*x0+y0*y0)
                rho1=sqrt(x1*x1+y3*y3)
                rho2=sqrt(x2*x2+y4*y4)
                rho3=sqrt(x3*x3+y3*y3)
                rho4=sqrt(x4*x4+y4*y4)

                Hb_t0=Hb_t_0(i,j)

                !define spatial transition functions along y
                f1=trans(rho0,rho3,rho4)

                if (unibs_background.eq.1) then
                  F_t_np1=2*gb_tx_np1(i,j)*f1
     &                 +2*gb_tx_np1(i,j)*cbulk*(1-f1)
                else if (unibs_background.eq.2) then
                  write(*,*) "Schwarzschild slices in Gullstrand-Painleve coords 
     &                          not yet implemented"
                  stop
                end if

                !--------------------------------------------------------------------------
                ! define time transition function
                ! used to smoothly join the t=0 gauge near-boundary value to the
                ! target gauge near-boundary value 
                ! (xi2 must typically be very small, compared to xi1, 
                ! so that the target gauge is reached quickly)
                !--------------------------------------------------------------------------
                f0=trans(rho0,rho1,rho2)
                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4

                Hb_t_np1(i,j)=F_t_np1+(Hb_t0-F_t_np1)*exp(-g0)

              end if
            end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! otherwise
        !--------------------------------------------------------------

        write(*,*) 'hb_t_evo : error, gauge,',gauge,' unknown'

        return
        end

c-----------------------------------------------------------------------
        subroutine hb_i_evo(res,
     &                      gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                      gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                      gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                      gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                      gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                      gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                      gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                      Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                      Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                      Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                      L,x,y,dt,chr,ex,
     &                      phys_bdy,ghost_width,Nx,Ny,app_dim,
     &                      Hb_t_0,Hb_x_0,Hb_y_0,
     &                      gauge,t_n,x1,x2,x3,x4,y3,y4,xi1,xi2,
     &                      cbulk,unibs_background)

        implicit none

        integer unibs_background
        integer Nx,Ny,app_dim,gauge
        integer phys_bdy(2*app_dim),ghost_width(2*app_dim)
        real*8 res(Nx,Ny),t_n,t_np1
        real*8 Hb_t_0(Nx,Ny),Hb_x_0(Nx,Ny),Hb_y_0(Nx,Ny)
        real*8 Hb_t_np1(Nx,Ny),Hb_x_np1(Nx,Ny),Hb_y_np1(Nx,Ny)
        real*8 Hb_t_n(Nx,Ny),Hb_x_n(Nx,Ny),Hb_y_n(Nx,Ny)
        real*8 Hb_t_nm1(Nx,Ny),Hb_x_nm1(Nx,Ny),Hb_y_nm1(Nx,Ny)
        real*8 chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny),dt
        real*8 x1,x2,x3,x4,y3,y4,xi1,xi2,cbulk
        real*8 gb_tt_np1(Nx,Ny),gb_tx_np1(Nx,Ny),gb_ty_np1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xy_np1(Nx,Ny),gb_yy_np1(Nx,Ny)
        real*8 gb_zz_np1(Nx,Ny)
        real*8 gb_tt_n(Nx,Ny),gb_tx_n(Nx,Ny),gb_ty_n(Nx,Ny)
        real*8 gb_xx_n(Nx,Ny),gb_xy_n(Nx,Ny),gb_yy_n(Nx,Ny)
        real*8 gb_zz_n(Nx,Ny)
        real*8 gb_tt_nm1(Nx,Ny),gb_tx_nm1(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_nm1(Nx,Ny),gb_xy_nm1(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 gb_zz_nm1(Nx,Ny)

        real*8 Hb_t0,Hb_x0,Hb_y0

        real*8 F_t_np1,F_x_np1,F_y_np1
        real*8 F_t_np1_xfixslices
        real*8 F_x_np1_xfixslices
        real*8 F_y_np1_xfixslices


        real*8 x0,y0,rho0,rho1,rho2,rho3,rho4

        integer i,j,k,i1,j1
        real*8 dx,dy
        real*8 f0,f1,f1x,f1y,g0,trans

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j,i1,j1/0,0,0,0/

        data Hb_t0,Hb_x0,Hb_y0/0.0,0.0,0.0/
        data F_t_np1,F_x_np1,F_y_np1/0.0,0.0,0.0/
        data F_t_np1_xfixslices/0.0/
        data F_x_np1_xfixslices/0.0/
        data F_y_np1_xfixslices/0.0/
        data x0,y0/0.0,0.0/
        data f0,f1,f1x,f1y,g0/0.0,0.0,0.0,0.0,0.0/
        data dx,dy/0.0,0.0/

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        t_np1=t_n+dt

        ! initialize output variables
        do i=1,Nx
          do j=1,Ny
            res(i,j)=0
            if (chr(i,j).eq.ex) then
              Hb_x_np1(i,j)=0
              Hb_y_np1(i,j)=0
            end if
          end do
        end do

        !--------------------------------------------------------------
        ! gauge 0
        !--------------------------------------------------------------

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 1
        !
        ! set target gauge near AdS4D bdy of perturbed Schwarzschild-AdS4D
        ! on slices at fixed x away from AdS5D bdy at x=1, then smoothly connect this along x
        ! with target gauge near AdS5D-bdy of perturbed uniform black string-AdS5D
        !--------------------------------------------------------------

        if (gauge.eq.1) then
          do i=2,Nx-1
            do j=2,Ny-1
              if (chr(i,j).ne.ex) then
                x0=x(i)
                y0=y(j)

                Hb_x0=Hb_x_0(i,j)
                Hb_y0=Hb_y_0(i,j)


                !define spatial transition functions along y
                f1y=trans(y0,y3,y4)

                if (unibs_background.eq.1) then
                  !stable gauge for perturbations of Schwarzschild-AdS4 
                  !in Cartesian Kerr-Schild coordinates
                  F_x_np1_xfixslices=0
                  F_y_np1_xfixslices=(3.0d0/2.0d0)*gb_yy_np1(i,j)*f1y

                  !define spatial transition functions along x
                  !smoothly join target gauge near the boundary x=1 with the 
                  !F_i_np1_xfixslices in the bulk
                  f1x=trans(x0,x3,x4)
                  F_x_np1=((1.0d0/2.0d0)*(1-y0**2)**2*gb_yy_np1(i,j)
     &               +2*gb_xx_np1(i,j))*f1x
     &             +((1.0d0/2.0d0)*(1-y0**2)**2*gb_yy_np1(i,j)
     &               +2*gb_xx_np1(i,j))*F_x_np1_xfixslices*(1-f1x)
                  F_y_np1=2*gb_xy_np1(i,j)*f1x
     &               +2*gb_xy_np1(i,j)*F_y_np1_xfixslices*(1-f1x)
                else if (unibs_background.eq.2) then
                  write(*,*) "Schwarzschild slices in Gullstrand-Painleve coords 
     &                          not yet implemented"
                  stop
                end if

                !--------------------------------------------------------------------------
                ! define time transition function
                ! used to smoothly join the t=0 gauge near-boundary value to the
                ! target, stable gauge near-boundary value 
                ! (xi2 must typically be very small, compared to xi1, 
                ! so that the target gauge is reached quickly)
                !--------------------------------------------------------------------------
                f0=trans(x0,x1,x2)
                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4

                Hb_x_np1(i,j)=F_x_np1+(Hb_x0-F_x_np1)*exp(-g0)
                Hb_y_np1(i,j)=F_y_np1+(Hb_x0-F_y_np1)*exp(-g0)

              end if
            end do
          end do
          return
        end if


        !--------------------------------------------------------------
        ! gauge 2
        !
        ! sets target gauge for perturbed uniform black string-AdS5D outside 
        ! of a elliptic region that includes the AdS5D boundary and the AdS4D boundary
        !--------------------------------------------------------------

        if (gauge.eq.2) then
          do i=2,Nx-1
            do j=2,Ny-1
              if (chr(i,j).ne.ex) then
                x0=x(i)
                y0=y(j)
                rho0=sqrt(x0*x0+y0*y0)
                rho1=sqrt(x1*x1+y3*y3)
                rho2=sqrt(x2*x2+y4*y4)
                rho3=sqrt(x3*x3+y3*y3)
                rho4=sqrt(x4*x4+y4*y4)

                Hb_x0=Hb_x_0(i,j)
                Hb_y0=Hb_y_0(i,j)


                !define spatial transition functions along y
                f1=trans(rho0,rho3,rho4)

                if (unibs_background.eq.1) then
                  F_x_np1=((1.0d0/2.0d0)*(1-y0**2)**2*gb_yy_np1(i,j)
     &                        +2*gb_xx_np1(i,j))*f1
     &                   +((1.0d0/2.0d0)*(1-y0**2)**2*gb_yy_np1(i,j)
     &                        +2*gb_xx_np1(i,j))*cbulk*(1-f1)
                  F_y_np1=2*gb_xy_np1(i,j)*f1
     &                 +2*gb_xy_np1(i,j)*cbulk*(1-f1)
                else if (unibs_background.eq.2) then
                  write(*,*) "Schwarzschild slices in Gullstrand-Painleve coords 
     &                          not yet implemented"
                  stop
                end if

                !--------------------------------------------------------------------------
                ! define time transition function
                ! used to smoothly join the t=0 gauge near-boundary value to the
                ! target gauge near-boundary value 
                ! (xi2 must typically be very small, compared to xi1, 
                ! so that the target gauge is reached quickly)
                !--------------------------------------------------------------------------
                f0=trans(rho0,rho1,rho2)
                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4

                Hb_x_np1(i,j)=F_x_np1+(Hb_x0-F_x_np1)*exp(-g0)
                Hb_y_np1(i,j)=F_y_np1+(Hb_y0-F_y_np1)*exp(-g0)

              end if
            end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! otherwise
        !--------------------------------------------------------------

        write(*,*) 'hb_i_evo : error, gauge,',gauge,' unknown'

        return
        end

c-----------------------------------------------------------------------
