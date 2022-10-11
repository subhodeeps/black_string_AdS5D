c----------------------------------------------------------------------
c in cartesian coordinates t,x,y,z for x,y,z in [-1,1]
c
c routines associated with the source functions
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c smooth (C2) transition function from 0 at x=rho1, to 1 at x=rho2
c i.e. trans=0 everywhere for rho1=1.0,rho2=1.0, 
c      trans=1 everywhere for rho1=0.0,rho2=0.0
c----------------------------------------------------------------------
        real*8 function trans(x,rho1,rho2)
        implicit none
        real*8 x,rho1,rho2

        real*8 xb

        ! initialize fixed-size variables
        data xb/0.0/

        !--------------------------------------------------------------

        xb=(rho2-x)/(rho2-rho1)

        if (x.ge.rho2) then
          trans=1
        else if (x.ge.rho1) then
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
c rho1,rho2,xi1,xi2: real constants, gauge dependent
c
c current schemes:
c
c Hb_t:
c
c gauge = 0 : fixed gauge
c gauge = 1 : AdS-asymptotics gauge
c
c Hb_i:
c
c gauge = 0 : fixed gauge
c gauge = 1 : AdS-asymptotics gauge
c
c----------------------------------------------------------------------
        subroutine hb_t_evo(res,
     &                      gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                      gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                      gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                      gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                      gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                      gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                      gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                      gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                      gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                      gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                      Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                      Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                      Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                      Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                      phi1_np1,phi1_n,phi1_nm1,
     &                      L,x,y,z,dt,chr,ex,
     &                      phys_bdy,ghost_width,Nx,Ny,Nz,
     &                      Hb_t_0,Hb_x_0,Hb_y_0,
     &                      Hb_z_0,
     &                      gauge,t_n,rho1,rho2,rho3,rho4,xi1,xi2,
     &                      c1,c2,c3,cbulk)
        implicit none
        integer Nx,Ny,Nz,gauge
        integer phys_bdy(6),ghost_width(6)
        real*8 c1,c2,c3
        real*8 res(Nx,Ny,Nz),t_n,t_np1
        real*8 Hb_t_0(Nx,Ny,Nz),Hb_x_0(Nx,Ny,Nz),Hb_y_0(Nx,Ny,Nz)
        real*8 Hb_z_0(Nx,Ny,Nz)
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_x_np1(Nx,Ny,Nz),Hb_y_np1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz)
        real*8 Hb_t_n(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz)
        real*8 Hb_z_n(Nx,Ny,Nz)
        real*8 Hb_t_nm1(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_nm1(Nx,Ny,Nz)
        real*8 chr(Nx,Ny,Nz),ex,L
        real*8 x(Nx),y(Ny),z(Nz),dt,rho1,rho2,rho3,rho4,xi1,xi2,cbulk
        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tx_np1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_xx_np1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_yy_np1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz)
        real*8 gb_tt_n(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz)
        real*8 gb_ty_n(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz)
        real*8 gb_tz_n(Nx,Ny,Nz)
        real*8 gb_xz_n(Nx,Ny,Nz)
        real*8 gb_yz_n(Nx,Ny,Nz)
        real*8 gb_yy_n(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz)
        real*8 gb_tt_nm1(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_nm1(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xy_nm1(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_zz_nm1(Nx,Ny,Nz)
        real*8 gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yz_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)

        real*8 Hb_t0,Hb_x0,Hb_y0
        real*8 Hb_z0

        real*8 F_t_np1,F_x_np1,F_y_np1
        real*8 F_z_np1
        real*8 G_t_np1,G_x_np1,G_y_np1
        real*8 G_z_np1

        real*8 g0uu11,gbtt,gbtx,gbty,gbtz
        real*8 gbxx,gbxy,gbxz,gbyy,gbyz,gbzz
        real*8 alphasq0,alpha0,betax0,betay0,betaz0

        real*8 x0,y0,z0,rho0

        integer i,j,k,i1,j1
        real*8 dx,dy,dz
        real*8 f0,f1,f2,g0,g1,g2,trans

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j,k,i1,j1/0,0,0,0,0/

        data Hb_t0,Hb_x0,Hb_y0/0.0,0.0,0.0/
        data Hb_z0/0.0/
        data F_t_np1,F_x_np1,F_y_np1/0.0,0.0,0.0/
        data F_z_np1/0.0/
        data G_t_np1,G_x_np1,G_y_np1/0.0,0.0,0.0/
        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/
        data f0,f1,f2,g0,g1,g2/0.0,0.0,0.0,0.0,0.0,0.0/
        data dx,dy,dz/0.0,0.0,0.0/

        !--------------------------------------------------------------
 
        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        t_np1=t_n+dt

        ! initialize output variables
        do i=1,Nx
          do j=1,Ny
           do k=1,Nz
            res(i,j,k)=0
            if (chr(i,j,k).eq.ex) then
              Hb_t_np1(i,j,k)=0
            end if
           end do
          end do
        end do

        !--------------------------------------------------------------
        ! gauge 0
        !--------------------------------------------------------------

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 1 
        !--------------------------------------------------------------

        if (gauge.eq.1) then
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              if (chr(i,j,k).ne.ex) then
                x0=x(i)
                y0=y(j)
                rho0=sqrt(x0**2+y0**2)

                F_t_np1=
     &                 +gb_tx_np1(i,j,k)*x0/rho0*(200)/100
     &                 +gb_ty_np1(i,j,k)*y0/rho0*(200)/100

                Hb_t_np1(i,j,k)=F_t_np1
              end if
             end do
            end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! COMMENTS FROM PREVIOUS CODES
        ! gauge 2 (works for stabilizing gb_tt,gb_tx,gb_ty,psi)
        !         (rho_cd=-1 stabilizes left corner of gb_tt,gb_tx,gb_ty,psi)
        !         (rho_cd=0  stabilizes gb_yy)
        !         (rho_cd=0  stabilizes gb_xx asymmetry)
        !--------------------------------------------------------------

        if (gauge.eq.2) then
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              if (chr(i,j,k).ne.ex) then
                x0=x(i)
                y0=y(j)
                rho0=sqrt(x0**2+y0**2)

                F_t_np1=
     &                 +gb_tx_np1(i,j,k)*2*x0/rho0
     &                 +gb_ty_np1(i,j,k)*2*y0/rho0

                Hb_t_np1(i,j,k)=F_t_np1
              end if
             end do
            end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! gave the stable runs used in 2011.12970
        ! gauge 3
        !--------------------------------------------------------------

        if (gauge.eq.3) then
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              if (chr(i,j,k).ne.ex) then
                x0=x(i)
                y0=y(j)
                z0=z(k)
                rho0=sqrt(x0**2+y0**2+z0**2)
!             write(*,*) "DEBUG from gauge.f"
!             write(*,*) "i,j,k,x0,y0,z0,dx=",i,j,k,x0,y0,z0,dx
!             write(*,*) "rho0=",rho0

!        write(*,*) 'DEBUG from gauge.f'
!        write(*,*) 'L,x0,y0,z0,rho0,dx=',L,x0,y0,z0,rho0,dx

                Hb_t0=Hb_t_0(i,j,k)

                f0=trans(rho0,rho1,rho2)
                f1=trans(rho0,rho3,rho4)
             !time transition with transition near the AdS boundary   
             !use when using initial data that DOES NOT preserve the Hamiltonian and momentum constraints
!                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4
             !time transition with NO transition near the AdS boundary
             ! use when using initial data that preserves the Hamiltonian and momentum constraints
                g0=(t_np1/(xi1*(1-f0)))**4
!        write(*,*) 'f0,f1,g0=',f0,f1,g0
!        write(*,*) 'F_t_np1,Hb_t_np1(i,j,k)=',F_t_np1,Hb_t_np1(i,j,k)

        if (rho0.ne.0.0d0) then        
                F_t_np1=gb_tx_np1(i,j,k)*1.5d0*x0/rho0*f1
     &                 +gb_ty_np1(i,j,k)*1.5d0*y0/rho0*f1
     &                 +gb_tz_np1(i,j,k)*1.5d0*z0/rho0*f1
     &                 +gb_tx_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_ty_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_tz_np1(i,j,k)*cbulk*z0*(1-f1)
        else
                F_t_np1=gb_tx_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_ty_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_tz_np1(i,j,k)*cbulk*z0*(1-f1)
        end if

                Hb_t_np1(i,j,k)=F_t_np1+(Hb_t0-F_t_np1)*exp(-g0)

!        write(*,*) 'DEBUG from gauge.f'
!        write(*,*) 'L,x0,y0,z0,rho0,dx=',L,x0,y0,z0,rho0,dx
!        write(*,*) 'F_t_np1,Hb_t_np1(i,j,k)=',F_t_np1,Hb_t_np1(i,j,k)
!        write(*,*) 'gb_tx_np1(i,j,k)=',gb_tx_np1(i,j,k)

              end if
            end do
           end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! gauge 4
        !--------------------------------------------------------------

        if (gauge.eq.4) then
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              if (chr(i,j,k).ne.ex) then
                x0=x(i)
                y0=y(j)
                z0=z(k)
                rho0=sqrt(x0**2+y0**2+z0**2)
                Hb_t0=Hb_t_0(i,j,k)

                f0=trans(rho0,rho1,rho2)
                f1=trans(rho0,rho3,rho4)
                f2=trans(rho0,0.0d0,c2)
                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4
                g1=(t_np1/xi1)**4
                g2=(t_np1/c3)**4

                gbtt =gb_tt_np1(i,j,k)
                gbtx =gb_tx_np1(i,j,k)
                gbty =gb_ty_np1(i,j,k)
                gbtz =gb_tz_np1(i,j,k)
                gbxx =gb_xx_np1(i,j,k)
                gbxy =gb_xy_np1(i,j,k)
                gbxz =gb_xz_np1(i,j,k)
                gbyy =gb_yy_np1(i,j,k)
                gbyz =gb_yz_np1(i,j,k)
                gbzz=gb_zz_np1(i,j,k)
                g0uu11=(-gbxz**2*(4+gbyy)-4*gbyz**2+4*gbyy*gbzz
     &               +16*(4+gbyy+gbzz)+2*(-16*gbyy+gbxz**2*(8+3*gbyy)
     &               +8*gbyz**2
     &               -8*(2+gbyy)*gbzz)*rho0**2+(-3*gbxz**2*(8+5*gbyy)
     &               +8*(-3*gbyz**2+2*gbzz+gbyy*(2+3*gbzz)))*rho0**4
     &               +4*(gbxz**2*(4+5*gbyy)
     &                +4*gbyz**2-4*gbyy*gbzz)*rho0**6
     &               -(gbxz**2*(4+15*gbyy)+4*gbyz**2
     &                -4*gbyy*gbzz)*rho0**8
     &               +6*gbxz**2*gbyy*rho0**10-gbxz**2*gbyy*rho0**12
     &               +2*gbxy*gbxz*gbyz*(-1+rho0**2)**6
     &               -gbxy**2*(-1+rho0**2)**4*(4+gbzz*(-1+rho0**2)**2)
     &               +gbxx*(-1+rho0**2)**2*(-gbyz**2*(-1+rho0**2)**4
     &               +4*(4+gbzz*(-1+rho0**2)**2)
     &               +gbyy*(-1+rho0**2)**2*(4+gbzz*(-1+rho0**2)**2)))
     &               /((-1+rho0**2)**6
     &                *(gbyz*(-gbtx*(gbtz*gbxy+gbty*gbxz)
     &                +gbtx**2*gbyz+(gbty*gbtz*(4+gbxx*(-1+rho0**2)**2))
     &               /(-1+rho0**2)**2
     &               +((gbtt*(-1+rho0**2)**2
     &               -(1+rho0**2)**2)*(gbxy*gbxz*(-1+rho0**2)**2
     &               -gbyz*(4+gbxx*(-1+rho0**2)**2)))/(-1+rho0**2)**4)
     &               +gbxz*(-gbty*gbtz*gbxy+gbty**2*gbxz-gbtx*gbty*gbyz
     &               +(gbtx*gbtz*(4+gbyy*(-1+rho0**2)**2))
     &               /(-1+rho0**2)**2
     &               +((gbtt*(-1+rho0**2)**2
     &                -(1+rho0**2)**2)*(gbxy*gbyz*(-1+rho0**2)**2
     &               -gbxz*(4+gbyy*(-1+rho0**2)**2)))/(-1+rho0**2)**4)
     &               +(1/((-1+rho0**2)**4))*gbtz*((-1+rho0**2)**2
     &               *(-gbty*gbxy*gbxz*(-1+rho0**2)**2
     &               -gbtx*gbxy*gbyz*(-1+rho0**2)**2
     &                +gbty*gbyz*(4+gbxx*(-1+rho0**2)**2)
     &               +gbtx*gbxz*(4+gbyy*(-1+rho0**2)**2))
     &               -gbtz*(-gbxy**2*(-1+rho0**2)**4
     &               +4*(4+gbyy*(-1+rho0**2)**2)
     &               +gbxx*(-1+rho0**2)**2*(4+gbyy*(-1+rho0**2)**2)))
     &               -(1/((-1+rho0**2)**8))*(4+gbzz*(-1+rho0**2)**2)
     &                *(gbty**2*(4+gbxx)
     &                -gbxy**2+gbxx*gbyy
     &               -2*(gbty**2*(8+3*gbxx)-gbxy**2+gbxx*gbyy)*rho0**2
     &               +(16+3*gbty**2*(8+5*gbxx)+gbxy**2-8*gbyy
     &                -gbxx*(8+gbyy))*rho0**4
     &               -4*(gbty**2*(4+5*gbxx)+gbxy**2-gbxx*gbyy)*rho0**6
     &               +(gbty**2*(4+15*gbxx)+gbxy**2
     &                -gbxx*(-4+gbyy)+4*gbyy)*rho0**8
     &               +2*(gbxy**2-gbxx*(3*gbty**2+gbyy))*rho0**10
     &               +(-gbxy**2+gbxx*(gbty**2+gbyy))*rho0**12
     &               -2*gbtx*gbty*gbxy*(-1+rho0**2)**6
     &               +4*(4+gbxx+gbyy+8*rho0**2)
     &               +gbtx**2*(-1+rho0**2)**4*(4+gbyy*(-1+rho0**2)**2)
     &               -gbtt*(-1+rho0**2)**2*(-gbxy**2*(-1+rho0**2)**4
     &               +4*(4+gbyy*(-1+rho0**2)**2)
     &               +gbxx*(-1+rho0**2)**2*(4+gbyy*(-1+rho0**2)**2)))))
                alphasq0=-1/g0uu11
                if (alphasq0.ge.0) then
                  alpha0=sqrt(alphasq0)
                else
                  alpha0=1e-8
                end if
                betax0=gb_tx_np1(i,j,k)
                betay0=gb_ty_np1(i,j,k)
                betaz0=gb_tz_np1(i,j,k)

        if (rho0.ne.0.0d0) then
                F_t_np1=
     &                 +gb_tx_np1(i,j,k)*1.5d0*x0/rho0*f1
     &                 +gb_ty_np1(i,j,k)*1.5d0*y0/rho0*f1
     &                 +gb_tz_np1(i,j,k)*1.5d0*z0/rho0*f1
     &                 +gb_tx_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_ty_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_tz_np1(i,j,k)*cbulk*z0*(1-f1)
         else
                F_t_np1=gb_tx_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_ty_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_tz_np1(i,j,k)*cbulk*z0*(1-f1)
         end if

                G_t_np1=-c1*alpha0*log(1.0d0/alpha0)*(1-f2)

                Hb_t_np1(i,j,k)=Hb_t0*exp(-g0)
     &                       +(F_t_np1)*(1-exp(-g0))
     &                       +(G_t_np1)*(1-exp(-g1))
              end if
             end do
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
     &                      gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                      gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                      gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                      gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                      gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                      gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                      gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                      Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                      Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                      Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                      Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                      phi1_np1,phi1_n,phi1_nm1,
     &                      L,x,y,z,dt,chr,ex,
     &                      phys_bdy,ghost_width,Nx,Ny,Nz,
     &                      Hb_t_0,Hb_x_0,Hb_y_0,
     &                      Hb_z_0,
     &                      gauge,t_n,rho1,rho2,rho3,rho4,xi1,xi2,
     &                      c1,c2,c3,cbulk)
        implicit none
        integer Nx,Ny,Nz,gauge
        integer phys_bdy(6),ghost_width(6)
        real*8 c1,c2,c3
        real*8 res(Nx,Ny,Nz),t_n,t_np1
        real*8 Hb_t_0(Nx,Ny,Nz),Hb_x_0(Nx,Ny,Nz),Hb_y_0(Nx,Ny,Nz)
        real*8 Hb_z_0(Nx,Ny,Nz)
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_x_np1(Nx,Ny,Nz),Hb_y_np1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz)
        real*8 Hb_t_n(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz)
        real*8 Hb_z_n(Nx,Ny,Nz)
        real*8 Hb_t_nm1(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_nm1(Nx,Ny,Nz)
        real*8 chr(Nx,Ny,Nz),ex,L
        real*8 x(Nx),y(Ny),z(Nz),dt,rho1,rho2,rho3,rho4,xi1,xi2,cbulk
        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tx_np1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_xx_np1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_yy_np1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz)
        real*8 gb_tt_n(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz)
        real*8 gb_ty_n(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz)
        real*8 gb_tz_n(Nx,Ny,Nz)
        real*8 gb_xz_n(Nx,Ny,Nz)
        real*8 gb_yz_n(Nx,Ny,Nz)
        real*8 gb_yy_n(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz)
        real*8 gb_tt_nm1(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_nm1(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_xy_nm1(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_zz_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)

        real*8 Hb_t0,Hb_x0,Hb_y0
        real*8 Hb_z0

        real*8 F_t_np1,F_x_np1,F_y_np1
        real*8 F_z_np1
        real*8 G_t_np1,G_x_np1,G_y_np1
        real*8 G_z_np1

        real*8 g0uu11,gbtt,gbtx,gbty,gbtz
        real*8 gbxx,gbxy,gbxz,gbyy,gbyz,gbzz
        real*8 alphasq0,alpha0,betax0,betay0,betaz0

        real*8 x0,y0,z0,rho0

        integer i,j,k,i1,j1
        real*8 dx,dy,dz
        real*8 f0,f1,f2,g0,g1,g2,trans

        real*8 PI
        parameter (PI=3.141592653589793d0)

        ! initialize fixed-size variables
        data i,j,k,i1,j1/0,0,0,0,0/

        data Hb_t0,Hb_x0,Hb_y0/0.0,0.0,0.0/
        data Hb_z0/0.0/
        data F_t_np1,F_x_np1,F_y_np1/0.0,0.0,0.0/
        data F_z_np1/0.0/
        data G_t_np1,G_x_np1,G_y_np1/0.0,0.0,0.0/
        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/
        data f0,f1,f2,g0,g1,g2/0.0,0.0,0.0,0.0,0.0,0.0/
        data dx,dy,dz/0.0,0.0,0.0/

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        t_np1=t_n+dt

        ! initialize output variables
        do i=1,Nx
          do j=1,Ny
           do k=1,Nz
            res(i,j,k)=0
            if (chr(i,j,k).eq.ex) then
              Hb_x_np1(i,j,k)=0
              Hb_y_np1(i,j,k)=0
              Hb_z_np1(i,j,k)=0
            end if
           end do
          end do
        end do

        !--------------------------------------------------------------
        ! gauge 0
        !--------------------------------------------------------------

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 1
        ! (had tried metricindep=22, zeroed init nm1 (except 22), and set
        !  F_x_np1=gb_xx_np1(i,j,k)*x0/rho0*(150)/100
        !  but did not work; something is not right with EFE_cartesian.nb)
        !--------------------------------------------------------------

        if (gauge.eq.1) then
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              if (chr(i,j,k).ne.ex) then
                x0=x(i)
                y0=y(j)
                rho0=sqrt(x0**2+y0**2)

                F_x_np1=
     &                 gb_xx_np1(i,j,k)*x0**2*rho0*(200)/100
                F_y_np1=
     &                 gb_yy_np1(i,j,k)*y0/rho0*(-100)/100

                Hb_x_np1(i,j,k)=F_x_np1
                Hb_y_np1(i,j,k)=F_y_np1
              end if
             end do
            end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! gauge 2 (works for stabilizing gb_tt,gb_tx,gb_ty,psi)
        !         (rho_cd=-1 stabilizes left corner of gb_tt,gb_tx,gb_ty,psi)
        !         (rho_cd=0  stabilizes gb_yy)
        !         (rho_cd=0  stabilizes gb_xx asymmetry)
        !--------------------------------------------------------------

        if (gauge.eq.2) then
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              if (chr(i,j,k).ne.ex) then
                x0=x(i)
                y0=y(j)
                rho0=sqrt(x0**2+y0**2)

                F_x_np1=
     &                 +gb_xx_np1(i,j,k)*x0**2*(1-y0**2)*2  !stabilizes gb_xx
                F_y_np1=
     &                 +gb_yy_np1(i,j,k)*2*y0/rho0
     &                 +gb_xy_np1(i,j,k)*4*x0/rho0

                Hb_x_np1(i,j,k)=F_x_np1
                Hb_y_np1(i,j,k)=F_y_np1
              end if
             end do
            end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! gauge 3
        !--------------------------------------------------------------

        if (gauge.eq.3) then
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              if (chr(i,j,k).ne.ex) then
                x0=x(i)
                y0=y(j)
                z0=z(k)
                rho0=sqrt(x0**2+y0**2+z0**2)
                Hb_x0=Hb_x_0(i,j,k)
                Hb_y0=Hb_y_0(i,j,k)
                Hb_z0=Hb_z_0(i,j,k)

                f0=trans(rho0,rho1,rho2)
                f1=trans(rho0,rho3,rho4)
             !time transition with transition near the AdS boundary   
             !use when using initial data that DOES NOT preserve the Hamiltonian and momentum constraints
!                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4
             !time transition with NO transition near the AdS boundary
             ! use when using initial data that preserves the Hamiltonian and momentum constraints
                g0=(t_np1/(xi1*(1-f0)))**4

               if (rho0.ne.0.0d0) then

                F_x_np1=gb_xx_np1(i,j,k)*1.5d0*x0/rho0*f1
     &                 +gb_xy_np1(i,j,k)*1.5d0*y0/rho0*f1
     &                 +gb_xz_np1(i,j,k)*1.5d0*z0/rho0*f1
     &                 +gb_xx_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_xy_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xz_np1(i,j,k)*cbulk*z0*(1-f1)
                F_y_np1=gb_yy_np1(i,j,k)*1.5d0*y0/rho0*f1
     &                 +gb_xy_np1(i,j,k)*1.5d0*x0/rho0*f1
     &                 +gb_yz_np1(i,j,k)*1.5d0*z0/rho0*f1
     &                 +gb_yy_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xy_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_yz_np1(i,j,k)*cbulk*z0*(1-f1)
                F_z_np1=gb_yz_np1(i,j,k)*1.5d0*y0/rho0*f1
     &                 +gb_xz_np1(i,j,k)*1.5d0*x0/rho0*f1
     &                 +gb_zz_np1(i,j,k)*1.5d0*z0/rho0*f1
     &                 +gb_yz_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xz_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_zz_np1(i,j,k)*cbulk*z0*(1-f1)
               else 
                F_x_np1=gb_xx_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_xy_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xz_np1(i,j,k)*cbulk*z0*(1-f1)
                F_y_np1=gb_yy_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xy_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_yz_np1(i,j,k)*cbulk*z0*(1-f1)
                F_z_np1=gb_yz_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xz_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_zz_np1(i,j,k)*cbulk*z0*(1-f1)
               end if

                Hb_x_np1(i,j,k)=F_x_np1+(Hb_x0-F_x_np1)*exp(-g0)
                Hb_y_np1(i,j,k)=F_y_np1+(Hb_y0-F_y_np1)*exp(-g0)
                Hb_z_np1(i,j,k)=F_z_np1+(Hb_z0-F_z_np1)*exp(-g0)
              end if
            end do
           end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! gauge 4
        !--------------------------------------------------------------

        if (gauge.eq.4) then
          do i=2,Nx-1
            do j=2,Ny-1
             do k=2,Nz-1
              if (chr(i,j,k).ne.ex) then
                x0=x(i)
                y0=y(j)
                z0=z(k)
                rho0=sqrt(x0**2+y0**2+z0**2)
                Hb_x0=Hb_x_0(i,j,k)
                Hb_y0=Hb_y_0(i,j,k)
                Hb_z0=Hb_z_0(i,j,k)

                f0=trans(rho0,rho1,rho2)
                f1=trans(rho0,rho3,rho4)
                f2=trans(rho0,0.0d0,c2)
                g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4
                g1=(t_np1/xi1)**4
                g2=(t_np1/c3)**4

                gbtt =gb_tt_np1(i,j,k)
                gbtx =gb_tx_np1(i,j,k)
                gbty =gb_ty_np1(i,j,k)
                gbtz =gb_tz_np1(i,j,k)
                gbxx =gb_xx_np1(i,j,k)
                gbxy =gb_xy_np1(i,j,k)
                gbxz =gb_xz_np1(i,j,k)
                gbyy =gb_yy_np1(i,j,k)
                gbyz =gb_yz_np1(i,j,k)
                gbzz=gb_zz_np1(i,j,k)
                g0uu11=(-gbxz**2*(4+gbyy)-4*gbyz**2+4*gbyy*gbzz
     &               +16*(4+gbyy+gbzz)+2*(-16*gbyy+gbxz**2*(8+3*gbyy)
     &               +8*gbyz**2
     &               -8*(2+gbyy)*gbzz)*rho0**2+(-3*gbxz**2*(8+5*gbyy)
     &               +8*(-3*gbyz**2+2*gbzz+gbyy*(2+3*gbzz)))*rho0**4
     &               +4*(gbxz**2*(4+5*gbyy)
     &                +4*gbyz**2-4*gbyy*gbzz)*rho0**6
     &               -(gbxz**2*(4+15*gbyy)+4*gbyz**2
     &                -4*gbyy*gbzz)*rho0**8
     &               +6*gbxz**2*gbyy*rho0**10-gbxz**2*gbyy*rho0**12
     &               +2*gbxy*gbxz*gbyz*(-1+rho0**2)**6
     &               -gbxy**2*(-1+rho0**2)**4*(4+gbzz*(-1+rho0**2)**2)
     &               +gbxx*(-1+rho0**2)**2*(-gbyz**2*(-1+rho0**2)**4
     &               +4*(4+gbzz*(-1+rho0**2)**2)
     &               +gbyy*(-1+rho0**2)**2*(4+gbzz*(-1+rho0**2)**2)))
     &               /((-1+rho0**2)**6
     &                *(gbyz*(-gbtx*(gbtz*gbxy+gbty*gbxz)
     &                +gbtx**2*gbyz+(gbty*gbtz*(4+gbxx*(-1+rho0**2)**2))
     &               /(-1+rho0**2)**2
     &               +((gbtt*(-1+rho0**2)**2
     &               -(1+rho0**2)**2)*(gbxy*gbxz*(-1+rho0**2)**2
     &               -gbyz*(4+gbxx*(-1+rho0**2)**2)))/(-1+rho0**2)**4)
     &               +gbxz*(-gbty*gbtz*gbxy+gbty**2*gbxz-gbtx*gbty*gbyz
     &               +(gbtx*gbtz*(4+gbyy*(-1+rho0**2)**2))
     &               /(-1+rho0**2)**2
     &               +((gbtt*(-1+rho0**2)**2
     &                -(1+rho0**2)**2)*(gbxy*gbyz*(-1+rho0**2)**2
     &               -gbxz*(4+gbyy*(-1+rho0**2)**2)))/(-1+rho0**2)**4)
     &               +(1/((-1+rho0**2)**4))*gbtz*((-1+rho0**2)**2
     &               *(-gbty*gbxy*gbxz*(-1+rho0**2)**2
     &               -gbtx*gbxy*gbyz*(-1+rho0**2)**2
     &                +gbty*gbyz*(4+gbxx*(-1+rho0**2)**2)
     &               +gbtx*gbxz*(4+gbyy*(-1+rho0**2)**2))
     &               -gbtz*(-gbxy**2*(-1+rho0**2)**4
     &               +4*(4+gbyy*(-1+rho0**2)**2)
     &               +gbxx*(-1+rho0**2)**2*(4+gbyy*(-1+rho0**2)**2)))
     &               -(1/((-1+rho0**2)**8))*(4+gbzz*(-1+rho0**2)**2)
     &                *(gbty**2*(4+gbxx)
     &                -gbxy**2+gbxx*gbyy
     &               -2*(gbty**2*(8+3*gbxx)-gbxy**2+gbxx*gbyy)*rho0**2
     &               +(16+3*gbty**2*(8+5*gbxx)+gbxy**2-8*gbyy
     &                -gbxx*(8+gbyy))*rho0**4
     &               -4*(gbty**2*(4+5*gbxx)+gbxy**2-gbxx*gbyy)*rho0**6
     &               +(gbty**2*(4+15*gbxx)+gbxy**2
     &                -gbxx*(-4+gbyy)+4*gbyy)*rho0**8
     &               +2*(gbxy**2-gbxx*(3*gbty**2+gbyy))*rho0**10
     &               +(-gbxy**2+gbxx*(gbty**2+gbyy))*rho0**12
     &               -2*gbtx*gbty*gbxy*(-1+rho0**2)**6
     &               +4*(4+gbxx+gbyy+8*rho0**2)
     &               +gbtx**2*(-1+rho0**2)**4*(4+gbyy*(-1+rho0**2)**2)
     &               -gbtt*(-1+rho0**2)**2*(-gbxy**2*(-1+rho0**2)**4
     &               +4*(4+gbyy*(-1+rho0**2)**2)
     &               +gbxx*(-1+rho0**2)**2*(4+gbyy*(-1+rho0**2)**2)))))
                alphasq0=-1/g0uu11
                if (alphasq0.ge.0) then
                  alpha0=sqrt(alphasq0)
                else
                  alpha0=1e-8
                end if
                betax0=gb_tx_np1(i,j,k)
                betay0=gb_ty_np1(i,j,k)
                betaz0=gb_tz_np1(i,j,k)

               if (rho0.ne.0.0d0) then
                F_x_np1=gb_xx_np1(i,j,k)*1.5d0*x0/rho0*f1
     &                 +gb_xy_np1(i,j,k)*1.5d0*y0/rho0*f1
     &                 +gb_xz_np1(i,j,k)*1.5d0*z0/rho0*f1
     &                 +gb_xx_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_xy_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xz_np1(i,j,k)*cbulk*z0*(1-f1)
                F_y_np1=gb_yy_np1(i,j,k)*1.5d0*y0/rho0*f1
     &                 +gb_xy_np1(i,j,k)*1.5d0*x0/rho0*f1
     &                 +gb_yz_np1(i,j,k)*1.5d0*z0/rho0*f1
     &                 +gb_yy_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xy_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_yz_np1(i,j,k)*cbulk*z0*(1-f1)
                F_z_np1=gb_xz_np1(i,j,k)*1.5d0*x0/rho0*f1
     &                 +gb_yz_np1(i,j,k)*1.5d0*y0/rho0*f1
     &                 +gb_zz_np1(i,j,k)*1.5d0*z0/rho0*f1
     &                 +gb_xz_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_yz_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_zz_np1(i,j,k)*cbulk*z0*(1-f1)
               else
                F_x_np1=gb_xx_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_xy_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xz_np1(i,j,k)*cbulk*z0*(1-f1)
                F_y_np1=gb_yy_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xy_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_yz_np1(i,j,k)*cbulk*z0*(1-f1)
                F_z_np1=gb_yz_np1(i,j,k)*cbulk*y0*(1-f1)
     &                 +gb_xz_np1(i,j,k)*cbulk*x0*(1-f1)
     &                 +gb_zz_np1(i,j,k)*cbulk*z0*(1-f1)
               end if

                G_x_np1=-c1*betax0/alpha0*(1-f2)
                G_y_np1=-c1*betay0/alpha0*(1-f2)
                G_z_np1=-c1*betaz0/alpha0*(1-f2)

                Hb_x_np1(i,j,k)=Hb_x0*exp(-g0)
     &                       +(F_x_np1)*(1-exp(-g0))
     &                       +(G_x_np1)*(1-exp(-g1))
                Hb_y_np1(i,j,k)=Hb_y0*exp(-g0)
     &                       +(F_y_np1)*(1-exp(-g0))
     &                       +(G_y_np1)*(1-exp(-g1))
                Hb_z_np1(i,j,k)=Hb_z0*exp(-g0)
     &                       +(F_z_np1)*(1-exp(-g0))
     &                       +(G_z_np1)*(1-exp(-g1))
              end if
             end do
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
