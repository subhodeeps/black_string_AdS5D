c----------------------------------------------------------------------
c calculates uniform AdS5D-black string metric with L=1 and its first and second derivatives in Cartesian Kerr-Schild coordinates (horizon-penetrating) 
c on the hypersurface S at omega_1=z=0, y>=0. 
c The other coordinates are time t, and the radial coordinate x.
c
c r0 is the radius parameter, i.e. r0=2*M0, where M0 is the BH mass parameter
c (r0 has no physical meaning, it is NOT the horizon radius)
c----------------------------------------------------------------------
        subroutine unibsads5_derivs_kerrschildcoords(
     &                  gunibsads5_ll,gunibsads5_uu,gunibsads5_ll_x,
     &                  gunibsads5_uu_x,gunibsads5_ll_xx,
     &                  Hunibsads5_l,
     &                  gammaunibsads5_ull,
     &                  nsym,x,y,dt,chr,L,ex,Nx,Ny,Nz,i,j,
     &                  ief_bh_r0,
     &                  calc_adv_quant)

        implicit none

        integer Nx,Ny,Nz
        integer i,j,k
        real*8  ief_bh_r0,M0
        logical calc_adv_quant

        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),dt,L

        integer a,b,c,d,e,f,g,h
        real*8 dx,dy,dz
        real*8 x0,y0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer nsym

        !--------------------------------------------------------------
        ! variables for tensor manipulations
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------

        real*8 gunibsads5_ll(4,4),gunibsads5_uu(4,4)
        real*8 gunibsads5_ll_x(4,4,4),gunibsads5_uu_x(4,4,4)
        real*8 gunibsads5_ll_xx(4,4,4,4)
        real*8 Hunibsads5_l(4)
        real*8 gunibsads5_ll(4,4),gunibsads5_uu(4,4)
        real*8 gunibsads5_ll_x(4,4,4),gunibsads5_uu_x(4,4,4)
        real*8 gunibsads5_ll_xx(4,4,4,4)
        real*8 dxqssph_dxcar(4,4)
        real*8 d2xqssph_dxcardxcar(4,4,4)
        real*8 d3xqssph_dxcardxcardxcar(4,4,4,4)
        real*8 gammaunibsads5_ull(4,4,4)
        real*8 boxunibsads5x_u(4)
        real*8 gammaunibsads5_ull_x(4,4,4,4)
        real*8 riemannunibsads5_ulll(4,4,4,4)
        real*8 ricciunibsads5_ll(4,4),ricciunibsads5_lu(4,4)
        real*8 ricciunibsads5
        real*8 einsteinunibsads5_ll(4,4)

        real*8 detg0_unibsads50

!----------------------------------------------------------------------

        
        dx=(x(2)-x(1))
        dy=(y(2)-y(1))

        x0=x(i)
        y0=y(j)

      ! Schwarzschild mass parameter
        M0=ief_bh_r0/2

!hard-coding metric components in asymptotically Cartesian Kerr-Schild coordinates (horizon penetrating)

        ! give values to the unibsads5 metric
        gunibsads5_ll(1,1)=-0.5*((1 + x0**2)**2*(ief_bh_r0*(-1 + y0**2)**3 
     &         + 2*y0*(1 + y0**2)**2))/((-1 + x0**2)**2*y0*(-1 + y0**2)**2)
        gunibsads5_ll(1,2)=0
        gunibsads5_ll(1,3)=-((ief_bh_r0*(1 + x0**2)**2*(-1 + y0**2))/((-1 + x0**2)**2*(y0 + y0**3)))
        gunibsads5_ll(2,2)=4/(-1 + x0**2)**2
        gunibsads5_ll(2,3)=0
        gunibsads5_ll(3,3)=(-2*(1 + x0**2)**2*(ief_bh_r0*(-1 + y0**2)**3 
     &         - 2*y0*(1 + y0**2)**2))/((-1 + x0**2)**2*y0*(-1 + y0**4)**2)

        gunibsads5_lzlz=(4*(1 + x0**2)**2)/((-1 + x0**2)**2*(-1 + y0**2)**2)

        do a=1,3
          do b=1,3
            gunibsads5_ll(a,b)=gunibsads5_ll(min(a,b),max(a,b))
          end do
        end do


        !first derivatives
        gunibsads5_ll_x(1,1,1)   =0
        gunibsads5_ll_x(1,1,2)   =(4*x0*(1 + x0**2)*(ief_bh_r0*(-1 + y0**2)**3 
     &         + 2*y0*(1 + y0**2)**2))/((-1 + x0**2)**3*y0*(-1 + y0**2)**2)
        gunibsads5_ll_x(1,1,3)   =-0.5*((1 + x0**2)**2*(1 + y0**2)*(-16*y0**3 
     &         + ief_bh_r0*(-1 + y0**2)**3))/((-1 + x0**2)**2*y0**2*(-1 + y0**2)**3)

        gunibsads5_ll_x(1,2,1)   =0
        gunibsads5_ll_x(1,2,2)   =0
        gunibsads5_ll_x(1,2,3)   =0

        gunibsads5_ll_x(1,3,1)   =0
        gunibsads5_ll_x(1,3,2)   =(8*ief_bh_r0*x0*(1 + x0**2)*(-1 + y0**2))
     &         /((-1 + x0**2)**3*(y0 + y0**3))
        gunibsads5_ll_x(1,3,3)   =(ief_bh_r0*(1 + x0**2)**2*(-1 - 4*y0**2 + y0**4))
     &         /((-1 + x0**2)**2*y0**2*(1 + y0**2)**2)

        gunibsads5_ll_x(2,2,1)   =0
        gunibsads5_ll_x(2,2,2)   =(-16*x0)/(-1 + x0**2)**3
        gunibsads5_ll_x(2,2,3)   =0

        gunibsads5_ll_x(2,3,1)   =0
        gunibsads5_ll_x(2,3,2)   =0
        gunibsads5_ll_x(2,3,3)   =0

        gunibsads5_ll_x(3,3,1)   =0
        gunibsads5_ll_x(3,3,2)   =(16*x0*(1 + x0**2)*(ief_bh_r0*(-1 + y0**2)**3 
     &         - 2*y0*(1 + y0**2)**2))/((-1 + x0**2)**3*y0*(-1 + y0**4)**2)
        gunibsads5_ll_x(3,3,3)   =(2*(1 + x0**2)**2*(-8*(y0 + y0**3)**3 + ief_bh_r0*
     &         (-1 + y0**2)**3*(-1 - 6*y0**2 + 3*y0**4)))/((-1 + x0**2)**2*y0**2*(-1 + y0**4)**3)

        gunibsads5_lzlz_x(1)     =0
        gunibsads5_lzlz_x(2)     =(-32*(x0 + x0**3))/((-1 + x0**2)**3*(-1 + y0**2)**2)
        gunibsads5_lzlz_x(3)     =(-16*(1 + x0**2)**2*y0)/((-1 + x0**2)**2*(-1 + y0**2)**3)

        !second derivatives
        gunibsads5_ll_xx(1,1,1,1)=0
        gunibsads5_ll_xx(1,1,1,2)=0
        gunibsads5_ll_xx(1,1,1,3)=0
        gunibsads5_ll_xx(1,1,2,2)=(-4*(1 + 8*x0**2 + 3*x0**4)*(ief_bh_r0*(-1 + y0**2)**3 
     &         + 2*y0*(1 + y0**2)**2))/((-1 + x0**2)**4*y0*(-1 + y0**2)**2)
        gunibsads5_ll_xx(1,1,2,3)=(4*x0*(1 + x0**2)*(1 + y0**2)*(-16*y0**3 
     &         + ief_bh_r0*(-1 + y0**2)**3))/((-1 + x0**2)**3*y0**2*(-1 + y0**2)**3)
        gunibsads5_ll_xx(1,1,3,3)=((1 + x0**2)**2*(ief_bh_r0*(-1 + y0**2)**4 
     &         - 8*(y0**3 + 8*y0**5 + 3*y0**7)))/((-1 + x0**2)**2*y0**3*(-1 + y0**2)**4)

        gunibsads5_ll_xx(1,2,1,1)=0
        gunibsads5_ll_xx(1,2,1,2)=0
        gunibsads5_ll_xx(1,2,1,3)=0
        gunibsads5_ll_xx(1,2,2,2)=0
        gunibsads5_ll_xx(1,2,2,3)=0
        gunibsads5_ll_xx(1,2,3,3)=0

        gunibsads5_ll_xx(1,3,1,1)=0
        gunibsads5_ll_xx(1,3,1,2)=0
        gunibsads5_ll_xx(1,3,1,3)=0
        gunibsads5_ll_xx(1,3,2,2)=(-8*ief_bh_r0*(1 + 8*x0**2 + 3*x0**4)*(-1 + y0**2))
     &         /((-1 + x0**2)**4*(y0 + y0**3))
        gunibsads5_ll_xx(1,3,2,3)=(-8*ief_bh_r0*x0*(1 + x0**2)*(-1 - 4*y0**2 + y0**4))
     &         /((-1 + x0**2)**3*(y0 + y0**3)**2)
        gunibsads5_ll_xx(1,3,3,3)=(-2*ief_bh_r0*(1 + x0**2)**2*(-1 - 3*y0**2 
     &         - 9*y0**4 + y0**6))/((-1 + x0**2)**2*(y0 + y0**3)**3)

        gunibsads5_ll_xx(2,2,1,1)=0
        gunibsads5_ll_xx(2,2,1,2)=0
        gunibsads5_ll_xx(2,2,1,3)=0
        gunibsads5_ll_xx(2,2,2,2)=(16*(1 + 5*x0**2))/(-1 + x0**2)**4
        gunibsads5_ll_xx(2,2,2,3)=0
        gunibsads5_ll_xx(2,2,3,3)=0

        gunibsads5_ll_xx(2,3,1,1)=0
        gunibsads5_ll_xx(2,3,1,2)=0
        gunibsads5_ll_xx(2,3,1,3)=0
        gunibsads5_ll_xx(2,3,2,2)=0
        gunibsads5_ll_xx(2,3,2,3)=0
        gunibsads5_ll_xx(2,3,3,3)=0

        gunibsads5_ll_xx(3,3,1,1)=0
        gunibsads5_ll_xx(3,3,1,2)=0
        gunibsads5_ll_xx(3,3,1,3)=0
        gunibsads5_ll_xx(3,3,2,2)=(-16*(1 + 8*x0**2 + 3*x0**4)*(ief_bh_r0*(-1 + y0**2)**3 
     &         - 2*y0*(1 + y0**2)**2))/((-1 + x0**2)**4*y0*(-1 + y0**4)**2)
        gunibsads5_ll_xx(3,3,2,3)=(-16*x0*(1 + x0**2)*(-8*(y0 + y0**3)**3 
     &         + ief_bh_r0*(-1 + y0**2)**3*(-1 - 6*y0**2 + 3*y0**4)))
     &          /((-1 + x0**2)**3*y0**2*(-1 + y0**4)**3)
        gunibsads5_ll_xx(3,3,3,3)=(-4*(1 + x0**2)**2*(-4*y0**3*(1 + y0**2)**4*(1 + 5*y0**2) 
     &         + ief_bh_r0*(-1 + y0**2)**4*(-1 - 4*y0**2 - 21*y0**4 + 6*y0**6)))
     &         /((-1 + x0**2)**2*y0**3*(-1 + y0**4)**4)

        gunibsads5_lzlz_xx(1,1)  =0
        gunibsads5_lzlz_xx(1,2)  =0
        gunibsads5_lzlz_xx(1,3)  =0
        gunibsads5_lzlz_xx(2,2)  =(32*(1 + 8*x0**2 + 3*x0**4))/((-1 + x0**2)**4*(-1 + y0**2)**2)
        gunibsads5_lzlz_xx(2,3)  =(128*x0*(1 + x0**2)*y0)/((-1 + x0**2)**3*(-1 + y0**2)**3)
        gunibsads5_lzlz_xx(3,3)  =(16*(1 + x0**2)**2*(1 + 5*y0**2))/((-1 + x0**2)**2*(-1 + y0**2)**4)


        do a=1,3
          do b=1,3
            gunibsads5_lzlz_xx(a,b)=gunibsads5_lzlz_xx(min(a,b),max(a,b))
             do c=1,3
               gunibsads5_ll_x(a,b,c)=gunibsads5_ll_x(min(a,b),max(a,b),c)
               do d=1,3
                gunibsads5_ll_xx(a,b,c,d)=
     &              gunibsads5_ll_xx(min(a,b),max(a,b),min(c,d),max(c,d))
               end do
          end do
        end do

        
        return
        end
!--------------------------------------------------------------------------------