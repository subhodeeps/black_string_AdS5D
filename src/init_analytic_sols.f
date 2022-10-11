c----------------------------------------------------------------------
c Initializes analytic solutions
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c computes the gb function to obtain an exact Schwarzschild-AdS black hole solution in horizon penetrating 
c coordinates as deviation from pure AdS in standard Cartesian coordinates.
c The Cartesian horizon penetrating coordinates are defined as follows: from the (compactified) spherical coordinates 
c (t,rho,theta,phi) we go to horizon penetrating spherical coords (u,rho,theta,phi), defined by 
c dt = du - 1/H dr, with H=U ((1 -rho_h)/(1 -rho))^n, where n=3 to obtain gb that go to 0 linearly in 1-rho
c near the AdS boundary. rho_h is the rho-value of the event horizon. 
c Then we go from (u,r,theta,phi) to (u,x,y,z) with 
c x=rho cos(theta),y=rho sin(theta) cos(phi),z=rho sin(theta) sin(phi)
c 
c Radius parameter r0=2*M0, where M0 is the BHmass (r0 has no physical meaning, it is NOT 
c the horizon radius.).
c We denote the horizon radius by r_h in non-compactified coordinates and rho_h in compactified coordinates.
c----------------------------------------------------------------------
        subroutine init_schwads4d_bh_adhoc_coords(ief_bh_r0,L,
     &                         gb_tt,gb_tx,
     &                         gb_ty,
     &                         gb_tz,
     &                         gb_xx,gb_xy,
     &                         gb_xz,
     &                         gb_yy,
     &                         gb_yz,
     &                         gb_zz,gb_tt_t,gb_tx_t,gb_ty_t,
     &                         gb_tz_t,
     &                         gb_xx_t,gb_xy_t,
     &                         gb_xz_t,
     &                         gb_yy_t,
     &                         gb_yz_t,
     &                         gb_zz_t,Hb_t,Hb_x,Hb_y,
     &                         Hb_z,
     &                         Hb_t_t,Hb_x_t,Hb_y_t,
     &                         Hb_z_t,
     &                         phys_bdy,
     &                         x,y,z,dt,chr,ex,Nx,Ny,Nz,regtype)
        implicit none

        integer Nx,Ny,Nz
        integer regtype
        integer phys_bdy(6)
        real*8 ief_bh_r0
        real*8 dt,ex,L
        real*8 chr(Nx,Ny,Nz)
        real*8 Hb_t(Nx,Ny,Nz),Hb_x(Nx,Ny,Nz)
        real*8 Hb_y(Nx,Ny,Nz)
        real*8 Hb_z(Nx,Ny,Nz)
        real*8 Hb_t_t(Nx,Ny,Nz),Hb_x_t(Nx,Ny,Nz)
        real*8 Hb_y_t(Nx,Ny,Nz)
        real*8 Hb_z_t(Nx,Ny,Nz)
        real*8 gb_tt(Nx,Ny,Nz),gb_tx(Nx,Ny,Nz)
        real*8 gb_ty(Nx,Ny,Nz)
        real*8 gb_tz(Nx,Ny,Nz)
        real*8 gb_xx(Nx,Ny,Nz),gb_xy(Nx,Ny,Nz)
        real*8 gb_xz(Nx,Ny,Nz)
        real*8 gb_yy(Nx,Ny,Nz),gb_zz(Nx,Ny,Nz)
        real*8 gb_yz(Nx,Ny,Nz)
        real*8 gb_tt_t(Nx,Ny,Nz),gb_tx_t(Nx,Ny,Nz)
        real*8 gb_ty_t(Nx,Ny,Nz),gb_zz_t(Nx,Ny,Nz)
        real*8 gb_tz_t(Nx,Ny,Nz)
        real*8 gb_xx_t(Nx,Ny,Nz),gb_xy_t(Nx,Ny,Nz)
        real*8 gb_xz_t(Nx,Ny,Nz)
        real*8 gb_yy_t(Nx,Ny,Nz)
        real*8 gb_yz_t(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz)

        integer n
        parameter (n=3)

        integer i,j,k

        real*8 rho0,f1,f0,cF0,C0,A0,B0,D0
        real*8 x0,y0,z0
        real*8 r_h,rho_h
        real*8 small
        parameter (small=1d-10)

        real*8 tfunction(Nx,Ny,Nz)

        ! initialize fixed-size variables
        data i,j,k/0,0,0/

        data rho0,f1,f0,cF0/0.0,0.0,0.0,0.0/
        data C0,A0,B0,D0/0.0,0.0,0.0,0.0/
        data x0,y0,z0/0.0,0.0,0.0/
        data r_h,rho_h/0.0,0.0/

        !--------------------------------------------------------------

        ! compute horizon global radius r_h and corresponding compactified rho_h
        !(NOTE ... here the x0,y0,rho0 are compactified (CODE)
        ! versions of the coordinates)

        r_h=-L**2
     &   /(3**(1.0d0/3.0d0))
     &   /((9*L**2*(ief_bh_r0/2)
     &    +sqrt(3.0d0)*sqrt(L**6+27*L**4*(ief_bh_r0/2)**2))
     &    **(1.0d0/3.0d0))
     &   +((9*L**2*(ief_bh_r0/2)
     &   +sqrt(3.0d0)*sqrt(L**6+27*L**4*(ief_bh_r0/2)**2))
     &    **(1.0d0/3.0d0))
     &   /(3**(2.0d0/3.0d0))

        rho_h=(-1 + sqrt(1 + r_h**2))/r_h

        ! initialize metric 
        do i=1,Nx
           do j=1,Ny
            do k=1,Nz
              gb_tt_t(i,j,k)=0
              gb_tx_t(i,j,k)=0
              gb_ty_t(i,j,k)=0
              gb_tz_t(i,j,k)=0
              gb_xx_t(i,j,k)=0
              gb_xy_t(i,j,k)=0
              gb_xz_t(i,j,k)=0
              gb_yy_t(i,j,k)=0
              gb_yz_t(i,j,k)=0
              gb_zz_t(i,j,k)=0
              Hb_t_t(i,j,k)=0
              Hb_x_t(i,j,k)=0
              Hb_y_t(i,j,k)=0
              Hb_z_t(i,j,k)=0
              if (chr(i,j,k).eq.ex) then
                 gb_tt(i,j,k)=0
                 gb_tx(i,j,k)=0
                 gb_ty(i,j,k)=0
                 gb_tz(i,j,k)=0
                 gb_xx(i,j,k)=0
                 gb_xy(i,j,k)=0
                 gb_xz(i,j,k)=0
                 gb_yy(i,j,k)=0
                 gb_yz(i,j,k)=0
                 gb_zz(i,j,k)=0
                 Hb_t(i,j,k)=0
                 Hb_x(i,j,k)=0
                 Hb_y(i,j,k)=0
                 Hb_z(i,j,k)=0
              else
                 x0=x(i)
                 y0=y(j)
                 z0=z(k)
                 rho0=sqrt(x0**2+y0**2+z0**2)

                 ! EF-like-near-horizon Schwarzschild-like-near-bdy coordinates

!!!2+1 version!!!!!!!!                 
!                 gb_tt(i,j,k)=(ief_bh_r0/2)*(1/rho0-rho0)
!
!                 gb_tx(i,j,k)=2*x0*(1+rho0**2)
!     &                      *((1-rho_h)/(1-rho0))**(-n)
!     &                      /rho0/(1-rho0**2)**2
!                 gb_ty(i,j,k)=2*y0*(1+rho0**2)
!     &                      *((1-rho_h)/(1-rho0))**(-n)
!     &                      /rho0/(1-rho0**2)**2
!                 gb_tz(i,j,k)=0
!                 gb_xx(i,j,k)=4*x0**2*(1+rho0**2)**2
!     &                      *(-1/(1+(-2+4/L**2)*rho0**2+rho0**4)
!     &                      +L**2*rho0*(1-((1-rho_h)/(1-rho0))**(-2*n))
!     &                      /(4*rho0**3+L**2*(1-rho0**2)**2
!     &                       *(rho0+(r0/2)*(-1+rho0**2))))
!     &                      /rho0**2/(1-rho0**2)**2
!                 gb_xy(i,j,k)=4*L**2*x0*y0*(1+rho0**2)**2
!     &                      *(-4*rho0**3+L**2*(1-rho0**2)**2
!     &                      *(-rho0-(r0/2)*(-1+rho0**2)
!     &                      *((1-rho_h)/(1-rho0))**(2*n)))
!     &                      *((1-rho_h)/(1-rho0))**(-2*n)
!     &                      /rho0**2/(1-rho0**2)**2
!     &                      /(4*rho0**2+L**2*(1-rho0**2)**2)
!     &                      /(4*rho0**3+L**2*(1-rho0**2)**2
!     &                       *(rho0+(r0/2)*(-1+rho0**2)))
!                 gb_xz(i,j,k)=0
!                 gb_yy(i,j,k)=4*y0**2*(1+rho0**2)**2
!     &                      *(-1/(1+(-2+4/L**2)*rho0**2+rho0**4)
!     &                      +L**2*rho0*(1-((1-rho_h)/(1-rho0))**(-2*n))
!     &                      /(4*rho0**3+L**2*(1-rho0**2)**2
!     &                       *(rho0+(r0/2)*(-1+rho0**2))))
!     &                      /rho0**2/(1-rho0**2)**2
!                 gb_yz(i,j,k)=0
!                 psi(i,j,k)=0
!!!!!!!!!!!!!!!!!!!!!!!

!!!CHECKED WITH Mathematica!!

                 gb_tt(i,j,k)=(ief_bh_r0/2)*(1/rho0-rho0)
                 gb_tx(i,j,k)=2*x0*(1+rho0**2)
     &                      *((-1+rho_h)/(-1+rho0))**(-n)
     &                      /rho0/(-1+rho0**2)**2
                 gb_ty(i,j,k)=2*y0*(1+rho0**2)
     &                      *((-1+rho_h)/(-1+rho0))**(-n)
     &                      /rho0/(-1+rho0**2)**2
                 gb_tz(i,j,k)=2*z0*(1+rho0**2)
     &                      *((-1+rho_h)/(-1+rho0))**(-n)
     &                      /rho0/(-1+rho0**2)**2
                 gb_xx(i,j,k)=(-((8*(-1+L**2)*(x0**2-y0**2-z0**2)
     &                   +8*rho0**2+4*L**2*(1+rho0**4))
     &                   /(4*rho0**2+L**2*(-1+rho0**2)**2))
     &                   +(4*(y0**2+z0**2
     &                    +(L**2*x0**2*rho0*(1+rho0**2)**2
     &                    *(-1+((-1+rho_h)/(-1+rho0))**(2*n))
     &                    *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                   /(4*rho0**3+L**2*(-1+rho0**2)**2
     &                   *(rho0+(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2)))))
     &                   /rho0**2)
     &                   /(-1+rho0**2)**2
                 gb_xy(i,j,k)=(4*L**2*x0*y0*(1+rho0**2)**2
     &                   *(-4*rho0**3+L**2*(-1+rho0**2)**2
     &                   *(-rho0-(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2)
     &                   *((-1+rho_h)/(-1+rho0))**(2*n)))
     &                   *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                   /(rho0**2*(-1+rho0**2)**2*(4*rho0**2
     &                   +L**2*(-1+rho0**2)**2)*(4*rho0**3
     &                   +L**2*(-1+rho0**2)**2
     &                   *(rho0+(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2))))
                 gb_xz(i,j,k)=(4*L**2*x0*z0*(1+rho0**2)**2
     &                   *(-4*rho0**3+L**2*(-1+rho0**2)**2
     &                   *(-rho0-(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2)
     &                   *((-1+rho_h)/(-1+rho0))**(2*n)))
     &                   *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                   /(rho0**2*(-1+rho0**2)**2*(4*rho0**2
     &                   +L**2*(-1+rho0**2)**2)*(4*rho0**3
     &                   +L**2*(-1+rho0**2)**2
     &                   *(rho0+(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2))))
           if (y0.ne.0.0d0) then
                 gb_yy(i,j,k)=(1/((-1+rho0**2)**2))
     &                   *((8*y0**2*(-x0**2+y0**2+z0**2)
     &                   -8*(y0**2+2*z0**2)*rho0**2
     &                   -4*L**2*(2*y0**4+z0**2*(-1+rho0**2)**2
     &                   +y0**2*(1-2*x0**2+2*z0**2+rho0**4)))
     &                   /((y0**2+z0**2)*(4*rho0**2
     &                    +L**2*(-1+rho0**2)**2))
     &                   +4*(z0**2/(y0**2+z0**2)
     &                   +(x0**2*y0**2)/((y0**2+z0**2)*rho0**2)
     &                   +(L**2*y0**2*(1+rho0**2)**2
     &                    *(-1+((-1+rho_h)/(-1+rho0))**(2*n))
     &                    *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                   /(4*rho0**4+L**2*rho0*(-1+rho0**2)**2
     &                 *(rho0+(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2)))))
          else
                 gb_yy(i,j,k)=0.0d0
          end if
                 gb_yz(i,j,k)=(4*L**2*y0*z0*(1+rho0**2)**2
     &                   *(-4*rho0**3+L**2*(-1+rho0**2)**2
     &                   *(-rho0-(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2)
     &                   *((-1+rho_h)/(-1+rho0))**(2*n)))
     &                   *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                   /(rho0**2*(-1+rho0**2)**2*(4*rho0**2
     &                   +L**2*(-1+rho0**2)**2)*(4*rho0**3
     &                   +L**2*(-1+rho0**2)**2
     &                   *(rho0+(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2))))
           if (z0.ne.0.0d0) then
                 gb_zz(i,j,k)=(1/((-1+rho0**2)**2))
     &                   *(
     &                    (
     &                    8*z0**2*(-x0**2+y0**2+z0**2)
     &                   -8*(2*y0**2+z0**2)*rho0**2
     &                   -4*L**2*(z0**2*(1-2*x0**2
     &                    +2*z0**2+rho0**4)
     &                    +y0**2*(2*z0**2+(-1+rho0**2)**2))
     &                    )
     &                   /((y0**2+z0**2)*(4*rho0**2
     &                    +L**2*(-1+rho0**2)**2))
     &                   +4*(y0**2/(y0**2+z0**2)
     &                   +(x0**2*z0**2)/((y0**2+z0**2)*rho0**2)
     &                   +(L**2*z0**2*(1+rho0**2)**2
     &                    *(-1+((-1+rho_h)/(-1+rho0))**(2*n))
     &                    *((-1+rho_h)/(-1+rho0))**(-2*n))
     &                   /(4*rho0**4+L**2*rho0*(-1+rho0**2)**2
     &                    *(rho0+(1.0d0/2.0d0)*ief_bh_r0*(-1+rho0**2))))
     &                   )
           else
                 gb_zz(i,j,k)=0.0d0
           end if


              end if
            end do
           end do
        end do

        ! y=0 axis regularization
!        call axi_reg_g(gb_tt,gb_tx,gb_ty,
!     &                 gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
!     &                 L,x,y,z,Nx,Ny,Nz,regtype)

!        call axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)

        return
        end


c----------------------------------------------------------------------
c sets the horizon radius in spherical (non-rotating at the boundary) 
c compactified coordinates for an exact Kerr-AdS black hole solution
c with radius parameter rbh=2*M0, where M0 is the BHmass (r0 has no physical meaning, it is NOT the horizon radius),
c and rotation parameter a_rot.
c----------------------------------------------------------------------
        subroutine set_kerrads4d_ahr(ief_bh_r0,a_rot,L,
     &                         AH_R,AH_xc,min_AH_R,max_AH_R,
     &                         AH_semiax,
     &                         AH_Nchi,AH_Nphi)

        implicit none

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer AH_Nchi,AH_Nphi
        real*8 ief_bh_r0,a_rot,M0,M0_min
        real*8 rblhor
        real*8 AH_R(AH_Nchi,AH_Nphi),AH_xc(3),AH_semiax(3)
        real*8 min_AH_R,max_AH_R
        real*8 AH_chi,AH_phi
        real*8 yipshor,rhohor
        real*8 L
        real*8 dahchi,dahphi
        real*8 min_AH_x,max_AH_x
        real*8 min_AH_y,max_AH_y
        real*8 min_AH_z,max_AH_z
        real*8 AH_x,AH_y,AH_z


        integer i,j,k

        ! initialize fixed-size variables
        data i,j,k/0,0,0/

        !--------------------------------------------------------------



      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "et_kerrads4d_ahr: ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if ((abs(M0).gt.10.0d0**(-10))
     &     .and.(M0.le.M0_min)) then
          write (*,*) "set_kerrads4d_ahr: ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if

      !event horizon radius in Boyer-Lindquist coordinates rotating at the boundary (non-spherical coordinates)
        rblhor=(Sqrt(-2*a_rot**2 - 2*L**2 + 
     -      (a_rot**4 + 14*a_rot**2*L**2 + L**4)/
     - (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 54*L**4*M0**2 +
     -          Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -   4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 + 
     -                54*L**4*M0**2)**2)/2.)**0.3333333333333333 + 
     -   (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 
     -  + 54*L**4*M0**2 + 
     -         Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -            4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                54*L**4*M0**2)**2)/2.)**0.3333333333333333) + 
     -    Sqrt(-4*a_rot**2 - 4*L**2 - 
     -      (a_rot**4 + 14*a_rot**2*L**2 + L**4)/
     -   (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 
     -  + 54*L**4*M0**2 +
     -          Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -             4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                 54*L**4*M0**2)**2)/2.)**0.3333333333333333 - 
     -   (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 
     -  + 54*L**4*M0**2 + 
     -         Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -            4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                54*L**4*M0**2)**2)/2.)**0.3333333333333333 + 
     -      (12*Sqrt(3.)*L**2*M0)/
     -       Sqrt(-2*a_rot**2 - 2*L**2 + 
     -         (a_rot**4 + 14*a_rot**2*L**2 + L**4)/
     -          (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 + L**6 + 
     -             54*L**4*M0**2 + 
     -             Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -                4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                    54*L**4*M0**2)**2)/2.)**0.3333333333333333 + 
     -    (a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  +L**6+54*L**4*M0**2+
     -            Sqrt(-4*(a_rot**4 + 14*a_rot**2*L**2 + L**4)**3 + 
     -               4*(a_rot**6 - 33*a_rot**4*L**2 - 33*a_rot**2*L**4 
     -  + L**6 + 
     -                   54*L**4*M0**2)**2)/2.)**0.3333333333333333)))/
     -  (2.*Sqrt(3.))

        dahchi=PI/(AH_Nchi-1)
        dahphi=2*PI/(AH_Nphi-1)

        min_AH_R=1
        max_AH_R=0

        min_AH_x=1
        max_AH_x=0
        min_AH_y=1
        max_AH_y=0
        min_AH_z=1
        max_AH_z=0

          do i=1,AH_Nchi
           do j=1,AH_Nphi

            !these are the angles obtained from the Cartesian coordinates of the code
            AH_chi=(i-1)*dahchi
            AH_phi=(j-1)*dahphi

            !position of event horizon in uncompactified spherical coordinates (non-rotating at the boundary)
            yipshor=(Sqrt(2.)*L*
     -       Sqrt(rblhor**2*(a_rot**2 + rblhor**2)))/
     -       Sqrt(2*L**2*rblhor**2 + a_rot**2*(L**2 - rblhor**2) 
     -       + a_rot**2*(L**2 + rblhor**2)*Cos(2*AH_chi))

            !position of event horizon in compactified spherical coordinates: rho=sqrt(x**2+y**2+z**2) where x,y,z are the code coordinates
            rhohor=(-1 + Sqrt(1 + yipshor**2))/
     -            yipshor

            AH_R(i,j)=rhohor

            if (AH_R(i,j).gt.max_AH_R) max_AH_R=AH_R(i,j)
            if (AH_R(i,j).lt.min_AH_R) min_AH_R=AH_R(i,j)

            AH_x=AH_R(i,j)*cos(AH_chi)+AH_xc(1)
            AH_y=AH_R(i,j)*sin(AH_chi)*cos(AH_phi)+AH_xc(2)
            AH_z=AH_R(i,j)*sin(AH_chi)*sin(AH_phi)+AH_xc(3)

            min_AH_x=min(AH_x,min_AH_x)
            min_AH_y=min(AH_y,min_AH_y)
            min_AH_z=min(AH_z,min_AH_z)

            max_AH_x=max(AH_x,max_AH_x)
            max_AH_y=max(AH_y,max_AH_y)
            max_AH_z=max(AH_z,max_AH_z)

            !sets semi-axes of ellipsoid approximation of AH:
            !if the AH can be approximated by an ellipsoid with axes along the x,y,z axes, the following gives the semi-axes of this ellipsoid
            AH_semiax(1)=(max_AH_x-min_AH_x)/2
            AH_semiax(2)=(max_AH_y-min_AH_y)/2
            AH_semiax(3)=(max_AH_z-min_AH_z)/2

           end do
          end do

        return
        end



c----------------------------------------------------------------------
c initializes the metric to an exact Kerr-AdS black hole solution
c with radius parameter rbh=2*M0, where M0 is the BHmass (r0 has no physical meaning, 
c it is NOT the horizon radius),
c and rotation parameter a_rot.
c----------------------------------------------------------------------
        subroutine init_kerrads4d_bh(ief_bh_r0,a_rot,L,
     &                         gb_tt,gb_tx,gb_ty,gb_tz,
     &                         gb_xx,gb_xy,gb_xz,
     &                         gb_yy,gb_yz,gb_zz,
     &                         gb_tt_t,gb_tx_t,gb_ty_t,gb_tz_t,
     &                         gb_xx_t,gb_xy_t,gb_xz_t,
     &                         gb_yy_t,gb_yz_t,
     &                         gb_zz_t,
     &                         Hb_t,Hb_x,Hb_y,Hb_z,
     &                         Hb_t_t,Hb_x_t,Hb_y_t,Hb_z_t,
     &                         phys_bdy,
     &                         x,y,z,dt,chr,ex,Nx,Ny,Nz,regtype,
     &                         kerrads_background)
        implicit none

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer Nx,Ny,Nz
        integer regtype
        integer kerrads_background
        integer phys_bdy(6)
        real*8 ief_bh_r0,a_rot,M0,M0_min
        real*8 rblhor
        real*8 dt,ex,L
        real*8 chr(Nx,Ny,Nz)
        real*8 Hb_t(Nx,Ny,Nz),Hb_x(Nx,Ny,Nz)
        real*8 Hb_y(Nx,Ny,Nz)
        real*8 Hb_z(Nx,Ny,Nz)
        real*8 Hb_t_t(Nx,Ny,Nz),Hb_x_t(Nx,Ny,Nz)
        real*8 Hb_y_t(Nx,Ny,Nz)
        real*8 Hb_z_t(Nx,Ny,Nz)
        real*8 gb_tt(Nx,Ny,Nz),gb_tx(Nx,Ny,Nz)
        real*8 gb_ty(Nx,Ny,Nz)
        real*8 gb_tz(Nx,Ny,Nz)
        real*8 gb_xx(Nx,Ny,Nz),gb_xy(Nx,Ny,Nz)
        real*8 gb_xz(Nx,Ny,Nz)
        real*8 gb_yy(Nx,Ny,Nz),gb_zz(Nx,Ny,Nz)
        real*8 gb_yz(Nx,Ny,Nz)
        real*8 gb_tt_t(Nx,Ny,Nz),gb_tx_t(Nx,Ny,Nz)
        real*8 gb_ty_t(Nx,Ny,Nz),gb_zz_t(Nx,Ny,Nz)
        real*8 gb_tz_t(Nx,Ny,Nz)
        real*8 gb_xx_t(Nx,Ny,Nz),gb_xy_t(Nx,Ny,Nz)
        real*8 gb_xz_t(Nx,Ny,Nz)
        real*8 gb_yy_t(Nx,Ny,Nz)
        real*8 gb_yz_t(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz)

        integer a,b,c,d,n
        parameter (n=3)

        integer i,j,k

        real*8 f1,f0,cF0,C0,A0,B0,D0
        real*8 x0,y0,z0,rho0,theta0,phi0
        real*8 U,Sigma,Deltatheta,Xi
        real*8 r_h,rho_h
        real*8 small
        parameter (small=1d-10)

        real*8 h_tt_kerrads_sph0
        real*8 h_trho_kerrads_sph0
        real*8 h_ttheta_kerrads_sph0
        real*8 h_tphi_kerrads_sph0
        real*8 h_rhorho_kerrads_sph0
        real*8 h_rhotheta_kerrads_sph0
        real*8 h_rhophi_kerrads_sph0
        real*8 h_thetatheta_kerrads_sph0
        real*8 h_thetaphi_kerrads_sph0
        real*8 h_phiphi_kerrads_sph0

        real*8 h_kerrads_ll_sph(4,4)
        real*8 h_kerrads_ll(4,4)

        real*8 dxsph_dxcar(4,4)

        logical is_nan

        real*8 tfunction(Nx,Ny,Nz)

        ! initialize fixed-size variables
        data i,j,k/0,0,0/

        data rho0,f1,f0,cF0/0.0,0.0,0.0,0.0/
        data C0,A0,B0,D0/0.0,0.0,0.0,0.0/
        data x0,y0,z0/0.0,0.0,0.0/
        data r_h,rho_h/0.0,0.0/

        !--------------------------------------------------------------

        !(NOTE ... here the x0,y0,z0,rho0 are compactified (CODE)
        ! versions of the coordinates)

      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &    + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &       + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "init_kerrads4d_bh: ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if ((abs(M0).gt.10.0d0**(-10))
     &     .and.(M0.le.M0_min)) then
          write (*,*) "init_kerrads4d_bh: ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger 
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if

        ! initialize metric
        do i=1,Nx
           do j=1,Ny
            do k=1,Nz

! Kerr-Schild Cartesian compactified coordinates (horizon-penetrating and non-rotating at the boundary): t0,x0,y0,z0

                 x0=x(i)
                 y0=y(j)
                 z0=z(k)

! Kerr-Schild spherical compactified coordinates (horizon-penetrating and non-rotating at the boundary):t0,rho0,theta0,phi0

                 rho0=sqrt(x0**2+y0**2+z0**2)
                 if (rho0.ne.0.0d0) then
                  theta0=acos(x0/rho0)
                 end if
                 if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
                  phi0=atan2(z0,y0)
                  if (phi0.lt.0) phi0=phi0+2*PI
                 end if

            !Kerr-AdS barred quantities if we evolve a perturbation of Kerr-AdS
            if (kerrads_background.eq.1) then
              gb_tt(i,j,k)=0
              gb_tx(i,j,k)=0
              gb_ty(i,j,k)=0
              gb_tz(i,j,k)=0
              gb_xx(i,j,k)=0
              gb_xy(i,j,k)=0
              gb_xz(i,j,k)=0
              gb_yy(i,j,k)=0
              gb_yz(i,j,k)=0
              gb_zz(i,j,k)=0
              Hb_t(i,j,k)=0
              Hb_x(i,j,k)=0
              Hb_y(i,j,k)=0
              Hb_z(i,j,k)=0

            !Kerr-AdS barred quantities if we evolve a perturbation of pure AdS
            else if (kerrads_background.eq.0) then

              if (chr(i,j,k).eq.ex) then
                 gb_tt(i,j,k)=0
                 gb_tx(i,j,k)=0
                 gb_ty(i,j,k)=0
                 gb_tz(i,j,k)=0
                 gb_xx(i,j,k)=0
                 gb_xy(i,j,k)=0
                 gb_xz(i,j,k)=0
                 gb_yy(i,j,k)=0
                 gb_yz(i,j,k)=0
                 gb_zz(i,j,k)=0
                 Hb_t(i,j,k)=0
                 Hb_x(i,j,k)=0
                 Hb_y(i,j,k)=0
                 Hb_z(i,j,k)=0
              else

!KerrAdS in Kerr-Schild coordinates (horizon-penetrating and non-rotating at the boundary) 
! from spherical coords

                if ((a_rot.gt.10.0d0**(-10)).and.
     -             (ief_bh_r0.gt.10.0d0**(-10))) then

                U=Sqrt(2.)/
     -  (L*Sqrt((-(a_rot**2*L**2) + 
     -        (4*L**2*rho0**2)/(-1 + rho0**2)**2 - 
     -        (4*a_rot**2*rho0**2*Sin(theta0)**2)/
     -         (-1 + rho0**2)**2 + 
     -        Sqrt(L**4*
     -           (a_rot**2 + 
     -              (4*rho0**2)/(-1 + rho0**2)**2)**2 + 
     -          (8*a_rot**2*L**2*rho0**2*
     -             (a_rot**2 - 2*L**2 - 
     -               (4*rho0**2)/(-1 + rho0**2)**2)*
     -             Sin(theta0)**2)/(-1 + rho0**2)**2 + 
     -          (16*a_rot**4*rho0**4*Sin(theta0)**4)/
     -           (-1 + rho0**2)**4))/
     -      (L**4*(a_rot**2 + 
     -            (4*rho0**2)/(-1 + rho0**2)**2)**2 + 
     -        (8*a_rot**2*L**2*rho0**2*
     -           (a_rot**2 - 2*L**2 - 
     -             (4*rho0**2)/(-1 + rho0**2)**2)*
     -           Sin(theta0)**2)/(-1 + rho0**2)**2 + 
     -        (16*a_rot**4*rho0**4*Sin(theta0)**4)/
     -         (-1 + rho0**2)**4)))

                Sigma=Sqrt(L**4*(a_rot**2 + 
     -        (4*rho0**2)/(-1 + rho0**2)**2)**2 + 
     -    (8*a_rot**2*L**2*rho0**2*
     -       (a_rot**2 - 2*L**2 - 
     -         (4*rho0**2)/(-1 + rho0**2)**2)*
     -       Sin(theta0)**2)/(-1 + rho0**2)**2 + 
     -    (16*a_rot**4*rho0**4*Sin(theta0)**4)/
     -     (-1 + rho0**2)**4)/L**2

                Deltatheta=-(a_rot**2*L**2 - 2*L**4 - 
     -     (4*L**2*rho0**2)/(-1 + rho0**2)**2 + 
     -     (4*a_rot**2*rho0**2*Sin(theta0)**2)/
     -      (-1 + rho0**2)**2 + 
     -     Sqrt(L**4*(a_rot**2 + 
     -           (4*rho0**2)/(-1 + rho0**2)**2)**2 + 
     -       (8*a_rot**2*L**2*rho0**2*
     -          (a_rot**2 - 2*L**2 - 
     -            (4*rho0**2)/(-1 + rho0**2)**2)*
     -          Sin(theta0)**2)/(-1 + rho0**2)**2 + 
     -       (16*a_rot**4*rho0**4*Sin(theta0)**4)/
     -        (-1 + rho0**2)**4))/(2.*L**4)

                Xi=1 - a_rot**2/L**2

        !the following assumes L=1
        h_tt_kerrads_sph0 =
     -   (2*Deltatheta**2*M0)/(U*Xi**2)
        h_trho_kerrads_sph0 =
     -   (8*Sqrt(2.)*Deltatheta*M0*rho0*(1 + rho0**2)*
     -    Sigma*(-a_rot**4 - 8*rho0**2 + 
     -      8*a_rot**2*rho0**2 - a_rot**4*rho0**2 - 
     -      a_rot**4*rho0**4 - 
     -      a_rot**4*rho0**2*Cos(4*theta0) - 
     -      2*Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)) + 
     -      a_rot**2*Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)) + 
     -      a_rot**2*Cos(2*theta0)*
     -       (-2 - 4*rho0**2 - 2*rho0**4 + 
     -         a_rot**2*(1 + rho0**2)**2 - 
     -         Sqrt(a_rot**4 + 16*rho0**4 - 
     -           16*a_rot**2*rho0**4 + 
     -           4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -           4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -            (1 + rho0**4)*Cos(2*theta0) + 
     -           2*a_rot**4*rho0**4*Cos(4*theta0))))*
     -    Sqrt(1/
     -      (-a_rot**2 + (4*rho0**2)/(-1 + rho0**2)**2 + 
     -        Sqrt((a_rot**4 + 16*rho0**4 - 
     -            16*a_rot**2*rho0**4 + 
     -            4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -            4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -             (1 + rho0**4)*Cos(2*theta0) + 
     -            2*a_rot**4*rho0**4*Cos(4*theta0))/
     -          (-1 + rho0**2)**4) - 
     -        (4*a_rot**2*rho0**2*Sin(theta0)**2)/
     -         (-1 + rho0**2)**2)))/
     -  ((-1 + rho0)*(1 + rho0)*U*Xi*
     -    Sqrt((a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0))/
     -      (-1 + rho0**2)**4)*
     -    (2 - a_rot**2 + 2*rho0**4 - a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0)))*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -      a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0))))
        h_ttheta_kerrads_sph0 =
     -   (-8*Sqrt(2.)*a_rot**2*Deltatheta*M0*rho0**2*Sigma*
     -    Sqrt(1/
     -      (-a_rot**2 + (4*rho0**2)/(-1 + rho0**2)**2 + 
     -        Sqrt((a_rot**4 + 16*rho0**4 - 
     -            16*a_rot**2*rho0**4 + 
     -            4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -            4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -             (1 + rho0**4)*Cos(2*theta0) + 
     -            2*a_rot**4*rho0**4*Cos(4*theta0))/
     -          (-1 + rho0**2)**4) - 
     -        (4*a_rot**2*rho0**2*Sin(theta0)**2)/
     -         (-1 + rho0**2)**2))*Sin(2*theta0))/
     -  (U*Xi*Sqrt((a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0))/
     -      (-1 + rho0**2)**4)*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -      a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0))))
        h_tphi_kerrads_sph0 =
     -   -((Deltatheta*M0*
     -      (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -        a_rot**2*rho0**4 + 
     -        2*a_rot**2*rho0**2*Cos(2*theta0) - 
     -        Sqrt(a_rot**4 + 16*rho0**4 - 
     -          16*a_rot**2*rho0**4 + 
     -          4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -          4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -           (1 + rho0**4)*Cos(2*theta0) + 
     -          2*a_rot**4*rho0**4*Cos(4*theta0))))/
     -    (a_rot*(-1 + rho0)**2*(1 + rho0)**2*U*Xi**2))
        h_rhorho_kerrads_sph0 =
     -   (64*M0*rho0**2*(-1 + rho0**2)**4*
     -    (1 + rho0**2)**2*Sigma**2*
     -    (a_rot**4 + 8*rho0**2 - 8*a_rot**2*rho0**2 + 
     -       a_rot**4*rho0**2 + a_rot**4*rho0**4 + 
     -       a_rot**4*rho0**2*Cos(4*theta0) + 
     -       2*Sqrt(a_rot**4 + 16*rho0**4 - 
     -          16*a_rot**2*rho0**4 + 
     -          4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -          4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -           (1 + rho0**4)*Cos(2*theta0) + 
     -          2*a_rot**4*rho0**4*Cos(4*theta0)) - 
     -       a_rot**2*Sqrt(a_rot**4 + 16*rho0**4 - 
     -          16*a_rot**2*rho0**4 + 
     -          4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -          4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -           (1 + rho0**4)*Cos(2*theta0) + 
     -          2*a_rot**4*rho0**4*Cos(4*theta0)) - 
     -       a_rot**2*Cos(2*theta0)*
     -        (-2 - 4*rho0**2 - 2*rho0**4 + 
     -          a_rot**2*(1 + rho0**2)**2 - 
     -          Sqrt(a_rot**4 + 16*rho0**4 - 
     -            16*a_rot**2*rho0**4 + 
     -            4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -            4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -             (1 + rho0**4)*Cos(2*theta0) + 
     -            2*a_rot**4*rho0**4*Cos(4*theta0))))**2)
     -   /(U*(a_rot**4 + 16*rho0**4 - 
     -      16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -      a_rot**4*rho0**8 - 
     -      4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -       (1 + rho0**4)*Cos(2*theta0) + 
     -      2*a_rot**4*rho0**4*Cos(4*theta0))*
     -    (-a_rot**2 + 4*rho0**2 - a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0)))*
     -    (2 - a_rot**2 + 2*rho0**4 - a_rot**2*rho0**4 + 
     -       2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -       Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)))**2*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -       a_rot**2*rho0**4 + 
     -       2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -       Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)))**2)
        h_rhotheta_kerrads_sph0 =
     -   (-64*a_rot**2*M0*rho0**3*(-1 + rho0**2)**5*
     -    (1 + rho0**2)*Sigma**2*
     -    (-a_rot**4 - 8*rho0**2 + 8*a_rot**2*rho0**2 - 
     -      a_rot**4*rho0**2 - a_rot**4*rho0**4 - 
     -      a_rot**4*rho0**2*Cos(4*theta0) - 
     -      2*Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)) + 
     -      a_rot**2*Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)) + 
     -      a_rot**2*Cos(2*theta0)*
     -       (-2 - 4*rho0**2 - 2*rho0**4 + 
     -         a_rot**2*(1 + rho0**2)**2 - 
     -         Sqrt(a_rot**4 + 16*rho0**4 - 
     -           16*a_rot**2*rho0**4 + 
     -           4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -           4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -            (1 + rho0**4)*Cos(2*theta0) + 
     -           2*a_rot**4*rho0**4*Cos(4*theta0))))*
     -    Sin(2*theta0))/
     -  (U*(a_rot**4 + 16*rho0**4 - 16*a_rot**2*rho0**4 + 
     -      4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -      4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -       (1 + rho0**4)*Cos(2*theta0) + 
     -      2*a_rot**4*rho0**4*Cos(4*theta0))*
     -    (-a_rot**2 + 4*rho0**2 - a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0)))*
     -    (2 - a_rot**2 + 2*rho0**4 - a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0)))*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -       a_rot**2*rho0**4 + 
     -       2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -       Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)))**2)
        h_rhophi_kerrads_sph0 =
     -   (4*Sqrt(2.)*M0*rho0*(1 + rho0**2)*Sigma*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -      a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) - 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0)))*
     -    (-a_rot**4 - 8*rho0**2 + 8*a_rot**2*rho0**2 - 
     -      a_rot**4*rho0**2 - a_rot**4*rho0**4 - 
     -      a_rot**4*rho0**2*Cos(4*theta0) - 
     -      2*Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)) + 
     -      a_rot**2*Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)) + 
     -      a_rot**2*Cos(2*theta0)*
     -       (-2 - 4*rho0**2 - 2*rho0**4 + 
     -         a_rot**2*(1 + rho0**2)**2 - 
     -         Sqrt(a_rot**4 + 16*rho0**4 - 
     -           16*a_rot**2*rho0**4 + 
     -           4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -           4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -            (1 + rho0**4)*Cos(2*theta0) + 
     -           2*a_rot**4*rho0**4*Cos(4*theta0))))*
     -    Sqrt(1/
     -      (-a_rot**2 + (4*rho0**2)/(-1 + rho0**2)**2 + 
     -        Sqrt((a_rot**4 + 16*rho0**4 - 
     -            16*a_rot**2*rho0**4 + 
     -            4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -            4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -             (1 + rho0**4)*Cos(2*theta0) + 
     -            2*a_rot**4*rho0**4*Cos(4*theta0))/
     -          (-1 + rho0**2)**4) - 
     -        (4*a_rot**2*rho0**2*Sin(theta0)**2)/
     -         (-1 + rho0**2)**2)))/
     -  (a_rot*(-1 + rho0**2)**3*U*Xi*
     -    Sqrt((a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0))/
     -      (-1 + rho0**2)**4)*
     -    (-2 + a_rot**2 - 2*rho0**4 + a_rot**2*rho0**4 - 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) - 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0)))*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -      a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0))))
        h_thetatheta_kerrads_sph0 =
     -   (64*a_rot**4*M0*rho0**4*(-1 + rho0**2)**6*
     -    Sigma**2*Sin(2*theta0)**2)/
     -  (U*(a_rot**4 + 16*rho0**4 - 16*a_rot**2*rho0**4 + 
     -      4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -      4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -       (1 + rho0**4)*Cos(2*theta0) + 
     -      2*a_rot**4*rho0**4*Cos(4*theta0))*
     -    (-a_rot**2 + 4*rho0**2 - a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0)))*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -       a_rot**2*rho0**4 + 
     -       2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -       Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)))**2)
        h_thetaphi_kerrads_sph0 =
     -   (4*Sqrt(2.)*a_rot*M0*rho0**2*Sigma*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -      a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) - 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0)))*
     -    Sqrt(1/
     -      (-a_rot**2 + (4*rho0**2)/(-1 + rho0**2)**2 + 
     -        Sqrt((a_rot**4 + 16*rho0**4 - 
     -            16*a_rot**2*rho0**4 + 
     -            4*a_rot**4*rho0**4 + a_rot**4*rho0**8 - 
     -            4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -             (1 + rho0**4)*Cos(2*theta0) + 
     -            2*a_rot**4*rho0**4*Cos(4*theta0))/
     -          (-1 + rho0**2)**4) - 
     -        (4*a_rot**2*rho0**2*Sin(theta0)**2)/
     -         (-1 + rho0**2)**2))*Sin(2*theta0))/
     -  ((-1 + rho0**2)**2*U*Xi*
     -    Sqrt((a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0))/
     -      (-1 + rho0**2)**4)*
     -    (a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -      a_rot**2*rho0**4 + 
     -      2*a_rot**2*rho0**2*Cos(2*theta0) + 
     -      Sqrt(a_rot**4 + 16*rho0**4 - 
     -        16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -        a_rot**4*rho0**8 - 
     -        4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -         (1 + rho0**4)*Cos(2*theta0) + 
     -        2*a_rot**4*rho0**4*Cos(4*theta0))))
        h_phiphi_kerrads_sph0 =
     -   (M0*(a_rot**2 + 4*rho0**2 - 4*a_rot**2*rho0**2 + 
     -       a_rot**2*rho0**4 + 
     -       2*a_rot**2*rho0**2*Cos(2*theta0) - 
     -       Sqrt(a_rot**4 + 16*rho0**4 - 
     -         16*a_rot**2*rho0**4 + 4*a_rot**4*rho0**4 + 
     -         a_rot**4*rho0**8 - 
     -         4*a_rot**2*(-2 + a_rot**2)*rho0**2*
     -          (1 + rho0**4)*Cos(2*theta0) + 
     -         2*a_rot**4*rho0**4*Cos(4*theta0)))**2)/
     -  (2.*a_rot**2*(-1 + rho0**2)**4*U*Xi**2)

        else if (ief_bh_r0.gt.10.0d0**(-10)) then !Schwarzschild-AdS data

          h_tt_kerrads_sph0 =
     -  M0/rho0 - M0*rho0
          h_trho_kerrads_sph0 =
     -  (2*M0*(1 - 2*rho0**2 + rho0**4))/
     -  (rho0*Sqrt((-1 + rho0**2)**2)*(1 + rho0**2))
          h_ttheta_kerrads_sph0 =
     -  0
          h_tphi_kerrads_sph0 =
     -  0
          h_rhorho_kerrads_sph0 =
     -  (-4*M0*(-1 + rho0**2))/(rho0*(1 + rho0**2)**2)
          h_rhotheta_kerrads_sph0 =
     -  0
          h_rhophi_kerrads_sph0 =
     -  0
          h_thetatheta_kerrads_sph0 =
     -  0
          h_thetaphi_kerrads_sph0 =
     -  0
          h_phiphi_kerrads_sph0 =
     -  0

        else !pure AdS data
          h_tt_kerrads_sph0 =
     -  0
          h_trho_kerrads_sph0 =
     -  0
          h_ttheta_kerrads_sph0 =
     -  0
          h_tphi_kerrads_sph0 =
     -  0
          h_rhorho_kerrads_sph0 =
     -  0
          h_rhophi_kerrads_sph0 =
     -  0
          h_thetatheta_kerrads_sph0 =
     -  0
          h_thetaphi_kerrads_sph0 =
     -  0
          h_phiphi_kerrads_sph0 =
     -  0
        end if

        h_kerrads_ll_sph(1,1)=h_tt_kerrads_sph0
        h_kerrads_ll_sph(1,2)=h_trho_kerrads_sph0
        h_kerrads_ll_sph(1,3)=h_ttheta_kerrads_sph0
        h_kerrads_ll_sph(1,4)=h_tphi_kerrads_sph0
        h_kerrads_ll_sph(2,2)=h_rhorho_kerrads_sph0
        h_kerrads_ll_sph(2,3)=h_rhotheta_kerrads_sph0
        h_kerrads_ll_sph(2,4)=h_rhophi_kerrads_sph0
        h_kerrads_ll_sph(3,3)=h_thetatheta_kerrads_sph0
        h_kerrads_ll_sph(3,4)=h_thetaphi_kerrads_sph0
        h_kerrads_ll_sph(4,4)=h_phiphi_kerrads_sph0

        do a=1,3
          do b=a+1,4
            h_kerrads_ll_sph(b,a)=h_kerrads_ll_sph(a,b)
          end do
        end do

        !define transformation matrix between spherical to Cartesian coordinates, 
        !e.g. dxsph_dxcar(2,3)=drho/dtheta

        dxsph_dxcar(1,1)=1
        dxsph_dxcar(1,2)=0
        dxsph_dxcar(1,3)=0
        dxsph_dxcar(1,4)=0

        dxsph_dxcar(2,1)=0
        dxsph_dxcar(2,2)=x0/rho0
        dxsph_dxcar(2,3)=y0/rho0
        dxsph_dxcar(2,4)=z0/rho0

        dxsph_dxcar(3,1)=0
        dxsph_dxcar(3,2)=-(Sqrt(rho0**2 - x0**2)/rho0**2)
        dxsph_dxcar(3,3)=(x0*y0)/(rho0**2*Sqrt(rho0**2 - x0**2))
        dxsph_dxcar(3,4)=(x0*z0)/(rho0**2*Sqrt(rho0**2 - x0**2))

        dxsph_dxcar(4,1)=0
        dxsph_dxcar(4,2)=0
        dxsph_dxcar(4,3)=-(z0/(y0**2 + z0**2))
        dxsph_dxcar(4,4)=y0/(y0**2 + z0**2)

        do a=1,4
         do b=1,4
          h_kerrads_ll(a,b)=0
          do c=1,4
           do d=1,4
            h_kerrads_ll(a,b)=h_kerrads_ll(a,b)+
     &         dxsph_dxcar(c,a)*dxsph_dxcar(d,b)*h_kerrads_ll_sph(c,d)
           end do
          end do
         end do
        end do

        !some of the dxsph_dxcar diverge at y=z=0 so we need to consider this case separately
       if ((abs(y0).lt.10.0d0**(-10)).and.
     &     (abs(z0).lt.10.0d0**(-10))) then

        h_kerrads_ll(1,3)=0
        h_kerrads_ll(1,4)=0
        h_kerrads_ll(2,3)=0
        h_kerrads_ll(2,4)=0
        h_kerrads_ll(3,3)=0
        h_kerrads_ll(3,4)=0
        h_kerrads_ll(4,4)=0

        !re-impose symmetry on indices
        do a=1,3
         do b=a+1,4
          h_kerrads_ll(b,a)=h_kerrads_ll(a,b)
         end do
        end do

       end if

        gb_tt(i,j,k)=h_kerrads_ll(1,1)
        gb_tx(i,j,k)=h_kerrads_ll(1,2)
        gb_ty(i,j,k)=h_kerrads_ll(1,3)
        gb_tz(i,j,k)=h_kerrads_ll(1,4)
        gb_xx(i,j,k)=h_kerrads_ll(2,2)
        gb_xy(i,j,k)=h_kerrads_ll(2,3)
        gb_xz(i,j,k)=h_kerrads_ll(2,4)
        gb_yy(i,j,k)=h_kerrads_ll(3,3)
        gb_yz(i,j,k)=h_kerrads_ll(3,4)
        gb_zz(i,j,k)=h_kerrads_ll(4,4)


              if (is_nan(gb_tt(i,j,k)).or.is_nan(gb_tx(i,j,k))
     &        .or.is_nan(gb_ty(i,j,k))
     &        .or.is_nan(gb_tz(i,j,k))
     &        .or.is_nan(gb_xx(i,j,k)).or.is_nan(gb_xy(i,j,k))
     &        .or.is_nan(gb_xz(i,j,k))
     &        .or.is_nan(gb_yy(i,j,k))
     &        .or.is_nan(gb_yz(i,j,k)).or.is_nan(gb_zz(i,j,k)) ) then

!      if ( (abs(x0-(-0.9375)).lt.10.0d0**(-10)).and.
!     &  (abs(y0-(-0.0)).lt.10.0d0**(-10)).and.
!     &  (abs(z0-(-0.0)).lt.10.0d0**(-10)) ) then

       write (*,*) 'L,i,j,k,x0,y0,z0,rho0=',
     &      L,i,j,k,x0,y0,z0,rho0
       write (*,*) ' U=',U
       write (*,*) ' Sigma=',Sigma
       write (*,*) ' Deltatheta=',Deltatheta
       write (*,*) ' Xi=',Xi
       do a=1,4
        do b=1,4
         write (*,*) "a,b,h_kerrads_ll_sph(a,b)="
     &                ,a,b,h_kerrads_ll_sph(a,b)
        end do
       end do
       do a=1,4
        do b=1,4
         write (*,*) "a,b,h_kerrads_ll(a,b)="
     &                ,a,b,h_kerrads_ll(a,b)
        end do
       end do
       do a=1,4
        do b=1,4
         write (*,*) "a,b,dxsph_dxcar(a,b)="
     &                ,a,b,dxsph_dxcar(a,b)
        end do
       end do
       write (*,*) ' gb_tt=',gb_tt(i,j,k)
       write (*,*) ' gb_tx=',gb_tx(i,j,k)
       write (*,*) ' gb_ty=',gb_ty(i,j,k)
       write (*,*) ' gb_tz=',gb_tz(i,j,k)
       write (*,*) ' gb_xx=',gb_xx(i,j,k)
       write (*,*) ' gb_xy=',gb_xy(i,j,k)
       write (*,*) ' gb_xz=',gb_xz(i,j,k)
       write (*,*) ' gb_yy=',gb_yy(i,j,k)
       write (*,*) ' gb_yz=',gb_yz(i,j,k)
       write (*,*) ' gb_zz=',gb_zz(i,j,k)
       stop
        end if


              end if !closes condition on chr(i,j,k).eq.ex
             end if !closes condition on kerrads_background.eq.0
            end do
           end do
        end do

        ! y=0 axis regularization
!        call axi_reg_g(gb_tt,gb_tx,gb_ty,
!     &                 gb_xx,gb_xy,gb_yy,psi,tfunction,chr,ex,
!     &                 L,x,y,z,Nx,Ny,Nz,regtype)

!        call axi_reg_Hb(Hb_t,Hb_x,Hb_y,chr,ex,L,x,y,z,Nx,Ny,Nz,regtype)

        return
        end

        