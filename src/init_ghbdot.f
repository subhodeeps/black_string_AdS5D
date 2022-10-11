c----------------------------------------------------------------------
c in polar coordinates t==t, x==rho, y==chi/PI
c
c this routine sets d(gb)dt at t=0
c----------------------------------------------------------------------
        subroutine init_ghbdot(gb_tt_n,gb_tx_n,gb_ty_n,
     &                         gb_tz_n,
     &                         gb_xx_n,gb_xy_n,
     &                         gb_xz_n,
     &                         gb_yy_n,
     &                         gb_yz_n,
     &                         gb_zz_n,gb_tt_t_n,gb_tx_t_n,gb_ty_t_n,
     &                         gb_tz_t_n,
     &                         gb_xx_t_n,gb_xy_t_n,
     &                         gb_xz_t_n,
     &                         gb_yy_t_n,
     &                         gb_yz_t_n,
     &                         gb_zz_t_n,Hb_t_n,Hb_x_n,Hb_y_n,
     &                         Hb_z_n,
     &                         Hb_t_t_n,Hb_x_t_n,Hb_y_t_n,
     &                         Hb_z_t_n,
     &                         L,phys_bdy,x,y,z,dt,chr,ex,
     &                         Nx,Ny,Nz,regtype)
        implicit none
        integer Nx,Ny,Nz
        integer regtype
        integer phys_bdy(6),ghost_width(6)
        real*8 dt,ex,L
        real*8 chr(Nx,Ny,Nz)
        real*8 gb_tt_n(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz)
        real*8 gb_ty_n(Nx,Ny,Nz)
        real*8 gb_tz_n(Nx,Ny,Nz)
        real*8 gb_xx_n(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz)
        real*8 gb_xz_n(Nx,Ny,Nz)
        real*8 gb_yy_n(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz)
        real*8 gb_yz_n(Nx,Ny,Nz)
        real*8 gb_tt_t_n(Nx,Ny,Nz),gb_tx_t_n(Nx,Ny,Nz)
        real*8 gb_ty_t_n(Nx,Ny,Nz),gb_zz_t_n(Nx,Ny,Nz)
        real*8 gb_tz_t_n(Nx,Ny,Nz)
        real*8 gb_xx_t_n(Nx,Ny,Nz),gb_xy_t_n(Nx,Ny,Nz)
        real*8 gb_xz_t_n(Nx,Ny,Nz)
        real*8 gb_yy_t_n(Nx,Ny,Nz)
        real*8 gb_yz_t_n(Nx,Ny,Nz)
        real*8 Hb_t_n(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz)
        real*8 Hb_y_n(Nx,Ny,Nz)
        real*8 Hb_z_n(Nx,Ny,Nz)
        real*8 Hb_t_t_n(Nx,Ny,Nz),Hb_x_t_n(Nx,Ny,Nz)
        real*8 Hb_y_t_n(Nx,Ny,Nz)
        real*8 Hb_z_t_n(Nx,Ny,Nz)
        real*8 x(Nx),y(Ny),z(Nz)

        !--------------------------------------------------------------
        ! the following are first time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------

        real*8 gb_tt_t, gb_tt_x, gb_tt_y, gb_tt_tt
        real*8 gb_tt_z
        real*8 gb_tt_xx,gb_tt_yy,gb_tt_tx,gb_tt_ty
        real*8 gb_tt_xy

        real*8 gb_tx_t, gb_tx_x, gb_tx_y, gb_tx_tt
        real*8 gb_tx_z
        real*8 gb_tx_xx,gb_tx_yy,gb_tx_tx,gb_tx_ty
        real*8 gb_tx_xy

        real*8 gb_ty_t, gb_ty_x, gb_ty_y, gb_ty_tt
        real*8 gb_ty_z
        real*8 gb_ty_xx,gb_ty_yy,gb_ty_tx,gb_ty_ty
        real*8 gb_ty_xy

        real*8 gb_tz_t, gb_tz_x, gb_tz_y, gb_tz_tt
        real*8 gb_tz_z
        real*8 gb_tz_xx,gb_tz_yy,gb_tz_tx,gb_tz_ty
        real*8 gb_tz_xy

        real*8 gb_xx_t, gb_xx_x, gb_xx_y, gb_xx_tt
        real*8 gb_xx_z
        real*8 gb_xx_xx,gb_xx_yy,gb_xx_tx,gb_xx_ty
        real*8 gb_xx_xy

        real*8 gb_xy_t, gb_xy_x, gb_xy_y, gb_xy_tt
        real*8 gb_xy_z
        real*8 gb_xy_xx,gb_xy_yy,gb_xy_tx,gb_xy_ty
        real*8 gb_xy_xy

        real*8 gb_xz_t, gb_xz_x, gb_xz_y, gb_xz_tt
        real*8 gb_xz_z
        real*8 gb_xz_xx,gb_xz_yy,gb_xz_tx,gb_xz_ty
        real*8 gb_xz_xy

        real*8 gb_yy_t, gb_yy_x, gb_yy_y, gb_yy_tt
        real*8 gb_yy_z
        real*8 gb_yy_xx,gb_yy_yy,gb_yy_tx,gb_yy_ty
        real*8 gb_yy_xy

        real*8 gb_yz_t, gb_yz_x, gb_yz_y, gb_yz_tt
        real*8 gb_yz_z
        real*8 gb_yz_xx,gb_yz_yy,gb_yz_tx,gb_yz_ty
        real*8 gb_yz_xy

        real*8 gb_zz_t, gb_zz_x, gb_zz_y, gb_zz_tt
        real*8 gb_zz_z
        real*8 gb_zz_xx,gb_zz_yy,gb_zz_tx,gb_zz_ty
        real*8 gb_zz_xy

        real*8 g0_tt_t, g0_tt_x, g0_tt_y, g0_tt_tt
        real*8 g0_tt_z
        real*8 g0_tt_xx,g0_tt_yy,g0_tt_tx,g0_tt_ty
        real*8 g0_tt_xy

        real*8 g0_tx_t, g0_tx_x, g0_tx_y, g0_tx_tt
        real*8 g0_tx_z
        real*8 g0_tx_xx,g0_tx_yy,g0_tx_tx,g0_tx_ty
        real*8 g0_tx_xy

        real*8 g0_ty_t, g0_ty_x, g0_ty_y, g0_ty_tt
        real*8 g0_ty_z
        real*8 g0_ty_xx,g0_ty_yy,g0_ty_tx,g0_ty_ty
        real*8 g0_ty_xy

        real*8 g0_tz_t, g0_tz_x, g0_tz_y, g0_tz_tt
        real*8 g0_tz_z
        real*8 g0_tz_xx,g0_tz_yy,g0_tz_tx,g0_tz_ty
        real*8 g0_tz_xy

        real*8 g0_xx_t, g0_xx_x, g0_xx_y, g0_xx_tt
        real*8 g0_xx_z
        real*8 g0_xx_xx,g0_xx_yy,g0_xx_tx,g0_xx_ty
        real*8 g0_xx_xy

        real*8 g0_xy_t, g0_xy_x, g0_xy_y, g0_xy_tt
        real*8 g0_xy_z
        real*8 g0_xy_xx,g0_xy_yy,g0_xy_tx,g0_xy_ty
        real*8 g0_xy_xy

        real*8 g0_xz_t, g0_xz_x, g0_xz_y, g0_xz_tt
        real*8 g0_xz_z
        real*8 g0_xz_xx,g0_xz_yy,g0_xz_tx,g0_xz_ty
        real*8 g0_xz_xy

        real*8 g0_yy_t, g0_yy_x, g0_yy_y, g0_yy_tt
        real*8 g0_yy_z
        real*8 g0_yy_xx,g0_yy_yy,g0_yy_tx,g0_yy_ty
        real*8 g0_yy_xy

        real*8 g0_yz_t, g0_yz_x, g0_yz_y, g0_yz_tt
        real*8 g0_yz_z
        real*8 g0_yz_xx,g0_yz_yy,g0_yz_tx,g0_yz_ty
        real*8 g0_yz_xy

        real*8 g0_zz_t, g0_zz_x, g0_zz_y, g0_zz_tt
        real*8 g0_zz_z
        real*8 g0_zz_xx,g0_zz_yy,g0_zz_tx,g0_zz_ty
        real*8 g0_zz_xy

        real*8 H0_t_t,H0_t_x,H0_t_y
        real*8 H0_t_z
        real*8 H0_x_t,H0_x_x,H0_x_y
        real*8 H0_x_z
        real*8 H0_y_t,H0_y_x,H0_y_y
        real*8 H0_y_z
        real*8 H0_z_t,H0_z_x,H0_z_y
        real*8 H0_z_z

        real*8 Hb_t0, Hb_x0, Hb_y0
        real*8 Hb_z0

        real*8 Hb_t_t,Hb_t_x,Hb_t_y
        real*8 Hb_t_z
        real*8 Hb_x_t,Hb_x_x,Hb_x_y
        real*8 Hb_x_z
        real*8 Hb_y_t,Hb_y_x,Hb_y_y
        real*8 Hb_y_z
        real*8 Hb_z_t,Hb_z_x,Hb_z_y
        real*8 Hb_z_z

        real*8 gb_tt0,g0u_tt0,g0_tt0
        real*8 gb_tx0,g0u_tx0,g0_tx0
        real*8 gb_ty0,g0u_ty0,g0_ty0
        real*8 gb_tz0,g0u_tz0,g0_tz0
        real*8 gb_xx0,g0u_xx0,g0_xx0
        real*8 gb_xy0,g0u_xy0,g0_xy0
        real*8 gb_xz0,g0u_xz0,g0_xz0
        real*8 gb_yy0,g0u_yy0,g0_yy0
        real*8 gb_yz0,g0u_yz0,g0_yz0
        real*8 gb_zz0,g0u_zz0,g0_zz0
        real*8 m_g0_det0

        real*8 g0_tt_ads_xx,g0_tt_ads_xy,g0_tt_ads_yy
        real*8 g0_tt_ads_x,g0_tt_ads_y,g0_tt_ads0
        real*8 g0_xx_ads_xx,g0_xx_ads_xy,g0_xx_ads_yy
        real*8 g0_xx_ads_x,g0_xx_ads_y,g0_xx_ads0
        real*8 g0_xy_ads_xx,g0_xy_ads_xy,g0_xy_ads_yy
        real*8 g0_xy_ads_x,g0_xy_ads_y,g0_xy_ads0
        real*8 g0_yy_ads_xx,g0_yy_ads_xy,g0_yy_ads_yy
        real*8 g0_yy_ads_x,g0_yy_ads_y,g0_yy_ads0
        real*8 g0_zz_ads_xx,g0_zz_ads_xy,g0_zz_ads_yy
        real*8 g0_zz_ads_y,g0_zz_ads_x,g0_zz_ads0

        real*8 H0_t0,H0_x0,H0_y0
        real*8 H0_z0

        real*8 x0,y0,z0,rho0

        integer i,j,k,is,ie,js,je,ks,ke
        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz

        logical ltrace,extrap_int_boundaries
        parameter (ltrace=.false.,extrap_int_boundaries=.true.)

        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/

        data gb_tt_t, gb_tt_x, gb_tt_y, gb_tt_tt/0.0,0.0,0.0,0.0/
        data gb_tt_xx,gb_tt_yy,gb_tt_tx,gb_tt_ty/0.0,0.0,0.0,0.0/
        data gb_tt_xy/0.0/

        data gb_tx_t, gb_tx_x, gb_tx_y, gb_tx_tt/0.0,0.0,0.0,0.0/
        data gb_tx_xx,gb_tx_yy,gb_tx_tx,gb_tx_ty/0.0,0.0,0.0,0.0/
        data gb_tx_xy/0.0/

        data gb_ty_t, gb_ty_x, gb_ty_y, gb_ty_tt/0.0,0.0,0.0,0.0/
        data gb_ty_xx,gb_ty_yy,gb_ty_tx,gb_ty_ty/0.0,0.0,0.0,0.0/
        data gb_ty_xy/0.0/

        data gb_tz_t, gb_tz_x, gb_tz_y, gb_tz_tt/0.0,0.0,0.0,0.0/
        data gb_tz_xx,gb_tz_yy,gb_tz_tx,gb_tz_ty/0.0,0.0,0.0,0.0/
        data gb_tz_xy/0.0/

        data gb_xx_t, gb_xx_x, gb_xx_y, gb_xx_tt/0.0,0.0,0.0,0.0/
        data gb_xx_xx,gb_xx_yy,gb_xx_tx,gb_xx_ty/0.0,0.0,0.0,0.0/
        data gb_xx_xy/0.0/

        data gb_xy_t, gb_xy_x, gb_xy_y, gb_xy_tt/0.0,0.0,0.0,0.0/
        data gb_xy_xx,gb_xy_yy,gb_xy_tx,gb_xy_ty/0.0,0.0,0.0,0.0/
        data gb_xy_xy/0.0/

        data gb_xz_t, gb_xz_x, gb_xz_y, gb_xz_tt/0.0,0.0,0.0,0.0/
        data gb_xz_xx,gb_xz_yy,gb_xz_tx,gb_xz_ty/0.0,0.0,0.0,0.0/
        data gb_xz_xy/0.0/

        data gb_yy_t, gb_yy_x, gb_yy_y, gb_yy_tt/0.0,0.0,0.0,0.0/
        data gb_yy_xx,gb_yy_yy,gb_yy_tx,gb_yy_ty/0.0,0.0,0.0,0.0/
        data gb_yy_xy/0.0/

        data gb_yz_t, gb_yz_x, gb_yz_y, gb_yz_tt/0.0,0.0,0.0,0.0/
        data gb_yz_xx,gb_yz_yy,gb_yz_tx,gb_yz_ty/0.0,0.0,0.0,0.0/
        data gb_yz_xy/0.0/

        data gb_zz_t, gb_zz_x, gb_zz_y, gb_zz_tt/0.0,0.0,0.0,0.0/
        data gb_zz_xx,gb_zz_yy,gb_zz_tx,gb_zz_ty/0.0,0.0,0.0,0.0/
        data gb_zz_xy/0.0/

        data g0_tt_t, g0_tt_x, g0_tt_y, g0_tt_tt/0.0,0.0,0.0,0.0/
        data g0_tt_xx,g0_tt_yy,g0_tt_tx,g0_tt_ty/0.0,0.0,0.0,0.0/
        data g0_tt_xy/0.0/

        data g0_tx_t, g0_tx_x, g0_tx_y, g0_tx_tt/0.0,0.0,0.0,0.0/
        data g0_tx_xx,g0_tx_yy,g0_tx_tx,g0_tx_ty/0.0,0.0,0.0,0.0/
        data g0_tx_xy/0.0/

        data g0_ty_t, g0_ty_x, g0_ty_y, g0_ty_tt/0.0,0.0,0.0,0.0/
        data g0_ty_xx,g0_ty_yy,g0_ty_tx,g0_ty_ty/0.0,0.0,0.0,0.0/
        data g0_ty_xy/0.0/

        data g0_tz_t, g0_tz_x, g0_tz_y, g0_tz_tt/0.0,0.0,0.0,0.0/
        data g0_tz_xx,g0_tz_yy,g0_tz_tx,g0_tz_ty/0.0,0.0,0.0,0.0/
        data g0_tz_xy/0.0/

        data g0_xx_t, g0_xx_x, g0_xx_y, g0_xx_tt/0.0,0.0,0.0,0.0/
        data g0_xx_xx,g0_xx_yy,g0_xx_tx,g0_xx_ty/0.0,0.0,0.0,0.0/
        data g0_xx_xy/0.0/

        data g0_xy_t, g0_xy_x, g0_xy_y, g0_xy_tt/0.0,0.0,0.0,0.0/
        data g0_xy_xx,g0_xy_yy,g0_xy_tx,g0_xy_ty/0.0,0.0,0.0,0.0/
        data g0_xy_xy/0.0/

        data g0_xz_t, g0_xz_x, g0_xz_y, g0_xz_tt/0.0,0.0,0.0,0.0/
        data g0_xz_xx,g0_xz_yy,g0_xz_tx,g0_xz_ty/0.0,0.0,0.0,0.0/
        data g0_xz_xy/0.0/

        data g0_yy_t, g0_yy_x, g0_yy_y, g0_yy_tt/0.0,0.0,0.0,0.0/
        data g0_yy_xx,g0_yy_yy,g0_yy_tx,g0_yy_ty/0.0,0.0,0.0,0.0/
        data g0_yy_xy/0.0/

        data g0_yz_t, g0_yz_x, g0_yz_y, g0_yz_tt/0.0,0.0,0.0,0.0/
        data g0_yz_xx,g0_yz_yy,g0_yz_tx,g0_yz_ty/0.0,0.0,0.0,0.0/
        data g0_yz_xy/0.0/

        data g0_zz_t, g0_zz_x, g0_zz_y, g0_zz_tt/0.0,0.0,0.0,0.0/
        data g0_zz_xx,g0_zz_yy,g0_zz_tx,g0_zz_ty/0.0,0.0,0.0,0.0/
        data g0_zz_xy/0.0/

        data H0_t_t,H0_t_x,H0_t_y/0.0,0.0,0.0/
        data H0_x_t,H0_x_x,H0_x_y/0.0,0.0,0.0/
        data H0_y_t,H0_y_x,H0_y_y/0.0,0.0,0.0/
        data H0_z_t,H0_z_x,H0_z_y/0.0,0.0,0.0/

        data Hb_t0, Hb_x0, Hb_y0/0.0,0.0,0.0/
        data Hb_z0/0.0/

        data Hb_t_t,Hb_t_x,Hb_t_y/0.0,0.0,0.0/
        data Hb_x_t,Hb_x_x,Hb_x_y/0.0,0.0,0.0/
        data Hb_y_t,Hb_y_x,Hb_y_y/0.0,0.0,0.0/
        data Hb_z_t,Hb_z_x,Hb_z_y/0.0,0.0,0.0/

        data gb_tt0,g0u_tt0,g0_tt0/0.0,0.0,0.0/
        data gb_tx0,g0u_tx0,g0_tx0/0.0,0.0,0.0/
        data gb_ty0,g0u_ty0,g0_ty0/0.0,0.0,0.0/
        data gb_tz0,g0u_tz0,g0_tz0/0.0,0.0,0.0/
        data gb_xx0,g0u_xx0,g0_xx0/0.0,0.0,0.0/
        data gb_xy0,g0u_xy0,g0_xy0/0.0,0.0,0.0/
        data gb_xz0,g0u_xz0,g0_xz0/0.0,0.0,0.0/
        data gb_yy0,g0u_yy0,g0_yy0/0.0,0.0,0.0/
        data gb_yz0,g0u_yz0,g0_yz0/0.0,0.0,0.0/
        data gb_zz0,g0u_zz0,g0_zz0/0.0,0.0,0.0/
        data m_g0_det0/0.0/

        data g0_tt_ads_xx,g0_tt_ads_xy,g0_tt_ads_yy/0.0,0.0,0.0/
        data g0_tt_ads_x,g0_tt_ads_y,g0_tt_ads0/0.0,0.0,0.0/
        data g0_xx_ads_xx,g0_xx_ads_xy,g0_xx_ads_yy/0.0,0.0,0.0/
        data g0_xx_ads_x,g0_xx_ads_y,g0_xx_ads0/0.0,0.0,0.0/
        data g0_xy_ads_xx,g0_xy_ads_xy,g0_xy_ads_yy/0.0,0.0,0.0/
        data g0_xy_ads_x,g0_xy_ads_y,g0_xy_ads0/0.0,0.0,0.0/
        data g0_yy_ads_xx,g0_yy_ads_xy,g0_yy_ads_yy/0.0,0.0,0.0/
        data g0_yy_ads_x,g0_yy_ads_y,g0_yy_ads0/0.0,0.0,0.0/
        data g0_zz_ads_xx,g0_zz_ads_xy,g0_zz_ads_yy/0.0,0.0,0.0/
        data g0_zz_ads_y,g0_zz_ads_x,g0_zz_ads0/0.0,0.0,0.0/

        data H0_t0,H0_x0,H0_y0/0.0,0.0,0.0/
        data H0_z0/0.0/

        data x0,y0,z0/0.0,0.0,0.0/

        data dx,dy,dz/0.0,0.0,0.0/

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        do i=1,Nx
         do j=1,Ny
          do k=1,Nz
           gb_tt_t_n(i,j,k)=0
           gb_tx_t_n(i,j,k)=0
           gb_ty_t_n(i,j,k)=0
           gb_tz_t_n(i,j,k)=0
           gb_xx_t_n(i,j,k)=0
           gb_xy_t_n(i,j,k)=0
           gb_xz_t_n(i,j,k)=0
           gb_yy_t_n(i,j,k)=0
           gb_yz_t_n(i,j,k)=0
           gb_zz_t_n(i,j,k)=0
           Hb_t_t_n(i,j,k)=0
           Hb_x_t_n(i,j,k)=0
           Hb_y_t_n(i,j,k)=0
           Hb_z_t_n(i,j,k)=0
          end do
         end do
        end do

! if we use is=1,ie=Nx,etc. for some processes, we will get error from derivative stencils if we start with initial data with an excised region (e.g. Schwarzschild initial data). That's because the excised region will pass too close to some points of interior boundaries of bounding boxes, so we will have points at the boundary that have no unexcised neighbouring point to use to compute stencils. In principle we can ignore these errors because those points are taken into account by the neighbouring process anyway. If we use is=2,ie=Nx-1,etc. for all processes, we don't have these errors. The points that we are ignoring are those at the boundary of the entire grid (where grid functions are set to 0 anyway) and points next to interior boundaries of bounding boxes, which are again taken into account by nearby processes. So we can use is=2,ie=Nx-1,etc. for all processes without changing the output of the code and avoiding errors
!        is=1
!        js=1
!        ks=1
!        ie=Nx
!        je=Ny
!        ke=Nz
!        if (phys_bdy(1).eq.1) is=2
!        if (phys_bdy(3).eq.1) js=2
!        if (phys_bdy(5).eq.1) ks=2
!        if (phys_bdy(2).eq.1) ie=Nx-1
!        if (phys_bdy(4).eq.1) je=Ny-1
!        if (phys_bdy(6).eq.1) ke=Nz-1


        is=2
        js=2
        ks=2
        ie=Nx-1
        je=Ny-1
        ke=Nz-1


        do i=is,ie
         do j=js,je
          do k=ks,ke
          if ( chr(i,j,k).ne.ex ) then
           x0=x(i)
           y0=y(j)
           z0=z(k)
           rho0=sqrt(x0**2+y0**2+z0**2)


           ! calculate all needed derivatives 
           call df1_int(gb_tt_n,gb_tt_n,gb_tt_n,gb_tt_t,gb_tt_x,
     &          gb_tt_y,
     &          gb_tt_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_tt')

           call df1_int(gb_tx_n,gb_tx_n,gb_tx_n,gb_tx_t,gb_tx_x,
     &          gb_tx_y,
     &          gb_tx_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_tx')

           call df1_int(gb_ty_n,gb_ty_n,gb_ty_n,gb_ty_t,gb_ty_x,
     &          gb_ty_y,
     &          gb_ty_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_ty')

           call df1_int(gb_tz_n,gb_tz_n,gb_tz_n,gb_tz_t,gb_tz_x,
     &          gb_tz_y,
     &          gb_tz_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_tz')

           call df1_int(gb_xx_n,gb_xx_n,gb_xx_n,gb_xx_t,gb_xx_x,
     &          gb_xx_y,
     &          gb_xx_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_xx')

           call df1_int(gb_xy_n,gb_xy_n,gb_xy_n,gb_xy_t,gb_xy_x,
     &          gb_xy_y,
     &          gb_xy_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_xy')

           call df1_int(gb_xz_n,gb_xz_n,gb_xz_n,gb_xz_t,gb_xz_x,
     &          gb_xz_y,
     &          gb_xz_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_xz')

           call df1_int(gb_yy_n,gb_yy_n,gb_yy_n,gb_yy_t,gb_yy_x,
     &          gb_yy_y,
     &          gb_yy_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_yy')

           call df1_int(gb_yz_n,gb_yz_n,gb_yz_n,gb_yz_t,gb_yz_x,
     &          gb_yz_y,
     &          gb_yz_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_yz')

           call df1_int(gb_zz_n,gb_zz_n,gb_zz_n,gb_zz_t,gb_zz_x,
     &          gb_zz_y,
     &          gb_zz_z,
     &          x,y,z,dt,i,j,k,
     &          chr,ex,Nx,Ny,Nz,'gb_zz')

           ! set time derivatives; zero gb_xx_t,gb_xy_t,gb_xz_t,gb_yy_t,gb_yz_t,gb_zz_t from time-symmetry,
           ! possibly nonzero gb_tt_t,gb_tx_t,gb_ty_t,gb_tz_t from gauge freedom
           gb_tt_t_n(i,j,k)=0
           gb_tx_t_n(i,j,k)=0
           gb_ty_t_n(i,j,k)=0
           gb_tz_t_n(i,j,k)=0
           gb_xx_t_n(i,j,k)=0
           gb_xy_t_n(i,j,k)=0
           gb_xz_t_n(i,j,k)=0
           gb_yy_t_n(i,j,k)=0
           gb_yz_t_n(i,j,k)=0
           gb_zz_t_n(i,j,k)=0
           
          end if
          end do
         end do
        end do

        return
        end
