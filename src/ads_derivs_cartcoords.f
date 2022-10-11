c----------------------------------------------------------------------
c calculates Kerr-AdS metric and its first and second derivatives in the version of 
c the Kerr-Schild coordinates (horizon-penetrating) that is non-rotating at the boundary.
c
c r0 is the radius parameter, i.e. r0=2*M0, where M0 is the BHmass 
c (r0 has no physical meaning, it is NOT the horizon radius)
c----------------------------------------------------------------------
        subroutine ads_derivs_cartcoords(
     &                  gads_ll,gads_uu,
     &                  gads_ll_x,
     &                  gads_uu_x,gads_ll_xx,
     &                  Hads_l,
     &                  gammaads_ull,
     &					phi1ads,
     &					phi1ads_x,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                  calc_der,calc_adv_quant)

        implicit none

        integer Nx,Ny,Nz
        integer i,j,k
        real*8  ief_bh_r0,a_rot,M0,M0_min
        logical calc_der,calc_adv_quant

        logical is_nan

        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L

        integer a,b,c,d,e,f,g,h
        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0
        real*8 f0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------
        ! variables for tensor manipulations
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4),gads_ll_xx(4,4,4,4)
        real*8 Hads_l(4),A_l(4),A_l_x(4,4)
        real*8 gads_ll_sph(4,4),gads_uu_sph(4,4)
        real*8 gads_ll_sph_x(4,4,4),gads_uu_sph_x(4,4,4)
        real*8 gads_ll_sph_xx(4,4,4,4)
        real*8 dxsph_dxcar(4,4)
        real*8 d2xsph_dxcardxcar(4,4,4)
        real*8 d3xsph_dxcardxcardxcar(4,4,4,4)
        real*8 gammaads_ull(4,4,4)
        real*8 boxadsx_u(4)
        real*8 phi1ads,phi1ads_x(4),phi1ads_xx(4,4)
        real*8 gammaads_ull_x(4,4,4,4)
        real*8 riemannads_ulll(4,4,4,4)
        real*8 ricciads_ll(4,4),ricciads_lu(4,4),ricciads
        real*8 einsteinads_ll(4,4),setads_ll(4,4)
        real*8 grad_phi1ads_sq

        real*8 g0_tt_ads_sph_t
        real*8 g0_tt_ads_sph_tt
        real*8 g0_tt_ads_sph_trho
        real*8 g0_tt_ads_sph_ttheta
        real*8 g0_tt_ads_sph_tphi
        real*8 g0_tt_ads_sph_rho
        real*8 g0_tt_ads_sph_rhorho
        real*8 g0_tt_ads_sph_rhotheta
        real*8 g0_tt_ads_sph_theta
        real*8 g0_tt_ads_sph_thetatheta
        real*8 g0_tt_ads_sph_phi
        real*8 g0_tt_ads_sph_rhophi
        real*8 g0_tt_ads_sph_thetaphi
        real*8 g0_tt_ads_sph_phiphi
        real*8 g0_trho_ads_sph_t
        real*8 g0_trho_ads_sph_tt
        real*8 g0_trho_ads_sph_trho
        real*8 g0_trho_ads_sph_ttheta
        real*8 g0_trho_ads_sph_tphi
        real*8 g0_trho_ads_sph_rho
        real*8 g0_trho_ads_sph_rhorho
        real*8 g0_trho_ads_sph_rhotheta
        real*8 g0_trho_ads_sph_theta
        real*8 g0_trho_ads_sph_thetatheta
        real*8 g0_trho_ads_sph_phi
        real*8 g0_trho_ads_sph_rhophi
        real*8 g0_trho_ads_sph_thetaphi
        real*8 g0_trho_ads_sph_phiphi
        real*8 g0_ttheta_ads_sph_t
        real*8 g0_ttheta_ads_sph_tt
        real*8 g0_ttheta_ads_sph_trho
        real*8 g0_ttheta_ads_sph_ttheta
        real*8 g0_ttheta_ads_sph_tphi
        real*8 g0_ttheta_ads_sph_rho
        real*8 g0_ttheta_ads_sph_rhorho
        real*8 g0_ttheta_ads_sph_rhotheta
        real*8 g0_ttheta_ads_sph_theta
        real*8 g0_ttheta_ads_sph_thetatheta
        real*8 g0_ttheta_ads_sph_phi
        real*8 g0_ttheta_ads_sph_rhophi
        real*8 g0_ttheta_ads_sph_thetaphi
        real*8 g0_ttheta_ads_sph_phiphi
        real*8 g0_tphi_ads_sph_t
        real*8 g0_tphi_ads_sph_tt
        real*8 g0_tphi_ads_sph_trho
        real*8 g0_tphi_ads_sph_ttheta
        real*8 g0_tphi_ads_sph_tphi
        real*8 g0_tphi_ads_sph_rho
        real*8 g0_tphi_ads_sph_rhorho
        real*8 g0_tphi_ads_sph_rhotheta
        real*8 g0_tphi_ads_sph_theta
        real*8 g0_tphi_ads_sph_thetatheta
        real*8 g0_tphi_ads_sph_phi
        real*8 g0_tphi_ads_sph_rhophi
        real*8 g0_tphi_ads_sph_thetaphi
        real*8 g0_tphi_ads_sph_phiphi
        real*8 g0_rhorho_ads_sph_t
        real*8 g0_rhorho_ads_sph_tt
        real*8 g0_rhorho_ads_sph_trho
        real*8 g0_rhorho_ads_sph_ttheta
        real*8 g0_rhorho_ads_sph_tphi
        real*8 g0_rhorho_ads_sph_rho
        real*8 g0_rhorho_ads_sph_rhorho
        real*8 g0_rhorho_ads_sph_rhotheta
        real*8 g0_rhorho_ads_sph_theta
        real*8 g0_rhorho_ads_sph_thetatheta
        real*8 g0_rhorho_ads_sph_phi
        real*8 g0_rhorho_ads_sph_rhophi
        real*8 g0_rhorho_ads_sph_thetaphi
        real*8 g0_rhorho_ads_sph_phiphi
        real*8 g0_rhotheta_ads_sph_t
        real*8 g0_rhotheta_ads_sph_tt
        real*8 g0_rhotheta_ads_sph_trho
        real*8 g0_rhotheta_ads_sph_ttheta
        real*8 g0_rhotheta_ads_sph_tphi
        real*8 g0_rhotheta_ads_sph_rho
        real*8 g0_rhotheta_ads_sph_rhorho
        real*8 g0_rhotheta_ads_sph_rhotheta
        real*8 g0_rhotheta_ads_sph_theta
        real*8 g0_rhotheta_ads_sph_thetatheta
        real*8 g0_rhotheta_ads_sph_phi
        real*8 g0_rhotheta_ads_sph_rhophi
        real*8 g0_rhotheta_ads_sph_thetaphi
        real*8 g0_rhotheta_ads_sph_phiphi
        real*8 g0_rhophi_ads_sph_t
        real*8 g0_rhophi_ads_sph_tt
        real*8 g0_rhophi_ads_sph_trho
        real*8 g0_rhophi_ads_sph_ttheta
        real*8 g0_rhophi_ads_sph_tphi
        real*8 g0_rhophi_ads_sph_rho
        real*8 g0_rhophi_ads_sph_rhorho
        real*8 g0_rhophi_ads_sph_rhotheta
        real*8 g0_rhophi_ads_sph_theta
        real*8 g0_rhophi_ads_sph_thetatheta
        real*8 g0_rhophi_ads_sph_phi
        real*8 g0_rhophi_ads_sph_rhophi
        real*8 g0_rhophi_ads_sph_thetaphi
        real*8 g0_rhophi_ads_sph_phiphi
        real*8 g0_thetatheta_ads_sph_t
        real*8 g0_thetatheta_ads_sph_tt
        real*8 g0_thetatheta_ads_sph_trho
        real*8 g0_thetatheta_ads_sph_ttheta
        real*8 g0_thetatheta_ads_sph_tphi
        real*8 g0_thetatheta_ads_sph_rho
        real*8 g0_thetatheta_ads_sph_rhorho
        real*8 g0_thetatheta_ads_sph_rhotheta
        real*8 g0_thetatheta_ads_sph_theta
        real*8 g0_thetatheta_ads_sph_thetatheta
        real*8 g0_thetatheta_ads_sph_phi
        real*8 g0_thetatheta_ads_sph_rhophi
        real*8 g0_thetatheta_ads_sph_thetaphi
        real*8 g0_thetatheta_ads_sph_phiphi
        real*8 g0_thetaphi_ads_sph_t
        real*8 g0_thetaphi_ads_sph_tt
        real*8 g0_thetaphi_ads_sph_trho
        real*8 g0_thetaphi_ads_sph_ttheta
        real*8 g0_thetaphi_ads_sph_tphi
        real*8 g0_thetaphi_ads_sph_rho
        real*8 g0_thetaphi_ads_sph_rhorho
        real*8 g0_thetaphi_ads_sph_rhotheta
        real*8 g0_thetaphi_ads_sph_theta
        real*8 g0_thetaphi_ads_sph_thetatheta
        real*8 g0_thetaphi_ads_sph_phi
        real*8 g0_thetaphi_ads_sph_rhophi
        real*8 g0_thetaphi_ads_sph_thetaphi
        real*8 g0_thetaphi_ads_sph_phiphi
        real*8 g0_phiphi_ads_sph_t
        real*8 g0_phiphi_ads_sph_tt
        real*8 g0_phiphi_ads_sph_trho
        real*8 g0_phiphi_ads_sph_ttheta
        real*8 g0_phiphi_ads_sph_tphi
        real*8 g0_phiphi_ads_sph_rho
        real*8 g0_phiphi_ads_sph_rhorho
        real*8 g0_phiphi_ads_sph_rhotheta
        real*8 g0_phiphi_ads_sph_theta
        real*8 g0_phiphi_ads_sph_thetatheta
        real*8 g0_phiphi_ads_sph_phi
        real*8 g0_phiphi_ads_sph_rhophi
        real*8 g0_phiphi_ads_sph_thetaphi
        real*8 g0_phiphi_ads_sph_phiphi

        real*8 g0_tt_ads_sph0,g0_rhorho_ads_sph0
        real*8 g0_trho_ads_sph0,g0_ttheta_ads_sph0,g0_tphi_ads_sph0
        real*8 g0_rhotheta_ads_sph0,g0_thetatheta_ads_sph0
        real*8 g0_phiphi_ads_sph0
        real*8 g0_rhophi_ads_sph0,g0_thetaphi_ads_sph0

        real*8 detg0_ads_sph0
        real*8 g0u_tt_ads_sph0,g0u_rhorho_ads_sph0
        real*8 g0u_trho_ads_sph0,g0u_ttheta_ads_sph0,g0u_tphi_ads_sph0
        real*8 g0u_rhotheta_ads_sph0,g0u_thetatheta_ads_sph0
        real*8 g0u_phiphi_ads_sph0
        real*8 g0u_rhophi_ads_sph0,g0u_thetaphi_ads_sph0


        real*8 g0_tt_ads_t
        real*8 g0_tt_ads_tt
        real*8 g0_tt_ads_tx
        real*8 g0_tt_ads_ty
        real*8 g0_tt_ads_tz
        real*8 g0_tt_ads_x,g0_tt_ads_xx,g0_tt_ads_xy
        real*8 g0_tt_ads_y,g0_tt_ads_yy
        real*8 g0_tt_ads_z
        real*8 g0_tt_ads_xz
        real*8 g0_tt_ads_yz
        real*8 g0_tt_ads_zz
        real*8 g0_tx_ads_t
        real*8 g0_tx_ads_tt
        real*8 g0_tx_ads_tx
        real*8 g0_tx_ads_ty
        real*8 g0_tx_ads_tz
        real*8 g0_tx_ads_x,g0_tx_ads_xx,g0_tx_ads_xy
        real*8 g0_tx_ads_y,g0_tx_ads_yy
        real*8 g0_tx_ads_z
        real*8 g0_tx_ads_xz
        real*8 g0_tx_ads_yz
        real*8 g0_tx_ads_zz
        real*8 g0_ty_ads_t
        real*8 g0_ty_ads_tt
        real*8 g0_ty_ads_tx
        real*8 g0_ty_ads_ty
        real*8 g0_ty_ads_tz
        real*8 g0_ty_ads_x,g0_ty_ads_xx,g0_ty_ads_xy
        real*8 g0_ty_ads_y,g0_ty_ads_yy
        real*8 g0_ty_ads_z
        real*8 g0_ty_ads_xz
        real*8 g0_ty_ads_yz
        real*8 g0_ty_ads_zz
        real*8 g0_tz_ads_t
        real*8 g0_tz_ads_tt
        real*8 g0_tz_ads_tx
        real*8 g0_tz_ads_ty
        real*8 g0_tz_ads_tz
        real*8 g0_tz_ads_x,g0_tz_ads_xx,g0_tz_ads_xy
        real*8 g0_tz_ads_y,g0_tz_ads_yy
        real*8 g0_tz_ads_z
        real*8 g0_tz_ads_xz
        real*8 g0_tz_ads_yz
        real*8 g0_tz_ads_zz
        real*8 g0_xx_ads_t
        real*8 g0_xx_ads_tt
        real*8 g0_xx_ads_tx
        real*8 g0_xx_ads_ty
        real*8 g0_xx_ads_tz
        real*8 g0_xx_ads_x,g0_xx_ads_xx,g0_xx_ads_xy
        real*8 g0_xx_ads_y,g0_xx_ads_yy
        real*8 g0_xx_ads_z
        real*8 g0_xx_ads_xz
        real*8 g0_xx_ads_yz
        real*8 g0_xx_ads_zz
        real*8 g0_xy_ads_t
        real*8 g0_xy_ads_tt
        real*8 g0_xy_ads_tx
        real*8 g0_xy_ads_ty
        real*8 g0_xy_ads_tz
        real*8 g0_xy_ads_x,g0_xy_ads_xx,g0_xy_ads_xy
        real*8 g0_xy_ads_y,g0_xy_ads_yy
        real*8 g0_xy_ads_z
        real*8 g0_xy_ads_xz
        real*8 g0_xy_ads_yz
        real*8 g0_xy_ads_zz
        real*8 g0_xz_ads_t
        real*8 g0_xz_ads_tt
        real*8 g0_xz_ads_tx
        real*8 g0_xz_ads_ty
        real*8 g0_xz_ads_tz
        real*8 g0_xz_ads_x,g0_xz_ads_xx,g0_xz_ads_xy
        real*8 g0_xz_ads_y,g0_xz_ads_yy
        real*8 g0_xz_ads_z
        real*8 g0_xz_ads_xz
        real*8 g0_xz_ads_yz
        real*8 g0_xz_ads_zz
        real*8 g0_yy_ads_t
        real*8 g0_yy_ads_tt
        real*8 g0_yy_ads_tx
        real*8 g0_yy_ads_ty
        real*8 g0_yy_ads_tz
        real*8 g0_yy_ads_x,g0_yy_ads_xx,g0_yy_ads_xy
        real*8 g0_yy_ads_y,g0_yy_ads_yy
        real*8 g0_yy_ads_z
        real*8 g0_yy_ads_xz
        real*8 g0_yy_ads_yz
        real*8 g0_yy_ads_zz
        real*8 g0_yz_ads_t
        real*8 g0_yz_ads_tt
        real*8 g0_yz_ads_tx
        real*8 g0_yz_ads_ty
        real*8 g0_yz_ads_tz
        real*8 g0_yz_ads_x,g0_yz_ads_xx,g0_yz_ads_xy
        real*8 g0_yz_ads_y,g0_yz_ads_yy
        real*8 g0_yz_ads_z
        real*8 g0_yz_ads_xz
        real*8 g0_yz_ads_yz
        real*8 g0_yz_ads_zz
        real*8 g0_zz_ads_t
        real*8 g0_zz_ads_tt
        real*8 g0_zz_ads_tx
        real*8 g0_zz_ads_ty
        real*8 g0_zz_ads_tz
        real*8 g0_zz_ads_x,g0_zz_ads_xx,g0_zz_ads_xy
        real*8 g0_zz_ads_y,g0_zz_ads_yy
        real*8 g0_zz_ads_z
        real*8 g0_zz_ads_xz
        real*8 g0_zz_ads_yz
        real*8 g0_zz_ads_zz

        real*8 g0_tt_ads0,g0_xx_ads0
        real*8 g0_tx_ads0,g0_ty_ads0,g0_tz_ads0
        real*8 g0_xy_ads0,g0_yy_ads0,g0_zz_ads0
        real*8 g0_xz_ads0,g0_yz_ads0

        real*8 detg0_ads0
        real*8 g0u_tt_ads0,g0u_xx_ads0
        real*8 g0u_tx_ads0,g0u_ty_ads0,g0u_tz_ads0
        real*8 g0u_xy_ads0,g0u_yy_ads0,g0u_zz_ads0
        real*8 g0u_xz_ads0,g0u_yz_ads0

!!!!!!!!!!DEBUG DERIVATIVE STENCILS!!!!!!!!!!!
        real*8 testf1(Nx,Ny,Nz),testf2(Nx,Ny,Nz),testf3(Nx,Ny,Nz)
        real*8 testf1_t,testf1_x,testf1_y,testf1_z
        real*8 testf2_t,testf2_x,testf2_y,testf2_z
        real*8 testf3_t,testf3_x,testf3_y,testf3_z
        real*8 testf1_tt,testf1_tx,testf1_ty
        real*8 testf1_xx,testf1_xy,testf1_yy
        real*8 testf1_tz,testf1_xz,testf1_yz,testf1_zz
        real*8 testf2_tt,testf2_tx,testf2_ty
        real*8 testf2_xx,testf2_xy,testf2_yy
        real*8 testf2_tz,testf2_xz,testf2_yz,testf2_zz
        real*8 testf3_tt,testf3_tx,testf3_ty
        real*8 testf3_xx,testf3_xy,testf3_yy
        real*8 testf3_tz,testf3_xz,testf3_yz,testf3_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------

        if ((calc_adv_quant).and.(.not.calc_der)) then
         write (*,*) "Error: cannot compute 
     -    advanced tensorial quantities 
     -    without computing metric derivatives"
         stop
        end if
        
        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        x0=x(i)
        y0=y(j)
        z0=z(k)
        rho0=sqrt(x0**2+y0**2+z0**2)
        if (rho0.ne.0.0d0) then
          theta0=acos(x0/rho0)
        end if
        if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
          phi0=atan2(z0,y0)
          if (phi0.lt.0) phi0=phi0+2*PI
        end if


!!DEBUG!!!
!
!        g0_tt_ads0 =0
!        g0_tx_ads0 =0
!        g0_ty_ads0 =0
!        g0_tz_ads0 =0
!        g0_xx_ads0 =0
!        g0_xy_ads0 =0
!        g0_xz_ads0 =0
!        g0_yy_ads0 =0
!        g0_yz_ads0 =0
!        g0_zz_ads0=0
!
!        g0u_tt_ads0 =0
!        g0u_tx_ads0 =0
!        g0u_ty_ads0 =0
!        g0u_tz_ads0 =0
!        g0u_xx_ads0 =0
!        g0u_xy_ads0 =0
!        g0u_xz_ads0 =0
!        g0u_yy_ads0 =0
!        g0u_yz_ads0 =0
!        g0u_zz_ads0=0
!
!        g0_tt_ads_x  =0
!        g0_tt_ads_y  =0
!        g0_tt_ads_z  =0
!        g0_tt_ads_xx =0
!        g0_tt_ads_xy =0
!        g0_tt_ads_xz =0
!        g0_tt_ads_yy =0
!        g0_tt_ads_yz =0
!        g0_tt_ads_zz =0
!
!        g0_tx_ads_x  =0
!        g0_tx_ads_y  =0
!        g0_tx_ads_z  =0
!        g0_tx_ads_xx =0
!        g0_tx_ads_xy =0
!        g0_tx_ads_xz =0
!        g0_tx_ads_yy =0
!        g0_tx_ads_yz =0
!        g0_tx_ads_zz =0
!
!        g0_ty_ads_x  =0
!        g0_ty_ads_y  =0
!        g0_ty_ads_z  =0
!        g0_ty_ads_xx =0
!        g0_ty_ads_xy =0
!        g0_ty_ads_xz =0
!        g0_ty_ads_yy =0
!        g0_ty_ads_yz =0
!        g0_ty_ads_zz =0
!
!        g0_tz_ads_x  =0
!        g0_tz_ads_y  =0
!        g0_tz_ads_z  =0
!        g0_tz_ads_xx =0
!        g0_tz_ads_xy =0
!        g0_tz_ads_xz =0
!        g0_tz_ads_yy =0
!        g0_tz_ads_yz =0
!        g0_tz_ads_zz =0
!
!        g0_xx_ads_x  =0
!        g0_xx_ads_y  =0
!        g0_xx_ads_z  =0
!        g0_xx_ads_xx =0
!        g0_xx_ads_xy =0
!        g0_xx_ads_xz =0
!        g0_xx_ads_yy =0
!        g0_xx_ads_yz =0
!        g0_xx_ads_zz =0
!
!        g0_xy_ads_x  =0
!        g0_xy_ads_y  =0
!        g0_xy_ads_z  =0
!        g0_xy_ads_xx =0
!        g0_xy_ads_xy =0
!        g0_xy_ads_xz =0
!        g0_xy_ads_yy =0
!        g0_xy_ads_yz =0
!        g0_xy_ads_zz =0
!
!        g0_xz_ads_x  =0
!        g0_xz_ads_y  =0
!        g0_xz_ads_z  =0
!        g0_xz_ads_xx =0
!        g0_xz_ads_xy =0
!        g0_xz_ads_xz =0
!        g0_xz_ads_yy =0
!        g0_xz_ads_yz =0
!        g0_xz_ads_zz =0
!
!        g0_yy_ads_x  =0
!        g0_yy_ads_y  =0
!        g0_yy_ads_z  =0
!        g0_yy_ads_xx =0
!        g0_yy_ads_xy =0
!        g0_yy_ads_xz =0
!        g0_yy_ads_yy =0
!        g0_yy_ads_yz =0
!        g0_yy_ads_zz =0
!
!        g0_yz_ads_x  =0
!        g0_yz_ads_y  =0
!        g0_yz_ads_z  =0
!        g0_yz_ads_xx =0
!        g0_yz_ads_xy =0
!        g0_yz_ads_xz =0
!        g0_yz_ads_yy =0
!        g0_yz_ads_yz =0
!        g0_yz_ads_zz =0
!
!        g0_zz_ads_x  =0
!        g0_zz_ads_y  =0
!        g0_zz_ads_z  =0
!        g0_zz_ads_xx =0
!        g0_zz_ads_xy =0
!        g0_zz_ads_xz =0
!        g0_zz_ads_yy =0
!        g0_zz_ads_yz =0
!        g0_zz_ads_zz =0
!
!        Hads_l(1)=0
!        Hads_l(2)=0
!        Hads_l(3)=0
!        Hads_l(4)=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        g0_tt_ads_sph0 =
     -   -((4*rho0**2 + L**2*(-1 + rho0**2)**2)/
     -    (L**2*(-1 + rho0**2)**2))
        g0_trho_ads_sph0 =
     -   0
        g0_ttheta_ads_sph0 =
     -   0
        g0_tphi_ads_sph0 =
     -   0
        g0_rhorho_ads_sph0 =
     -   (4*L**2*(1 + rho0**2)**2)/
     -  ((-1 + rho0**2)**2*
     -    (4*rho0**2 + L**2*(-1 + rho0**2)**2))
        g0_rhotheta_ads_sph0 =
     -   0
        g0_rhophi_ads_sph0 =
     -   0
        g0_thetatheta_ads_sph0 =
     -   (4*rho0**2)/(-1 + rho0**2)**2
        g0_thetaphi_ads_sph0 =
     -   0
        g0_phiphi_ads_sph0 =
     -   (4*rho0**2*Sin(theta0)**2)/(-1 + rho0**2)**2

             ! give values to the metric inverse
        call calc_g0uu(g0_tt_ads_sph0,g0_trho_ads_sph0,
     &         g0_ttheta_ads_sph0,g0_tphi_ads_sph0,
     &         g0_rhorho_ads_sph0,g0_rhotheta_ads_sph0,
     &         g0_rhophi_ads_sph0,
     &         g0_thetatheta_ads_sph0,g0_thetaphi_ads_sph0,
     &         g0_phiphi_ads_sph0,
     &         g0u_tt_ads_sph0,g0u_trho_ads_sph0,
     &         g0u_ttheta_ads_sph0,g0u_tphi_ads_sph0,
     &         g0u_rhorho_ads_sph0,g0u_rhotheta_ads_sph0,
     &         g0u_rhophi_ads_sph0,
     &         g0u_thetatheta_ads_sph0,g0u_thetaphi_ads_sph0,
     &         g0u_phiphi_ads_sph0,detg0_ads_sph0)

       if (calc_der) then
        g0_tt_ads_sph_t  =
     -   0
        g0_tt_ads_sph_rho  =
     -   (8*(rho0 + rho0**3))/(L**2*(-1 + rho0**2)**3)
        g0_tt_ads_sph_theta  =
     -   0
        g0_tt_ads_sph_phi  =
     -   0
        g0_tt_ads_sph_tt =
     -   0
        g0_tt_ads_sph_trho =
     -   0
        g0_tt_ads_sph_ttheta =
     -   0
        g0_tt_ads_sph_tphi =
     -   0
        g0_tt_ads_sph_rhorho =
     -    (-8*(1 + 8*rho0**2 + 3*rho0**4))/
     -  (L**2*(-1 + rho0**2)**4)
        g0_tt_ads_sph_rhotheta =
     -   0
        g0_tt_ads_sph_rhophi =
     -   0
        g0_tt_ads_sph_thetatheta =
     -   0
        g0_tt_ads_sph_thetaphi =
     -   0
        g0_tt_ads_sph_phiphi =
     -   0

        g0_trho_ads_sph_t  =
     -   0
        g0_trho_ads_sph_rho  =
     -   0
        g0_trho_ads_sph_theta  =
     -   0
        g0_trho_ads_sph_phi  =
     -   0
        g0_trho_ads_sph_tt =
     -   0
        g0_trho_ads_sph_trho =
     -   0
        g0_trho_ads_sph_ttheta =
     -   0
        g0_trho_ads_sph_tphi =
     -   0
        g0_trho_ads_sph_rhorho =
     -   0
        g0_trho_ads_sph_rhotheta =
     -   0
        g0_trho_ads_sph_rhophi =
     -   0
        g0_trho_ads_sph_thetatheta =
     -   0
        g0_trho_ads_sph_thetaphi =
     -   0
        g0_trho_ads_sph_phiphi =
     -   0

        g0_ttheta_ads_sph_t  =
     -   0
        g0_ttheta_ads_sph_rho  =
     -   0
        g0_ttheta_ads_sph_theta  =
     -   0
        g0_ttheta_ads_sph_phi  =
     -   0
        g0_ttheta_ads_sph_tt =
     -   0
        g0_ttheta_ads_sph_trho =
     -   0
        g0_ttheta_ads_sph_ttheta =
     -   0
        g0_ttheta_ads_sph_tphi =
     -   0
        g0_ttheta_ads_sph_rhorho =
     -   0
        g0_ttheta_ads_sph_rhotheta =
     -   0
        g0_ttheta_ads_sph_rhophi =
     -   0
        g0_ttheta_ads_sph_thetatheta =
     -   0
        g0_ttheta_ads_sph_thetaphi =
     -   0
        g0_ttheta_ads_sph_phiphi =
     -   0

        g0_tphi_ads_sph_t  =
     -   0
        g0_tphi_ads_sph_rho  =
     -   0
        g0_tphi_ads_sph_theta  =
     -   0
        g0_tphi_ads_sph_phi  =
     -   0
        g0_tphi_ads_sph_tt =
     -   0
        g0_tphi_ads_sph_trho =
     -   0
        g0_tphi_ads_sph_ttheta =
     -   0
        g0_tphi_ads_sph_tphi =
     -   0
        g0_tphi_ads_sph_rhorho =
     -   0
        g0_tphi_ads_sph_rhotheta =
     -   0
        g0_tphi_ads_sph_rhophi =
     -   0
        g0_tphi_ads_sph_thetatheta =
     -   0
        g0_tphi_ads_sph_thetaphi =
     -   0
        g0_tphi_ads_sph_phiphi =
     -   0

        g0_rhorho_ads_sph_t  =
     -   0
        g0_rhorho_ads_sph_rho  =
     -   (-16*L**2*rho0*(1 + rho0**2)*
     -    (-2 + 8*rho0**2 + 2*rho0**4 + 
     -      L**2*(-1 + rho0**2)**2*(3 + rho0**2)))/
     -  ((-1 + rho0**2)**3*
     -    (4*rho0**2 + L**2*(-1 + rho0**2)**2)**2)
        g0_rhorho_ads_sph_theta  =
     -   0
        g0_rhorho_ads_sph_phi  =
     -   0
        g0_rhorho_ads_sph_tt =
     -   0
        g0_rhorho_ads_sph_trho =
     -   0
        g0_rhorho_ads_sph_ttheta =
     -   0
        g0_rhorho_ads_sph_tphi =
     -   0
        g0_rhorho_ads_sph_rhorho =
     -   (16*(L**6*(-1 + rho0**2)**4*
     -       (3 + 39*rho0**2 + 33*rho0**4 + 5*rho0**6)
     -       + 8*L**2*rho0**2*
     -       (3 - 12*rho0**2 + 26*rho0**4 + 
     -         28*rho0**6 + 3*rho0**8) + 
     -      2*L**4*(-1 + rho0**2)**2*
     -       (-1 - 22*rho0**2 + 80*rho0**4 + 
     -         78*rho0**6 + 9*rho0**8)))/
     -  ((-1 + rho0**2)**4*
     -    (4*rho0**2 + L**2*(-1 + rho0**2)**2)**3)
        g0_rhorho_ads_sph_rhotheta =
     -   0
        g0_rhorho_ads_sph_rhophi =
     -   0
        g0_rhorho_ads_sph_thetatheta =
     -   0
        g0_rhorho_ads_sph_thetaphi =
     -   0
        g0_rhorho_ads_sph_phiphi =
     -   0

        g0_rhotheta_ads_sph_t  =
     -   0
        g0_rhotheta_ads_sph_rho  =
     -   0
        g0_rhotheta_ads_sph_theta  =
     -   0
        g0_rhotheta_ads_sph_phi  =
     -   0
        g0_rhotheta_ads_sph_tt =
     -   0
        g0_rhotheta_ads_sph_trho =
     -   0
        g0_rhotheta_ads_sph_ttheta =
     -   0
        g0_rhotheta_ads_sph_tphi =
     -   0
        g0_rhotheta_ads_sph_rhorho =
     -   0
        g0_rhotheta_ads_sph_rhotheta =
     -   0
        g0_rhotheta_ads_sph_rhophi =
     -   0
        g0_rhotheta_ads_sph_thetatheta =
     -   0
        g0_rhotheta_ads_sph_thetaphi =
     -   0
        g0_rhotheta_ads_sph_phiphi =
     -   0

        g0_rhophi_ads_sph_t  =
     -   0
        g0_rhophi_ads_sph_rho  =
     -   0
        g0_rhophi_ads_sph_theta  =
     -   0
        g0_rhophi_ads_sph_phi  =
     -   0
        g0_rhophi_ads_sph_tt =
     -   0
        g0_rhophi_ads_sph_trho =
     -   0
        g0_rhophi_ads_sph_ttheta =
     -   0
        g0_rhophi_ads_sph_tphi =
     -   0
        g0_rhophi_ads_sph_rhorho =
     -   0
        g0_rhophi_ads_sph_rhotheta =
     -   0
        g0_rhophi_ads_sph_rhophi =
     -   0
        g0_rhophi_ads_sph_thetatheta =
     -   0
        g0_rhophi_ads_sph_thetaphi =
     -   0
        g0_rhophi_ads_sph_phiphi =
     -   0

        g0_thetatheta_ads_sph_t  =
     -   0
        g0_thetatheta_ads_sph_rho  =
     -   (-8*(rho0 + rho0**3))/(-1 + rho0**2)**3
        g0_thetatheta_ads_sph_theta  =
     -   0
        g0_thetatheta_ads_sph_phi  =
     -   0
        g0_thetatheta_ads_sph_tt =
     -   0
        g0_thetatheta_ads_sph_trho =
     -   0
        g0_thetatheta_ads_sph_ttheta =
     -   0
        g0_thetatheta_ads_sph_tphi =
     -   0
        g0_thetatheta_ads_sph_rhorho =
     -   (8*(1 + 8*rho0**2 + 3*rho0**4))/(-1 + rho0**2)**4
        g0_thetatheta_ads_sph_rhotheta =
     -   0
        g0_thetatheta_ads_sph_rhophi =
     -   0
        g0_thetatheta_ads_sph_thetatheta =
     -   0
        g0_thetatheta_ads_sph_thetaphi =
     -   0
        g0_thetatheta_ads_sph_phiphi =
     -   0

        g0_thetaphi_ads_sph_t  =
     -   0
        g0_thetaphi_ads_sph_rho  =
     -   0
        g0_thetaphi_ads_sph_theta  =
     -   0
        g0_thetaphi_ads_sph_phi  =
     -   0
        g0_thetaphi_ads_sph_tt =
     -   0
        g0_thetaphi_ads_sph_trho =
     -   0
        g0_thetaphi_ads_sph_ttheta =
     -   0
        g0_thetaphi_ads_sph_tphi =
     -   0
        g0_thetaphi_ads_sph_rhorho =
     -   0
        g0_thetaphi_ads_sph_rhotheta =
     -   0
        g0_thetaphi_ads_sph_rhophi =
     -   0
        g0_thetaphi_ads_sph_thetatheta =
     -   0
        g0_thetaphi_ads_sph_thetaphi =
     -   0
        g0_thetaphi_ads_sph_phiphi =
     -   0

        g0_phiphi_ads_sph_t  =
     -   0
        g0_phiphi_ads_sph_rho  =
     -   (-8*(rho0 + rho0**3)*Sin(theta0)**2)/(-1 + rho0**2)**3
        g0_phiphi_ads_sph_theta  =
     -   (4*rho0**2*Sin(2*theta0))/(-1 + rho0**2)**2
        g0_phiphi_ads_sph_phi  =
     -   0
        g0_phiphi_ads_sph_tt =
     -   0
        g0_phiphi_ads_sph_trho =
     -   0
        g0_phiphi_ads_sph_ttheta =
     -   0
        g0_phiphi_ads_sph_tphi =
     -   0
        g0_phiphi_ads_sph_rhorho =
     -   (8*(1 + 8*rho0**2 + 3*rho0**4)*Sin(theta0)**2)/
     -  (-1 + rho0**2)**4
        g0_phiphi_ads_sph_rhotheta =
     -   (-8*(rho0 + rho0**3)*Sin(2*theta0))/(-1 + rho0**2)**3
        g0_phiphi_ads_sph_rhophi =
     -   0
        g0_phiphi_ads_sph_thetatheta =
     -   (8*rho0**2*Cos(2*theta0))/(-1 + rho0**2)**2
        g0_phiphi_ads_sph_thetaphi =
     -   0
        g0_phiphi_ads_sph_phiphi =
     -   0     

       end if

        ! give values to the AdS metric
        gads_ll_sph(1,1)=g0_tt_ads_sph0
        gads_ll_sph(1,2)=g0_trho_ads_sph0
        gads_ll_sph(1,3)=g0_ttheta_ads_sph0
        gads_ll_sph(1,4)=g0_tphi_ads_sph0
        gads_ll_sph(2,2)=g0_rhorho_ads_sph0
        gads_ll_sph(2,3)=g0_rhotheta_ads_sph0
        gads_ll_sph(2,4)=g0_rhophi_ads_sph0
        gads_ll_sph(3,3)=g0_thetatheta_ads_sph0
        gads_ll_sph(3,4)=g0_thetaphi_ads_sph0
        gads_ll_sph(4,4)=g0_phiphi_ads_sph0


        gads_uu_sph(1,1)=g0u_tt_ads_sph0
        gads_uu_sph(1,2)=g0u_trho_ads_sph0
        gads_uu_sph(1,3)=g0u_ttheta_ads_sph0
        gads_uu_sph(1,4)=g0u_tphi_ads_sph0
        gads_uu_sph(2,2)=g0u_rhorho_ads_sph0
        gads_uu_sph(2,3)=g0u_rhotheta_ads_sph0
        gads_uu_sph(2,4)=g0u_rhophi_ads_sph0
        gads_uu_sph(3,3)=g0u_thetatheta_ads_sph0
        gads_uu_sph(3,4)=g0u_thetaphi_ads_sph0
        gads_uu_sph(4,4)=g0u_phiphi_ads_sph0

       if (calc_der) then
        gads_ll_sph_x(1,1,1)   =g0_tt_ads_sph_t
        gads_ll_sph_x(1,1,2)   =g0_tt_ads_sph_rho
        gads_ll_sph_x(1,1,3)   =g0_tt_ads_sph_theta
        gads_ll_sph_x(1,1,4)   =g0_tt_ads_sph_phi
        gads_ll_sph_xx(1,1,1,1)=g0_tt_ads_sph_tt
        gads_ll_sph_xx(1,1,1,2)=g0_tt_ads_sph_trho
        gads_ll_sph_xx(1,1,1,3)=g0_tt_ads_sph_ttheta
        gads_ll_sph_xx(1,1,1,4)=g0_tt_ads_sph_tphi
        gads_ll_sph_xx(1,1,2,2)=g0_tt_ads_sph_rhorho
        gads_ll_sph_xx(1,1,2,3)=g0_tt_ads_sph_rhotheta
        gads_ll_sph_xx(1,1,2,4)=g0_tt_ads_sph_rhophi
        gads_ll_sph_xx(1,1,3,3)=g0_tt_ads_sph_thetatheta
        gads_ll_sph_xx(1,1,3,4)=g0_tt_ads_sph_thetaphi
        gads_ll_sph_xx(1,1,4,4)=g0_tt_ads_sph_phiphi

        gads_ll_sph_x(1,2,1)   =g0_trho_ads_sph_t
        gads_ll_sph_x(1,2,2)   =g0_trho_ads_sph_rho
        gads_ll_sph_x(1,2,3)   =g0_trho_ads_sph_theta
        gads_ll_sph_x(1,2,4)   =g0_trho_ads_sph_phi
        gads_ll_sph_xx(1,2,1,1)=g0_trho_ads_sph_tt
        gads_ll_sph_xx(1,2,1,2)=g0_trho_ads_sph_trho
        gads_ll_sph_xx(1,2,1,3)=g0_trho_ads_sph_ttheta
        gads_ll_sph_xx(1,2,1,4)=g0_trho_ads_sph_tphi
        gads_ll_sph_xx(1,2,2,2)=g0_trho_ads_sph_rhorho
        gads_ll_sph_xx(1,2,2,3)=g0_trho_ads_sph_rhotheta
        gads_ll_sph_xx(1,2,2,4)=g0_trho_ads_sph_rhophi
        gads_ll_sph_xx(1,2,3,3)=g0_trho_ads_sph_thetatheta
        gads_ll_sph_xx(1,2,3,4)=g0_trho_ads_sph_thetaphi
        gads_ll_sph_xx(1,2,4,4)=g0_trho_ads_sph_phiphi

        gads_ll_sph_x(1,3,1)   =g0_ttheta_ads_sph_t
        gads_ll_sph_x(1,3,2)   =g0_ttheta_ads_sph_rho
        gads_ll_sph_x(1,3,3)   =g0_ttheta_ads_sph_theta
        gads_ll_sph_x(1,3,4)   =g0_ttheta_ads_sph_phi
        gads_ll_sph_xx(1,3,1,1)=g0_ttheta_ads_sph_tt
        gads_ll_sph_xx(1,3,1,2)=g0_ttheta_ads_sph_trho
        gads_ll_sph_xx(1,3,1,3)=g0_ttheta_ads_sph_ttheta
        gads_ll_sph_xx(1,3,1,4)=g0_ttheta_ads_sph_tphi
        gads_ll_sph_xx(1,3,2,2)=g0_ttheta_ads_sph_rhorho
        gads_ll_sph_xx(1,3,2,3)=g0_ttheta_ads_sph_rhotheta
        gads_ll_sph_xx(1,3,2,4)=g0_ttheta_ads_sph_rhophi
        gads_ll_sph_xx(1,3,3,3)=g0_ttheta_ads_sph_thetatheta
        gads_ll_sph_xx(1,3,3,4)=g0_ttheta_ads_sph_thetaphi
        gads_ll_sph_xx(1,3,4,4)=g0_ttheta_ads_sph_phiphi

        gads_ll_sph_x(1,4,1)   =g0_tphi_ads_sph_t
        gads_ll_sph_x(1,4,2)   =g0_tphi_ads_sph_rho
        gads_ll_sph_x(1,4,3)   =g0_tphi_ads_sph_theta
        gads_ll_sph_x(1,4,4)   =g0_tphi_ads_sph_phi
        gads_ll_sph_xx(1,4,1,1)=g0_tphi_ads_sph_tt
        gads_ll_sph_xx(1,4,1,2)=g0_tphi_ads_sph_trho
        gads_ll_sph_xx(1,4,1,3)=g0_tphi_ads_sph_ttheta
        gads_ll_sph_xx(1,4,1,4)=g0_tphi_ads_sph_tphi
        gads_ll_sph_xx(1,4,2,2)=g0_tphi_ads_sph_rhorho
        gads_ll_sph_xx(1,4,2,3)=g0_tphi_ads_sph_rhotheta
        gads_ll_sph_xx(1,4,2,4)=g0_tphi_ads_sph_rhophi
        gads_ll_sph_xx(1,4,3,3)=g0_tphi_ads_sph_thetatheta
        gads_ll_sph_xx(1,4,3,4)=g0_tphi_ads_sph_thetaphi
        gads_ll_sph_xx(1,4,4,4)=g0_tphi_ads_sph_phiphi

        gads_ll_sph_x(2,2,1)   =g0_rhorho_ads_sph_t
        gads_ll_sph_x(2,2,2)   =g0_rhorho_ads_sph_rho
        gads_ll_sph_x(2,2,3)   =g0_rhorho_ads_sph_theta
        gads_ll_sph_x(2,2,4)   =g0_rhorho_ads_sph_phi
        gads_ll_sph_xx(2,2,1,1)=g0_rhorho_ads_sph_tt
        gads_ll_sph_xx(2,2,1,2)=g0_rhorho_ads_sph_trho
        gads_ll_sph_xx(2,2,1,3)=g0_rhorho_ads_sph_ttheta
        gads_ll_sph_xx(2,2,1,4)=g0_rhorho_ads_sph_tphi
        gads_ll_sph_xx(2,2,2,2)=g0_rhorho_ads_sph_rhorho
        gads_ll_sph_xx(2,2,2,3)=g0_rhorho_ads_sph_rhotheta
        gads_ll_sph_xx(2,2,2,4)=g0_rhorho_ads_sph_rhophi
        gads_ll_sph_xx(2,2,3,3)=g0_rhorho_ads_sph_thetatheta
        gads_ll_sph_xx(2,2,3,4)=g0_rhorho_ads_sph_thetaphi
        gads_ll_sph_xx(2,2,4,4)=g0_rhorho_ads_sph_phiphi

        gads_ll_sph_x(2,3,1)   =g0_rhotheta_ads_sph_t
        gads_ll_sph_x(2,3,2)   =g0_rhotheta_ads_sph_rho
        gads_ll_sph_x(2,3,3)   =g0_rhotheta_ads_sph_theta
        gads_ll_sph_x(2,3,4)   =g0_rhotheta_ads_sph_phi
        gads_ll_sph_xx(2,3,1,1)=g0_rhotheta_ads_sph_tt
        gads_ll_sph_xx(2,3,1,2)=g0_rhotheta_ads_sph_trho
        gads_ll_sph_xx(2,3,1,3)=g0_rhotheta_ads_sph_ttheta
        gads_ll_sph_xx(2,3,1,4)=g0_rhotheta_ads_sph_tphi
        gads_ll_sph_xx(2,3,2,2)=g0_rhotheta_ads_sph_rhorho
        gads_ll_sph_xx(2,3,2,3)=g0_rhotheta_ads_sph_rhotheta
        gads_ll_sph_xx(2,3,2,4)=g0_rhotheta_ads_sph_rhophi
        gads_ll_sph_xx(2,3,3,3)=g0_rhotheta_ads_sph_thetatheta
        gads_ll_sph_xx(2,3,3,4)=g0_rhotheta_ads_sph_thetaphi
        gads_ll_sph_xx(2,3,4,4)=g0_rhotheta_ads_sph_phiphi

        gads_ll_sph_x(2,4,1)   =g0_rhophi_ads_sph_t
        gads_ll_sph_x(2,4,2)   =g0_rhophi_ads_sph_rho
        gads_ll_sph_x(2,4,3)   =g0_rhophi_ads_sph_theta
        gads_ll_sph_x(2,4,4)   =g0_rhophi_ads_sph_phi
        gads_ll_sph_xx(2,4,1,1)=g0_rhophi_ads_sph_tt
        gads_ll_sph_xx(2,4,1,2)=g0_rhophi_ads_sph_trho
        gads_ll_sph_xx(2,4,1,3)=g0_rhophi_ads_sph_ttheta
        gads_ll_sph_xx(2,4,1,4)=g0_rhophi_ads_sph_tphi
        gads_ll_sph_xx(2,4,2,2)=g0_rhophi_ads_sph_rhorho
        gads_ll_sph_xx(2,4,2,3)=g0_rhophi_ads_sph_rhotheta
        gads_ll_sph_xx(2,4,2,4)=g0_rhophi_ads_sph_rhophi
        gads_ll_sph_xx(2,4,3,3)=g0_rhophi_ads_sph_thetatheta
        gads_ll_sph_xx(2,4,3,4)=g0_rhophi_ads_sph_thetaphi
        gads_ll_sph_xx(2,4,4,4)=g0_rhophi_ads_sph_phiphi

        gads_ll_sph_x(3,3,1)   =g0_thetatheta_ads_sph_t
        gads_ll_sph_x(3,3,2)   =g0_thetatheta_ads_sph_rho
        gads_ll_sph_x(3,3,3)   =g0_thetatheta_ads_sph_theta
        gads_ll_sph_x(3,3,4)   =g0_thetatheta_ads_sph_phi
        gads_ll_sph_xx(3,3,1,1)=g0_thetatheta_ads_sph_tt
        gads_ll_sph_xx(3,3,1,2)=g0_thetatheta_ads_sph_trho
        gads_ll_sph_xx(3,3,1,3)=g0_thetatheta_ads_sph_ttheta
        gads_ll_sph_xx(3,3,1,4)=g0_thetatheta_ads_sph_tphi
        gads_ll_sph_xx(3,3,2,2)=g0_thetatheta_ads_sph_rhorho
        gads_ll_sph_xx(3,3,2,3)=g0_thetatheta_ads_sph_rhotheta
        gads_ll_sph_xx(3,3,2,4)=g0_thetatheta_ads_sph_rhophi
        gads_ll_sph_xx(3,3,3,3)=g0_thetatheta_ads_sph_thetatheta
        gads_ll_sph_xx(3,3,3,4)=g0_thetatheta_ads_sph_thetaphi
        gads_ll_sph_xx(3,3,4,4)=g0_thetatheta_ads_sph_phiphi

        gads_ll_sph_x(3,4,1)   =g0_thetaphi_ads_sph_t
        gads_ll_sph_x(3,4,2)   =g0_thetaphi_ads_sph_rho
        gads_ll_sph_x(3,4,3)   =g0_thetaphi_ads_sph_theta
        gads_ll_sph_x(3,4,4)   =g0_thetaphi_ads_sph_phi
        gads_ll_sph_xx(3,4,1,1)=g0_thetaphi_ads_sph_tt
        gads_ll_sph_xx(3,4,1,2)=g0_thetaphi_ads_sph_trho
        gads_ll_sph_xx(3,4,1,3)=g0_thetaphi_ads_sph_ttheta
        gads_ll_sph_xx(3,4,1,4)=g0_thetaphi_ads_sph_tphi
        gads_ll_sph_xx(3,4,2,2)=g0_thetaphi_ads_sph_rhorho
        gads_ll_sph_xx(3,4,2,3)=g0_thetaphi_ads_sph_rhotheta
        gads_ll_sph_xx(3,4,2,4)=g0_thetaphi_ads_sph_rhophi
        gads_ll_sph_xx(3,4,3,3)=g0_thetaphi_ads_sph_thetatheta
        gads_ll_sph_xx(3,4,3,4)=g0_thetaphi_ads_sph_thetaphi
        gads_ll_sph_xx(3,4,4,4)=g0_thetaphi_ads_sph_phiphi

        gads_ll_sph_x(4,4,1)   =g0_phiphi_ads_sph_t
        gads_ll_sph_x(4,4,2)   =g0_phiphi_ads_sph_rho
        gads_ll_sph_x(4,4,3)   =g0_phiphi_ads_sph_theta
        gads_ll_sph_x(4,4,4)   =g0_phiphi_ads_sph_phi
        gads_ll_sph_xx(4,4,1,1)=g0_phiphi_ads_sph_tt
        gads_ll_sph_xx(4,4,1,2)=g0_phiphi_ads_sph_trho
        gads_ll_sph_xx(4,4,1,3)=g0_phiphi_ads_sph_ttheta
        gads_ll_sph_xx(4,4,1,4)=g0_phiphi_ads_sph_tphi
        gads_ll_sph_xx(4,4,2,2)=g0_phiphi_ads_sph_rhorho
        gads_ll_sph_xx(4,4,2,3)=g0_phiphi_ads_sph_rhotheta
        gads_ll_sph_xx(4,4,2,4)=g0_phiphi_ads_sph_rhophi
        gads_ll_sph_xx(4,4,3,3)=g0_phiphi_ads_sph_thetatheta
        gads_ll_sph_xx(4,4,3,4)=g0_phiphi_ads_sph_thetaphi
        gads_ll_sph_xx(4,4,4,4)=g0_phiphi_ads_sph_phiphi

       end if


        do a=1,3
          do b=a+1,4
            gads_ll_sph(b,a)=gads_ll_sph(a,b)
            gads_uu_sph(b,a)=gads_uu_sph(a,b)
            if (calc_der) then
             do c=1,4
               gads_ll_sph_x(b,a,c)=gads_ll_sph_x(a,b,c)
               do d=1,4
                 gads_ll_sph_xx(b,a,c,d)=gads_ll_sph_xx(a,b,c,d)
                 gads_ll_sph_xx(c,d,b,a)=gads_ll_sph_xx(c,d,a,b)
               end do
             end do
            end if
          end do
        end do

        !define transformation matrix between spherical to Cartesian coordinates, 
        !e.g. dxsph_dxcar(2,3)=drho/dtheta, d2xsph_dxcardxcar(2,3,4)=d^2rho/(dtheta dphi), 
        ! d3xsph_dxcardxcardxcar(2,3,4,2)=d^2rho/(dtheta dphi drho)

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

       if (calc_der) then
        d2xsph_dxcardxcar(1,1,1)=
     -   0
        d2xsph_dxcardxcar(1,1,2)=
     -   0
        d2xsph_dxcardxcar(1,1,3)=
     -   0
        d2xsph_dxcardxcar(1,1,4)=
     -   0
        d2xsph_dxcardxcar(1,2,2)=
     -   0
        d2xsph_dxcardxcar(1,2,3)=
     -   0
        d2xsph_dxcardxcar(1,2,4)=
     -   0
        d2xsph_dxcardxcar(1,3,3)=
     -   0
        d2xsph_dxcardxcar(1,3,4)=
     -   0
        d2xsph_dxcardxcar(1,4,4)=
     -   0

        d2xsph_dxcardxcar(2,1,1)=
     -   0
        d2xsph_dxcardxcar(2,1,2)=
     -   0
        d2xsph_dxcardxcar(2,1,3)=
     -   0
        d2xsph_dxcardxcar(2,1,4)=
     -   0
        d2xsph_dxcardxcar(2,2,2)=
     -   (rho0**2 - x0**2)/rho0**3
        d2xsph_dxcardxcar(2,2,3)=
     -   -((x0*y0)/rho0**3)
        d2xsph_dxcardxcar(2,2,4)=
     -   -((x0*z0)/rho0**3)
        d2xsph_dxcardxcar(2,3,3)=
     -   (rho0**2 - y0**2)/rho0**3
        d2xsph_dxcardxcar(2,3,4)=
     -   -((y0*z0)/rho0**3)
        d2xsph_dxcardxcar(2,4,4)=
     -   (rho0**2 - z0**2)/rho0**3

        d2xsph_dxcardxcar(3,1,1)=
     -   0
        d2xsph_dxcardxcar(3,1,2)=
     -   0
        d2xsph_dxcardxcar(3,1,3)=
     -   0
        d2xsph_dxcardxcar(3,1,4)=
     -   0
        d2xsph_dxcardxcar(3,2,2)=
     -   (2*x0*Sqrt(rho0**2 - x0**2))/rho0**4
        d2xsph_dxcardxcar(3,2,3)=
     -   ((rho0**2 - 2*x0**2)*y0)/
     -  (rho0**4*Sqrt(rho0**2 - x0**2))
        d2xsph_dxcardxcar(3,2,4)=
     -   ((rho0**2 - 2*x0**2)*z0)/
     -  (rho0**4*Sqrt(rho0**2 - x0**2))
        d2xsph_dxcardxcar(3,3,3)=
     -   (x0*(rho0**4 + 2*x0**2*y0**2 - 
     -      rho0**2*(x0**2 + 3*y0**2)))/
     -  (rho0**4*(rho0**2 - x0**2)**1.5)
        d2xsph_dxcardxcar(3,3,4)=
     -   ((-3*rho0**2*x0 + 2*x0**3)*y0*z0)/
     -  (rho0**4*(rho0**2 - x0**2)**1.5)
        d2xsph_dxcardxcar(3,4,4)=
     -   (x0*(rho0**4 + 2*x0**2*z0**2 - 
     -      rho0**2*(x0**2 + 3*z0**2)))/
     -  (rho0**4*(rho0**2 - x0**2)**1.5)

        d2xsph_dxcardxcar(4,1,1)=
     -   0
        d2xsph_dxcardxcar(4,1,2)=
     -   0
        d2xsph_dxcardxcar(4,1,3)=
     -   0
        d2xsph_dxcardxcar(4,1,4)=
     -   0
        d2xsph_dxcardxcar(4,2,2)=
     -   0
        d2xsph_dxcardxcar(4,2,3)=
     -   0
        d2xsph_dxcardxcar(4,2,4)=
     -   0
        d2xsph_dxcardxcar(4,3,3)=
     -   (2*y0*z0)/(y0**2 + z0**2)**2
        d2xsph_dxcardxcar(4,3,4)=
     -   (-y0**2 + z0**2)/(y0**2 + z0**2)**2
        d2xsph_dxcardxcar(4,4,4)=
     -   (-2*y0*z0)/(y0**2 + z0**2)**2

        d3xsph_dxcardxcardxcar(1,1,1,1)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,1,2)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,1,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,1,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,3,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,3,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,3,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,4,4,4)=
     -   0



        d3xsph_dxcardxcardxcar(2,1,1,1)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,1,2)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,1,3)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,1,4)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(2,2,2,2)=
     -   (3*x0*(-rho0 + x0)*(rho0 + x0))/rho0**5
        d3xsph_dxcardxcardxcar(2,2,2,3)=
     -   -(((rho0**2 - 3*x0**2)*y0)/rho0**5)
        d3xsph_dxcardxcardxcar(2,2,2,4)=
     -   -(((rho0**2 - 3*x0**2)*z0)/rho0**5)
        d3xsph_dxcardxcardxcar(2,2,3,3)=
     -   -((x0*(rho0**2 - 3*y0**2))/rho0**5)
        d3xsph_dxcardxcardxcar(2,2,3,4)=
     -   (3*x0*y0*z0)/rho0**5
        d3xsph_dxcardxcardxcar(2,2,4,4)=
     -   -((x0*(rho0**2 - 3*z0**2))/rho0**5)
        d3xsph_dxcardxcardxcar(2,3,3,3)=
     -   (3*y0*(-rho0 + y0)*(rho0 + y0))/rho0**5
        d3xsph_dxcardxcardxcar(2,3,3,4)=
     -   -(((rho0**2 - 3*y0**2)*z0)/rho0**5)
        d3xsph_dxcardxcardxcar(2,3,4,4)=
     -   -((y0*(rho0**2 - 3*z0**2))/rho0**5)
        d3xsph_dxcardxcardxcar(2,4,4,4)=
     -   (3*z0*(-rho0 + z0)*(rho0 + z0))/rho0**5

        d3xsph_dxcardxcardxcar(3,1,1,1)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,1,2)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,1,3)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,1,4)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(3,2,2,2)=
     -   (2*(rho0**2 - 4*x0**2)*
     -   Sqrt(rho0**2 - x0**2))/rho0**6
        d3xsph_dxcardxcardxcar(3,2,2,3)=
     -   (-6*rho0**2*x0*y0 + 8*x0**3*y0)/
     -  (rho0**6*Sqrt(rho0**2 - x0**2))
        d3xsph_dxcardxcardxcar(3,2,2,4)=
     -   (-6*rho0**2*x0*z0 + 8*x0**3*z0)/
     -  (rho0**6*Sqrt(rho0**2 - x0**2))
        d3xsph_dxcardxcardxcar(3,2,3,3)=
     -   (rho0**6 - 8*x0**4*y0**2 - 
     -    3*rho0**4*(x0**2 + y0**2) + 
     -    2*rho0**2*x0**2*(x0**2 + 6*y0**2))/
     -  (rho0**6*(rho0**2 - x0**2)**1.5)
        d3xsph_dxcardxcardxcar(3,2,3,4)=
     -   ((-3*rho0**4 + 12*rho0**2*x0**2 - 8*x0**4)*y0*
     -    z0)/(rho0**6*(rho0**2 - x0**2)**1.5)
        d3xsph_dxcardxcardxcar(3,2,4,4)=
     -   (rho0**6 - 8*x0**4*z0**2 - 
     -    3*rho0**4*(x0**2 + z0**2) + 
     -    2*rho0**2*x0**2*(x0**2 + 6*z0**2))/
     -  (rho0**6*(rho0**2 - x0**2)**1.5)
        d3xsph_dxcardxcardxcar(3,3,3,3)=
     -   (x0*y0*(-9*rho0**6 + 8*x0**4*y0**2 + 
     -      15*rho0**4*(x0**2 + y0**2) - 
     -      2*rho0**2*x0**2*(3*x0**2 + 10*y0**2)))/
     -  (rho0**6*(rho0**2 - x0**2)**2.5)
        d3xsph_dxcardxcardxcar(3,3,3,4)=
     -   (x0*(-3*rho0**6 + 8*x0**4*y0**2 + 
     -      5*rho0**4*(x0**2 + 3*y0**2) - 
     -      2*rho0**2*x0**2*(x0**2 + 10*y0**2))*z0)/
     -  (rho0**6*(rho0**2 - x0**2)**2.5)
        d3xsph_dxcardxcardxcar(3,3,4,4)=
     -   (x0*y0*(-3*rho0**6 + 8*x0**4*z0**2 + 
     -      5*rho0**4*(x0**2 + 3*z0**2) - 
     -      2*rho0**2*x0**2*(x0**2 + 10*z0**2)))/
     -  (rho0**6*(rho0**2 - x0**2)**2.5)
        d3xsph_dxcardxcardxcar(3,4,4,4)=
     -   (x0*z0*(-9*rho0**6 + 8*x0**4*z0**2 + 
     -      15*rho0**4*(x0**2 + z0**2) - 
     -      2*rho0**2*x0**2*(3*x0**2 + 10*z0**2)))/
     -  (rho0**6*(rho0**2 - x0**2)**2.5)

        d3xsph_dxcardxcardxcar(4,1,1,1)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,1,2)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,1,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,1,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,3,3,3)=
     -   (2*z0*(-3*y0**2 + z0**2))/(y0**2 + z0**2)**3
        d3xsph_dxcardxcardxcar(4,3,3,4)=
     -   (2*y0*(y0**2 - 3*z0**2))/(y0**2 + z0**2)**3
        d3xsph_dxcardxcardxcar(4,3,4,4)=
     -   (-2*(-3*y0**2*z0 + z0**3))/(y0**2 + z0**2)**3
        d3xsph_dxcardxcardxcar(4,4,4,4)=
     -   (-2*y0*(y0**2 - 3*z0**2))/(y0**2 + z0**2)**3


        do a=1,4
         do b=1,3
          do c=b+1,4
           d2xsph_dxcardxcar(a,c,b)=d2xsph_dxcardxcar(a,b,c)
          end do
         end do
        end do


        do a=1,4
         do b=1,4
          do c=1,4
           do d=1,4
            if ((max(b,c,d).eq.d).and.(max(b,c).eq.c)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,b,c,d)
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,b,c,d)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,b,c,d)
            else if ((max(b,c,d).eq.d).and.(max(b,c).eq.b)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,c,b,d)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,c,b,d)

            else if ((max(b,c,d).eq.c).and.(max(b,d).eq.d)) then
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,b,d,c)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,b,d,c)
            else if ((max(b,c,d).eq.c).and.(max(b,d).eq.b)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,d,b,c)
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,d,b,c)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,d,b,c)

            else if ((max(b,c,d).eq.b).and.(max(c,d).eq.d)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,c,d,b)
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,c,d,b)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,c,d,b)
            else if ((max(b,c,d).eq.b).and.(max(c,d).eq.c)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,d,c,b)
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,d,c,b)
            end if
           end do
          end do
         end do
        end do
       end if !closes condition on calc_der

        !compute Cartesian quantities in terms of spherical ones
        do a=1,4
         do b=1,4
          gads_ll(a,b)=0
          do c=1,4
            gads_ll_x(a,b,c)=0
           do d=1,4
            gads_ll_xx(a,b,c,d)=0
            gads_ll(a,b)=gads_ll(a,b)+
     &         dxsph_dxcar(c,a)*dxsph_dxcar(d,b)*gads_ll_sph(c,d)
            if (calc_der) then
              do e=1,4
               gads_ll_x(a,b,c)=gads_ll_x(a,b,c)+
     &       d2xsph_dxcardxcar(d,c,a)*dxsph_dxcar(e,b)
     &       *gads_ll_sph(d,e)+
     &       dxsph_dxcar(d,a)*d2xsph_dxcardxcar(e,c,b)
     &       *gads_ll_sph(d,e)
               do f=1,4
                gads_ll_x(a,b,c)=gads_ll_x(a,b,c)+
     &       dxsph_dxcar(d,a)*dxsph_dxcar(e,b)*dxsph_dxcar(f,c)
     &       *gads_ll_sph_x(d,e,f)
                gads_ll_xx(a,b,c,d)=gads_ll_xx(a,b,c,d)+
     &       d3xsph_dxcardxcardxcar(e,d,c,a)*dxsph_dxcar(f,b)
     &       *gads_ll_sph(e,f)+
     &       d2xsph_dxcardxcar(e,c,a)*d2xsph_dxcardxcar(f,d,b)
     &       *gads_ll_sph(e,f)+
     &       d2xsph_dxcardxcar(e,d,a)*d2xsph_dxcardxcar(f,c,b)
     &       *gads_ll_sph(e,f)+
     &       dxsph_dxcar(e,a)*d3xsph_dxcardxcardxcar(f,d,c,b)
     &       *gads_ll_sph(e,f)
                do g=1,4
             gads_ll_xx(a,b,c,d)=gads_ll_xx(a,b,c,d)+
     &       d2xsph_dxcardxcar(e,c,a)*dxsph_dxcar(f,b)*dxsph_dxcar(g,d)
     &       *gads_ll_sph_x(e,f,g)+
     &       dxsph_dxcar(e,a)*d2xsph_dxcardxcar(f,c,b)*dxsph_dxcar(g,d)
     &       *gads_ll_sph_x(e,f,g)+
     &       d2xsph_dxcardxcar(e,d,a)*dxsph_dxcar(f,b)*dxsph_dxcar(g,c)
     &       *gads_ll_sph_x(e,f,g)+
     &       dxsph_dxcar(e,a)*d2xsph_dxcardxcar(f,d,b)*dxsph_dxcar(g,c)
     &       *gads_ll_sph_x(e,f,g)+
     &       dxsph_dxcar(e,a)*dxsph_dxcar(f,b)*d2xsph_dxcardxcar(g,d,c)
     &       *gads_ll_sph_x(e,f,g)
                 do h=1,4
                  gads_ll_xx(a,b,c,d)=gads_ll_xx(a,b,c,d)+
     &       dxsph_dxcar(e,a)*dxsph_dxcar(f,b)*
     &       dxsph_dxcar(g,c)*dxsph_dxcar(h,d)*gads_ll_sph_xx(e,f,g,h)
                 end do
                end do
               end do
              end do
            end if
           end do
          end do
         end do
        end do

       !some of the dxsph_dxcar diverge at y=z=0 so we need to consider this case separately
       ! we only need to consider the quantities with a y or z index

       if ((abs(y0).lt.10.0d0**(-10)).and.
     &     (abs(z0).lt.10.0d0**(-10))) then
        gads_ll(1,2)=
     -   0
        gads_ll(1,3)=
     -   0
        gads_ll(1,4)=
     -   0
        gads_ll(2,2)=
     -   (4*L**2*(1 + x0**2)**2)/
     -  ((-1 + x0**2)**2*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2))
        gads_ll(2,3)=
     -   0
        gads_ll(2,4)=
     -   0
        gads_ll(3,3)=
     -   4/(-1 + x0**2)**2
        gads_ll(3,4)=
     -   0
        gads_ll(4,4)=
     -   4/(-1 + x0**2)**2

       if (calc_der) then
        gads_ll_x(1,1,2)=
     -   (8*(x0 + x0**3))/(L**2*(-1 + x0**2)**3)
        gads_ll_x(1,1,3)=
     -   0
        gads_ll_x(1,1,4)=
     -   0
        gads_ll_xx(1,1,1,2)=
     -   0
        gads_ll_xx(1,1,1,3)=
     -   0
        gads_ll_xx(1,1,1,4)=
     -   0
        gads_ll_xx(1,1,2,2)=
     -   (-8*(1 + 8*x0**2 + 3*x0**4))/
     -  (L**2*(-1 + x0**2)**4)
        gads_ll_xx(1,1,2,3)=
     -   0
        gads_ll_xx(1,1,2,4)=
     -   0
        gads_ll_xx(1,1,3,3)=
     -   (8*(1 + x0**2))/(L**2*(-1 + x0**2)**3)
        gads_ll_xx(1,1,3,4)=
     -   0
        gads_ll_xx(1,1,4,4)=
     -   (8*(1 + x0**2))/(L**2*(-1 + x0**2)**3)

        gads_ll_x(1,2,1)=
     -   0
        gads_ll_x(1,2,2)=
     -   0
        gads_ll_x(1,2,3)=
     -   0
        gads_ll_x(1,2,4)=
     -   0
        gads_ll_xx(1,2,1,1)=
     -   0
        gads_ll_xx(1,2,1,2)=
     -   0
        gads_ll_xx(1,2,1,3)=
     -   0
        gads_ll_xx(1,2,1,4)=
     -   0
        gads_ll_xx(1,2,2,2)=
     -   0
        gads_ll_xx(1,2,2,3)=
     -   0
        gads_ll_xx(1,2,2,4)=
     -   0
        gads_ll_xx(1,2,3,3)=
     -   0
        gads_ll_xx(1,2,3,4)=
     -   0
        gads_ll_xx(1,2,4,4)=
     -   0

        gads_ll_x(1,3,1)=
     -   0
        gads_ll_x(1,3,2)=
     -   0
        gads_ll_x(1,3,3)=
     -   0
        gads_ll_x(1,3,4)=
     -   0
        gads_ll_xx(1,3,1,1)=
     -   0
        gads_ll_xx(1,3,1,2)=
     -   0
        gads_ll_xx(1,3,1,3)=
     -   0
        gads_ll_xx(1,3,1,4)=
     -   0
        gads_ll_xx(1,3,2,2)=
     -   0
        gads_ll_xx(1,3,2,3)=
     -   0
        gads_ll_xx(1,3,2,4)=
     -   0
        gads_ll_xx(1,3,3,3)=
     -   0
        gads_ll_xx(1,3,3,4)=
     -   0
        gads_ll_xx(1,3,4,4)=
     -   0

        gads_ll_x(1,4,1)=
     -   0
        gads_ll_x(1,4,2)=
     -   0
        gads_ll_x(1,4,3)=
     -   0
        gads_ll_x(1,4,4)=
     -   0
        gads_ll_xx(1,4,1,1)=
     -   0
        gads_ll_xx(1,4,1,2)=
     -   0
        gads_ll_xx(1,4,1,3)=
     -   0
        gads_ll_xx(1,4,1,4)=
     -   0
        gads_ll_xx(1,4,2,2)=
     -   0
        gads_ll_xx(1,4,2,3)=
     -   0
        gads_ll_xx(1,4,2,4)=
     -   0
        gads_ll_xx(1,4,3,3)=
     -   0
        gads_ll_xx(1,4,3,4)=
     -   0
        gads_ll_xx(1,4,4,4)=
     -   0

        gads_ll_x(2,2,1)=
     -   0
        gads_ll_x(2,2,2)=
     -    (-16*L**2*x0*(1 + x0**2)*
     -    (L**2*(-1 + x0**2)**2*(3 + x0**2) + 
     -      2*(-1 + 4*x0**2 + x0**4)))/
     -  ((-1 + x0**2)**3*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2)**2)
        gads_ll_x(2,2,3)=
     -   0
        gads_ll_x(2,2,4)=
     -   0
        gads_ll_xx(2,2,1,1)=
     -   0
        gads_ll_xx(2,2,1,2)=
     -   0
        gads_ll_xx(2,2,1,3)=
     -   0
        gads_ll_xx(2,2,1,4)=
     -   0
        gads_ll_xx(2,2,2,2)=
     -   (16*(L**6*(-1 + x0**2)**4*
     -       (3 + 39*x0**2 + 33*x0**4 + 5*x0**6)
     -       + 8*L**2*x0**2*
     -       (3 - 12*x0**2 + 26*x0**4 + 
     -         28*x0**6 + 3*x0**8) + 
     -      2*L**4*(-1 + x0**2)**2*
     -       (-1 - 22*x0**2 + 80*x0**4 + 
     -         78*x0**6 + 9*x0**8)))/
     -  ((-1 + x0**2)**4*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2)**3)
        gads_ll_xx(2,2,2,3)=
     -   0
        gads_ll_xx(2,2,2,4)=
     -   0
        gads_ll_xx(2,2,3,3)=
     -   (-16*(-8*x0**2*(-1 + x0**2) + 
     -      8*L**2*x0**2*(-1 + 3*x0**2) + 
     -      L**4*(-1 + x0**2)**2*(1 + 6*x0**2 + x0**4)))
     -   /((-1 + x0**2)**3*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2)**2)
        gads_ll_xx(2,2,3,4)=
     -   0
        gads_ll_xx(2,2,4,4)=
     -   (-16*(-8*x0**2*(-1 + x0**2) + 
     -      8*L**2*x0**2*(-1 + 3*x0**2) + 
     -      L**4*(-1 + x0**2)**2*(1 + 6*x0**2 + x0**4)))
     -   /((-1 + x0**2)**3*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2)**2)

        gads_ll_x(2,3,1)=
     -   0
        gads_ll_x(2,3,2)=
     -   0
        gads_ll_x(2,3,3)=
     -   (16*(-1 + L**2)*x0)/
     -  ((-1 + x0**2)**2*
     -    (L**2 + 4*x0**2 - 2*L**2*x0**2 + L**2*x0**4))
        gads_ll_x(2,3,4)=
     -   0
        gads_ll_xx(2,3,1,1)=
     -   0
        gads_ll_xx(2,3,1,2)=
     -   0
        gads_ll_xx(2,3,1,3)=
     -   0
        gads_ll_xx(2,3,1,4)=
     -   0
        gads_ll_xx(2,3,2,2)=
     -   0
        gads_ll_xx(2,3,2,3)=
     -   (-16*(-1 + L**2)*
     -    (L**2 - 4*x0**2 + 5*L**2*x0**2 + 20*x0**4 - 
     -      13*L**2*x0**4 + 7*L**2*x0**6))/
     -  ((-1 + x0**2)**3*
     -    (L**2 + 4*x0**2 - 2*L**2*x0**2 + L**2*x0**4)**
     -     2)
        gads_ll_xx(2,3,2,4)=
     -   0
        gads_ll_xx(2,3,3,3)=
     -   0
        gads_ll_xx(2,3,3,4)=
     -   0
        gads_ll_xx(2,3,4,4)=
     -   0

        gads_ll_x(2,4,1)=
     -   0
        gads_ll_x(2,4,2)=
     -   0
        gads_ll_x(2,4,3)=
     -   0
        gads_ll_x(2,4,4)=
     -   (16*(-1 + L**2)*x0)/
     -  ((-1 + x0**2)**2*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2))
        gads_ll_xx(2,4,1,1)=
     -   0
        gads_ll_xx(2,4,1,2)=
     -   0
        gads_ll_xx(2,4,1,3)=
     -   0
        gads_ll_xx(2,4,1,4)=
     -   0
        gads_ll_xx(2,4,2,2)=
     -   0
        gads_ll_xx(2,4,2,3)=
     -   0
        gads_ll_xx(2,4,2,4)=
     -   (-16*(-1 + L**2)*
     -    (4*x0**2*(-1 + 5*x0**2) + 
     -      L**2*(-1 + x0**2)**2*(1 + 7*x0**2)))/
     -  ((-1 + x0**2)**3*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2)**2)
        gads_ll_xx(2,4,3,3)=
     -   0
        gads_ll_xx(2,4,3,4)=
     -   0
        gads_ll_xx(2,4,4,4)=
     -   0

        gads_ll_x(3,3,1)=
     -   0
        gads_ll_x(3,3,2)=
     -   (-16*x0)/(-1 + x0**2)**3
        gads_ll_x(3,3,3)=
     -   0
        gads_ll_x(3,3,4)=
     -   0
        gads_ll_xx(3,3,1,1)=
     -   0
        gads_ll_xx(3,3,1,2)=
     -   0
        gads_ll_xx(3,3,1,3)=
     -   0
        gads_ll_xx(3,3,1,4)=
     -   0
        gads_ll_xx(3,3,2,2)=
     -   (16*(1 + 5*x0**2))/(-1 + x0**2)**4
        gads_ll_xx(3,3,2,3)=
     -   0
        gads_ll_xx(3,3,2,4)=
     -   0
        gads_ll_xx(3,3,3,3)=
     -   (-16*(-2 + 6*x0**2 + 
     -      L**2*(3 - 4*x0**2 + x0**4)))/
     -  ((-1 + x0**2)**3*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2))
        gads_ll_xx(3,3,3,4)=
     -   0
        gads_ll_xx(3,3,4,4)=
     -   -16/(-1 + x0**2)**3

        gads_ll_x(3,4,1)=
     -   0
        gads_ll_x(3,4,2)=
     -   0
        gads_ll_x(3,4,3)=
     -   0
        gads_ll_x(3,4,4)=
     -   0
        gads_ll_xx(3,4,1,1)=
     -   0
        gads_ll_xx(3,4,1,2)=
     -   0
        gads_ll_xx(3,4,1,3)=
     -   0
        gads_ll_xx(3,4,1,4)=
     -   0
        gads_ll_xx(3,4,2,2)=
     -   0
        gads_ll_xx(3,4,2,3)=
     -   0
        gads_ll_xx(3,4,2,4)=
     -   0
        gads_ll_xx(3,4,3,3)=
     -   0
        gads_ll_xx(3,4,3,4)=
     -   (16*(-1 + L**2))/
     -  ((-1 + x0**2)**2*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2))
        gads_ll_xx(3,4,4,4)=
     -   0

        gads_ll_x(4,4,1)=
     -   0
        gads_ll_x(4,4,2)=
     -   (-16*x0)/(-1 + x0**2)**3
        gads_ll_x(4,4,3)=
     -   0
        gads_ll_x(4,4,4)=
     -   0
        gads_ll_xx(4,4,1,1)=
     -   0
        gads_ll_xx(4,4,1,2)=
     -   0
        gads_ll_xx(4,4,1,3)=
     -   0
        gads_ll_xx(4,4,1,4)=
     -   0
        gads_ll_xx(4,4,2,2)=
     -   (16*(1 + 5*x0**2))/(-1 + x0**2)**4
        gads_ll_xx(4,4,2,3)=
     -   0
        gads_ll_xx(4,4,2,4)=
     -   0
        gads_ll_xx(4,4,3,3)=
     -   -16/(-1 + x0**2)**3
        gads_ll_xx(4,4,3,4)=
     -   0
        gads_ll_xx(4,4,4,4)=
     -    (-16*(-2 + 6*x0**2 + 
     -      L**2*(3 - 4*x0**2 + x0**4)))/
     -  ((-1 + x0**2)**3*
     -    (4*x0**2 + L**2*(-1 + x0**2)**2))
       end if

        do a=1,3
          do b=a+1,4
            gads_ll(b,a)=gads_ll(a,b)
            if (calc_der) then
             do c=1,4
               gads_ll_x(b,a,c)=gads_ll_x(a,b,c)
             end do
            end if
          end do
        end do

        if (calc_der) then
          do a=1,4
            do b=1,4
              do c=1,4
                do d=1,4
                  gads_ll_xx(a,b,c,d)=
     &               gads_ll_xx(min(a,b),max(a,b),min(c,d),max(c,d))
                end do
              end do
            end do
          end do
        end if

       end if !closes condition on y=z=0

        call calc_g0uu(gads_ll(1,1),gads_ll(1,2),
     &         gads_ll(1,3),gads_ll(1,4),
     &         gads_ll(2,2),gads_ll(2,3),
     &         gads_ll(2,4),
     &         gads_ll(3,3),gads_ll(3,4),
     &         gads_ll(4,4),
     &         gads_uu(1,1),gads_uu(1,2),
     &         gads_uu(1,3),gads_uu(1,4),
     &         gads_uu(2,2),gads_uu(2,3),
     &         gads_uu(2,4),
     &         gads_uu(3,3),gads_uu(3,4),
     &         gads_uu(4,4),detg0_ads0)

        do a=1,3
          do b=a+1,4
            gads_uu(b,a)=gads_uu(a,b) 
          end do
        end do

     	!set values for the FULL bulk scalar field value of the analytic solution

     	phi1ads=0
        if (calc_der) then
         phi1ads_x(1)=0
         phi1ads_x(2)=0
         phi1ads_x(3)=0
         phi1ads_x(4)=0
         phi1ads_xx(1,1)=0
         phi1ads_xx(1,2)=0
         phi1ads_xx(1,3)=0
         phi1ads_xx(1,4)=0
         phi1ads_xx(2,2)=0
         phi1ads_xx(2,3)=0
         phi1ads_xx(2,4)=0
         phi1ads_xx(3,3)=0
         phi1ads_xx(3,4)=0

         do a=1,3
           do b=a+1,4
             phi1ads_xx(b,a)=phi1ads_xx(a,b)
           end do
         end do
        end if

        if (calc_der) then
         do a=1,4
           do b=1,4
             do c=1,4
               gads_uu_x(a,b,c)=
     &               -gads_ll_x(1,1,c)*gads_uu(a,1)*gads_uu(b,1)
     &               -gads_ll_x(1,2,c)*(gads_uu(a,1)*gads_uu(b,2)
     &                                +gads_uu(a,2)*gads_uu(b,1))
     &               -gads_ll_x(1,3,c)*(gads_uu(a,1)*gads_uu(b,3)
     &                                +gads_uu(a,3)*gads_uu(b,1))
     &               -gads_ll_x(1,4,c)*(gads_uu(a,1)*gads_uu(b,4)
     &                                +gads_uu(a,4)*gads_uu(b,1))
     &               -gads_ll_x(2,2,c)*gads_uu(a,2)*gads_uu(b,2)
     &               -gads_ll_x(2,3,c)*(gads_uu(a,2)*gads_uu(b,3)
     &                                +gads_uu(a,3)*gads_uu(b,2))
     &               -gads_ll_x(2,4,c)*(gads_uu(a,2)*gads_uu(b,4)
     &                                +gads_uu(a,4)*gads_uu(b,2))
     &               -gads_ll_x(3,3,c)*gads_uu(a,3)*gads_uu(b,3)
     &               -gads_ll_x(3,4,c)*(gads_uu(a,3)*gads_uu(b,4)
     &                                +gads_uu(a,4)*gads_uu(b,3))
     &               -gads_ll_x(4,4,c)*gads_uu(a,4)*gads_uu(b,4)
             end do
           end do
         end do

          ! give values to the Christoffel symbols
         do a=1,4
           do b=1,4
             do c=1,4
               gammaads_ull(a,b,c)=0
               do d=1,4
                 gammaads_ull(a,b,c)=gammaads_ull(a,b,c)
     &                           +0.5d0*gads_uu(a,d)
     &                                 *(gads_ll_x(c,d,b)
     &                                  -gads_ll_x(b,c,d)
     &                                  +gads_ll_x(d,b,c))
               end do
             end do
           end do
         end do

                ! calculate boxx^c at point i,j
               ! (boxadsx^c = -gads^ab gammaads^c_ab)
               do c=1,4
                 boxadsx_u(c)=-( gammaads_ull(c,1,1)*gads_uu(1,1)+
     &                        gammaads_ull(c,2,2)*gads_uu(2,2)+
     &                        gammaads_ull(c,3,3)*gads_uu(3,3)+
     &                        gammaads_ull(c,4,4)*gads_uu(4,4)+
     &                     2*(gammaads_ull(c,1,2)*gads_uu(1,2)+
     &                        gammaads_ull(c,1,3)*gads_uu(1,3)+
     &                        gammaads_ull(c,1,4)*gads_uu(1,4)+
     &                        gammaads_ull(c,2,3)*gads_uu(2,3)+
     &                        gammaads_ull(c,2,4)*gads_uu(2,4)+
     &                        gammaads_ull(c,3,4)*gads_uu(3,4)) )
               end do

                !compute Hads_l(a) in Cartesian coordinates
               ! (Hads_a = gads_ab boxadsx^b)
               do a=1,4
                 Hads_l(a)=boxadsx_u(1)*gads_ll(a,1)+
     &                     boxadsx_u(2)*gads_ll(a,2)+
     &                     boxadsx_u(3)*gads_ll(a,3)+
     &                     boxadsx_u(4)*gads_ll(a,4)
               end do
        end if

        if (calc_adv_quant) then
        ! calculate Christoffel symbol derivatives at point i,j
        !(gamma^a_bc,e = 1/2 g^ad_,e(g_bd,c  + g_cd,b  - g_bc,d)
        !              +   1/2 g^ad(g_bd,ce + g_cd,be - g_bc,de))
        do a=1,4
          do b=1,4
            do c=1,4
              do e=1,4
                gammaads_ull_x(a,b,c,e)=0
                do d=1,4
                  gammaads_ull_x(a,b,c,e)=gammaads_ull_x(a,b,c,e)
     &              +0.5d0*gads_uu_x(a,d,e)*(gads_ll_x(b,d,c)+
     &                     gads_ll_x(c,d,b)-gads_ll_x(b,c,d))
     &              +0.5d0*gads_uu(a,d)*(gads_ll_xx(b,d,c,e)+
     &                     gads_ll_xx(c,d,b,e)-gads_ll_xx(b,c,d,e))
                end do
              end do
            end do
          end do
        end do

        ! calculate riemann tensor at point i,j
        !(R^a_bcd =gamma^a_bd,c - gamma^a_bc,d
        !          +gamma^a_ce gamma^e_bd - gamma^a_de gamma^e_bc)
        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                riemannads_ulll(a,b,c,d)=
     &                gammaads_ull_x(a,b,d,c)-gammaads_ull_x(a,b,c,d)
                do e=1,4
                   riemannads_ulll(a,b,c,d)=riemannads_ulll(a,b,c,d)
     &               +gammaads_ull(a,c,e)*gammaads_ull(e,b,d)
     &               -gammaads_ull(a,d,e)*gammaads_ull(e,b,c)
                end do
              end do
            end do
          end do
        end do

        ! calculate Ricci tensor at point i,j
        !(R_bd = R^a_bad)
        do b=1,4
          do d=1,4
            ricciads_ll(b,d)=0
            do a=1,4
              ricciads_ll(b,d)=ricciads_ll(b,d)+riemannads_ulll(a,b,a,d)
            end do
          end do
        end do

        ! calculate raised Ricci tensor at point i,j
        !(R_a^b = R_ad g^db)
        do a=1,4
          do b=1,4
            ricciads_lu(a,b)=0
            do d=1,4
              ricciads_lu(a,b)=ricciads_lu(a,b)
     &         +ricciads_ll(a,d)*gads_uu(d,b)
            end do
          end do
        end do

        ! calculate Ricci scalar
        !(R = R_a^a)
        ricciads=0
        do a=1,4
          ricciads=ricciads+ricciads_lu(a,a)
        end do
  
        ! calculates Einstein tensor at point i,j
        !(G_ab = R_ab - 1/2 R g_ab)
        do a=1,4
          do b=1,4
            einsteinads_ll(a,b)=ricciads_ll(a,b)
     &       -0.5d0*ricciads*gads_ll(a,b)
          end do
        end do

        ! calculates stress-energy tensor at point i,j 
        !(T_ab = 2*phi1,a phi1,b - (phi1,c phi1,d) g^cd g_ab + ...)
        grad_phi1ads_sq=0
        do a=1,4
          do b=1,4
            grad_phi1ads_sq=grad_phi1ads_sq
     &                  +phi1ads_x(a)*phi1ads_x(b)*gads_uu(a,b)
          end do
        end do

        do a=1,4
          do b=1,4
            setads_ll(a,b)=
     &            phi1ads_x(a)*phi1ads_x(b)
     &           -gads_ll(a,b)*(grad_phi1ads_sq/2)
          end do
        end do

       end if


!!!!!!!!DEBUG!!!!!
!
!        if ((abs(gads_ll(1,1)-g0_tt_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(1,2)-g0_tx_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(1,3)-g0_ty_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(1,4)-g0_tz_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(2,2)-g0_xx_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(2,3)-g0_xy_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(2,4)-g0_xz_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(3,3)-g0_yy_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(3,4)-g0_yz_ads0).gt.10.0d0**(-10))
!     - .or.(abs(gads_ll(4,4)-g0_zz_ads0).gt.10.0d0**(-10))) then
!
!            write (*,*) "ERROR in AdS metric x0,y0,z0=",x0,y0,z0
!            write (*,*) "rho0,theta0,phi0=",rho0,theta0,phi0
!        do a=1,4
!         do b=1,4
!             write (*,*) 'a,b,gads_ll(a,b)=',
!     -        a,b,gads_ll(a,b)
!         end do
!        end do
!       write (*,*) ' g0_tt_ads0=',g0_tt_ads0
!       write (*,*) ' g0_tx_ads0=',g0_tx_ads0
!       write (*,*) ' g0_ty_ads0=',g0_ty_ads0
!       write (*,*) ' g0_tz_ads0=',g0_tz_ads0
!       write (*,*) ' g0_xx_ads0=',g0_xx_ads0
!       write (*,*) ' g0_xy_ads0=',g0_xy_ads0
!       write (*,*) ' g0_xz_ads0=',g0_xz_ads0
!       write (*,*) ' g0_yy_ads0=',g0_yy_ads0 
!       write (*,*) ' g0_yz_ads0=',g0_yz_ads0
!       write (*,*) ' g0_zz_ads0=',g0_zz_ads0
!        do a=1,4
!         do b=1,4
!             write (*,*) 'a,b,gads_ll_sph(a,b)=',
!     -        a,b,gads_ll_sph(a,b)
!         end do
!        end do
!       write (*,*) ' g0_tt_ads_sph0=',g0_tt_ads_sph0
!       write (*,*) ' g0_trho_ads_sph0=',g0_trho_ads_sph0
!       write (*,*) ' g0_ttheta_ads_sph0=',g0_ttheta_ads_sph0
!       write (*,*) ' g0_tphi_ads_sph0=',g0_tphi_ads_sph0
!       write (*,*) ' g0_rhorho_ads_sph0=',g0_rhorho_ads_sph0
!       write (*,*) ' g0_rhotheta_ads_sph0=',g0_rhotheta_ads_sph0
!       write (*,*) ' g0_rhophi_ads_sph0=',g0_rhophi_ads_sph0
!       write (*,*) ' g0_thetatheta_ads_sph0=',
!     -  g0_thetatheta_ads_sph0 
!       write (*,*) ' g0_thetaphi_ads_sph0=',g0_thetaphi_ads_sph0
!       write (*,*) ' g0_phiphi_ads_sph0=',g0_phiphi_ads_sph0


!        if ((abs(x0-(-1.0+5*dx)).lt.10.0d0**(-10))
!     -  .and.(abs(y0-(-1.0+12*dy)).lt.10.0d0**(-10))
!     -  .and.(abs(z0-(-1.0+10*dz)).lt.10.0d0**(-10))) then
!        write (*,*) "x0,y0,z0=",x0,y0,z0
!        write (*,*) "rho0,theta0,phi0=",rho0,theta0,phi0
!          do c=1,4
!            write (*,*) "c,dxsph_dxcar(1,c)="
!     -       ,c,dxsph_dxcar(1,c)
!          end do
!          do c=1,4
!            write (*,*) "c,dxsph_dxcar(2,c)="
!     -       ,c,dxsph_dxcar(2,c)
!          end do
!          do c=1,4
!            write (*,*) "c,dxsph_dxcar(3,c)="
!     -       ,c,dxsph_dxcar(3,c)
!          end do
!          do c=1,4
!            write (*,*) "c,dxsph_dxcar(4,c)="
!     -       ,c,dxsph_dxcar(4,c)
!          end do
!         do b=1,4
!          do c=1,4
!            write (*,*) "b,c,d2xsph_dxcardxcard(1,b,c)="
!     -       ,b,c,d2xsph_dxcardxcar(1,b,c)
!          end do
!         end do
!         do b=1,4
!          do c=1,4
!            write (*,*) "b,c,d2xsph_dxcardxcard(2,b,c)="
!     -       ,b,c,d2xsph_dxcardxcar(2,b,c)
!          end do
!         end do
!         do b=1,4
!          do c=1,4
!            write (*,*) "b,c,d2xsph_dxcardxcard(3,b,c)="
!     -       ,b,c,d2xsph_dxcardxcar(3,b,c)
!          end do
!         end do
!         do b=1,4
!          do c=1,4
!            write (*,*) "b,c,d2xsph_dxcardxcard(4,b,c)="
!     -       ,b,c,d2xsph_dxcardxcar(4,b,c)
!          end do
!         end do
!        do a=1,4
!         do b=1,4
!          do c=1,4
!            write (*,*) "a,b,c,d3xsph_dxcardxcardxcar(1,a,b,c)="
!     -       ,a,b,c,d3xsph_dxcardxcardxcar(1,a,b,c)
!          end do
!         end do
!        end do
!        do a=1,4
!         do b=1,4
!          do c=1,4
!            write (*,*) "a,b,c,d3xsph_dxcardxcardxcar(2,a,b,c)="
!     -       ,a,b,c,d3xsph_dxcardxcardxcar(2,a,b,c)
!          end do
!         end do
!        end do
!        do a=1,4
!         do b=1,4
!          do c=1,4
!            write (*,*) "a,b,c,d3xsph_dxcardxcardxcar(3,a,b,c)="
!     -       ,a,b,c,d3xsph_dxcardxcardxcar(3,a,b,c)
!          end do
!         end do
!        end do
!        do a=1,4
!         do b=1,4
!          do c=1,4
!            write (*,*) "a,b,c,d3xsph_dxcardxcardxcar(4,a,b,c)="
!     -       ,a,b,c,d3xsph_dxcardxcardxcar(4,a,b,c)
!          end do
!         end do
!        end do
!        do a=1,4
!         do b=1,4
!             write (*,*) "a,b,gads_ll(a,b)="
!     -                   ,a,b,gads_ll(a,b)
!             write (*,*) "a,b,gads_uu(a,b)="
!     -                   ,a,b,gads_uu(a,b)
!          do c=1,4
!             write (*,*) "a,b,c,gads_ll_x(a,b,c)="
!     -                   ,a,b,c,gads_ll_x(a,b,c)
!             write (*,*) "a,b,c,gads_ll_sph_x(a,b,c)="
!     -                   ,a,b,c,gads_ll_sph_x(a,b,c)
!
!           do d=1,4
!             write (*,*) "a,b,c,d,gads_ll_xx(a,b,c,d)="
!     -                   ,a,b,c,d,gads_ll_xx(a,b,c,d)
!           end do
!          end do
!         end do
!        end do
!		 do b=1,4
!          do c=1,4
!            write (*,*) "b,c,
!     -       einsteinads_ll(b,c)+ Lambda* gads_ll(a,b)
!     -		 setads_ll(b,c)="
!     -       ,b,c,
!     -		  einsteinads_ll(b,c)-3*gads_ll(b,c),
!     -        setads_ll(b,c)
!          end do
!         end do
!		 do b=1,4
!          do c=1,4
!            write (*,*) "b,c,ricciads_ll(b,c)="
!     -       ,b,c,ricciads_ll(b,c)
!          end do
!         end do
!         write (*,*) "ricciads=",ricciads
!          stop
!        end if
!
!        if ((abs(x0-(-0.9375d0)).lt.10.0d0**(-10))
!     -  .and.(abs(y0-(-0.0d0)).lt.10.0d0**(-10))
!     -  .and.(abs(z0-(-0.0d0)).lt.10.0d0**(-10))) then
!        do a=1,4
!         do b=1,4
!             write (*,*) "a,b,gads_ll(a,b)="
!     -                   ,a,b,gads_ll(a,b)
!             write (*,*) "a,b,gads_uu(a,b)="
!     -                   ,a,b,gads_uu(a,b)
!          do c=1,4
!             write (*,*) "a,b,c,gads_ll_x(a,b,c)="
!     -                   ,a,b,c,gads_ll_x(a,b,c)
!             write (*,*) "a,b,c,gads_ll_sph_x(a,b,c)="
!     -                   ,a,b,c,gads_ll_sph_x(a,b,c)
!
!           do d=1,4
!             write (*,*) "a,b,c,d,gads_ll_xx(a,b,c,d)="
!     -                   ,a,b,c,d,gads_ll_xx(a,b,c,d)
!           end do
!          end do
!         end do
!        end do
!        stop
!        end if
!
!        if ( (is_nan(gads_ll(1,1))).or.
!     -       (is_nan(gads_ll(1,2))).or.
!     -       (is_nan(gads_ll(1,3))).or.
!     -       (is_nan(gads_ll(1,4))).or.
!     -       (is_nan(gads_ll(2,2))).or.
!     -       (is_nan(gads_ll(2,3))).or.
!     -       (is_nan(gads_ll(2,4))).or.
!     -       (is_nan(gads_ll(3,3))).or.
!     -       (is_nan(gads_ll(3,4))).or.
!     -       (is_nan(gads_ll(4,4)))    ) then
!!
!
!
!        if ((abs(x0-(-0.9375d0)).lt.10.0d0**(-10))
!     -  .and.(abs(y0-(-0.0625d0)).lt.10.0d0**(-10))
!     -  .and.(abs(z0-(-0.0d0)).lt.10.0d0**(-10))) then
!
!        if ((abs(Hads_l(1)).gt.10.0d0**(-10))
!     & .or.(abs(Hads_l(2)-
!     &  (2*x0*(-L**2*(5+x0**4+2*y0**2+2*z0**2+(y0**2+z0**2)**2
!     &            +2*x0**2*(1+y0**2+z0**2))*(-1+rho0**2)
!     &            -4*(2+rho0**2+rho0**4)))
!     &            /((-1+rho0**2)*(1+rho0**2)
!     &            *(4*rho0**2+L**2*(-1+rho0**2)**2)))
!     &  .gt.10.0d0**(-10))
!     & .or.(abs(Hads_l(3)-
!     & (2*y0*(-L**2*(5+x0**4+2*y0**2+2*z0**2+(y0**2+z0**2)**2
!     &            +2*x0**2*(1+y0**2+z0**2))*(-1+rho0**2)
!     &            -4*(2+rho0**2+rho0**4)))
!     &            /((-1+rho0**2)*(1+rho0**2)
!     &            *(4*rho0**2+L**2*(-1+rho0**2)**2)))
!     &      .gt.10.0d0**(-10))
!     & .or.(abs(Hads_l(4)-(
!     &  -((4*z0)/(-1+rho0**2))-(2*z0)/(1+rho0**2)
!     &            +(4*z0*(4+L**2*(-3+rho0**2)))
!     &            /(4*rho0**2+L**2*(-1+rho0**2)**2)))
!     &  .gt.10.0d0**(-10))) then
!         write (*,*) "ERROR in source functions
!         write (*,*) "x0,y0,z0=",x0,y0,z0
!         write (*,*) "rho0,theta0,phi0=",rho0,theta0,phi0
!            do a=1,4
!             write (*,*) 'a,logsqrtminusdetgads_l(a)='
!     &        ,a,logsqrtminusdetgads_l(a)
!            end do
!            do a=1,4
!             write (*,*) 'a,Hads_l(a)=',a,Hads_l(a)
!            end do
!        end if
!        stop
!!!!!!!!!!!!

		
        return
        end
