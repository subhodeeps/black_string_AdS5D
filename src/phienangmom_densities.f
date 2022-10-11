c-----------------------------------------------------------------------
c Calculate energy density and angular momentum density of the (bulk) scalar field, that is, respectively
c sqrt(gamma) rhoEphi=sqrt(gamma)(-alpha T^t_t) and sqrt(gamma) rhoJphi=sqrt(gamma)(-(1/2))(-alpha T^t_phi),
c where gamma is the determinant of the metric on t=const. slices in Cartesian coordinates.
c These give a measure of the energy and angular momentum of the scalar field 
c when integrated over dx dy dz. 
c The reason why we actually output the square root of these quantities is that
c we want to compute a quantity proportional to energy and angular momentum 
c via the L^2-norm calculation in DV (we then have to square the result of the L^2-norm).
c-----------------------------------------------------------------------

       subroutine sqrt_phienangmom_densities(
     &                  sqrtphiendensity,sqrtphiangmomdensity,
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                  gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                  Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  AH_R,AH_xc,
     &                  AH_Ntheta,AH_Nphi,
     &                  x,y,z,dt,ct,chr,L,ex,Nx,Ny,Nz,
     &                  phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background)

        implicit none

        logical is_nan

        real*8  ief_bh_r0,a_rot,M0,M0_min
        real*8  kerrads_background

        logical calc_der,calc_adv_quant
        data calc_der/.true./
        data calc_adv_quant/.true./

        integer Nx,Ny,Nz
        integer i,j,k
        integer phys_bdy(6),ghost_width(6)
        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,ct,L
        real*8 lambda4
        integer AH_Ntheta,AH_Nphi
        real*8 AH_R(AH_Ntheta,AH_Nphi)
        real*8 AH_xc(3)
        real*8 xp,yp,zp,rhop,thetap,phip
        integer i0,j0
        real*8 dahtheta,dahphi
        real*8 ft,fp
        real*8 AH_Rp,rhostop
        real*8 k_u(4),m3_u(4)
        real*8 sqrtphiendensity(Nx,Ny,Nz)
        real*8 sqrtphiangmomdensity(Nx,Ny,Nz)
        real*8 phiendensity,phiangmomdensity

        logical ltrace
        parameter (ltrace=.false.)


        integer is,ie,js,je,ks,ke

        integer i1,j1,k1,a,b,c,d,e,f,g,h,p,q,r
        integer ic,jc,kc
        real*8 efe_ires(4,4)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0   
        real*8 Rad0
        real*8 drho_dRad,dRad_drho

        real*8 dxcar_dxqssph(4,4)

        real*8 gb_tt_np1(Nx,Ny,Nz),gb_tt_n(Nx,Ny,Nz),gb_tt_nm1(Nx,Ny,Nz)
        real*8 gb_tx_np1(Nx,Ny,Nz),gb_tx_n(Nx,Ny,Nz),gb_tx_nm1(Nx,Ny,Nz)
        real*8 gb_ty_np1(Nx,Ny,Nz),gb_ty_n(Nx,Ny,Nz),gb_ty_nm1(Nx,Ny,Nz)
        real*8 gb_tz_np1(Nx,Ny,Nz),gb_tz_n(Nx,Ny,Nz),gb_tz_nm1(Nx,Ny,Nz)
        real*8 gb_xx_np1(Nx,Ny,Nz),gb_xx_n(Nx,Ny,Nz),gb_xx_nm1(Nx,Ny,Nz)
        real*8 gb_xy_np1(Nx,Ny,Nz),gb_xy_n(Nx,Ny,Nz),gb_xy_nm1(Nx,Ny,Nz)
        real*8 gb_xz_np1(Nx,Ny,Nz),gb_xz_n(Nx,Ny,Nz),gb_xz_nm1(Nx,Ny,Nz)
        real*8 gb_yy_np1(Nx,Ny,Nz),gb_yy_n(Nx,Ny,Nz),gb_yy_nm1(Nx,Ny,Nz)
        real*8 gb_yz_np1(Nx,Ny,Nz),gb_yz_n(Nx,Ny,Nz),gb_yz_nm1(Nx,Ny,Nz)
        real*8 gb_zz_np1(Nx,Ny,Nz),gb_zz_n(Nx,Ny,Nz),gb_zz_nm1(Nx,Ny,Nz)
        real*8 Hb_t_np1(Nx,Ny,Nz),Hb_t_n(Nx,Ny,Nz),Hb_t_nm1(Nx,Ny,Nz)
        real*8 Hb_x_np1(Nx,Ny,Nz),Hb_x_n(Nx,Ny,Nz),Hb_x_nm1(Nx,Ny,Nz)
        real*8 Hb_y_np1(Nx,Ny,Nz),Hb_y_n(Nx,Ny,Nz),Hb_y_nm1(Nx,Ny,Nz)
        real*8 Hb_z_np1(Nx,Ny,Nz),Hb_z_n(Nx,Ny,Nz),Hb_z_nm1(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,theta,phi)
        !--------------------------------------------------------------
        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4)
        real*8 gads_ll_xx(4,4,4,4)
        real*8 Hads_l(4)
        real*8 phi1ads, phi1ads_x(4)

        real*8 g0_ll(4,4),g0_uu(4,4),detg0
        real*8 g0_ll_x(4,4,4),g0_uu_x(4,4,4),g0_ll_xx(4,4,4,4)
        real*8 h0_ll(4,4),h0_uu(4,4)
        real*8 h0_ll_x(4,4,4),h0_uu_x(4,4,4),h0_ll_xx(4,4,4,4)
        real*8 gamma_ull(4,4,4),gamma_ull_x(4,4,4,4)
        real*8 riemann_ulll(4,4,4,4)
        real*8 ricci_ll(4,4),ricci_lu(4,4),ricci
        real*8 einstein_ll(4,4),set_ll(4,4)
        real*8 A_l(4),A_l_x(4,4)
        real*8 phi10_x(4),phi10_xx(4,4)

        real*8 detg0_qssph0,detgamma,lapse
        real*8 g0_ll_qssph(4,4),g0_uu_qssph(4,4)
        real*8 set_ll_qssph(4,4),set_ul_qssph(4,4)
        real*8 set_ul(4,4)



        ! initialize fixed-size variables
        data i,j,k,is,ie,js,je,ks,ke/0,0,0,0,0,0,0,0,0/
        data ic,jc,kc/0,0,0/
        data i1,j1,k1,a,b,c,d,e,p,q,r/0,0,0,0,0,0,0,0,0,0,0/

        data dx,dy,dz/0.0,0.0,0.0/
        data x0,y0,z0,rho0/0.0,0.0,0.0,0.0/    

        data gads_ll,gads_uu/16*0.0,16*0.0/
        data gads_ll_x,gads_uu_x/64*0.0,64*0.0/
        data gads_ll_xx/256*0.0/
        data Hads_l/4*0.0/
        data phi1ads_x/4*0.0/

        data g0_ll,g0_uu/16*0.0,16*0.0/
        data g0_ll_x,g0_uu_x/64*0.0,64*0.0/
        data g0_ll_xx/256*0.0/
        data h0_ll,h0_uu/16*0.0,16*0.0/
        data h0_ll_x,h0_uu_x/64*0.0,64*0.0/
        data h0_ll_xx/256*0.0/
        data gamma_ull/64*0.0/
        data gamma_ull_x/256*0.0/
        data riemann_ulll/256*0.0/
        data ricci_ll,ricci_lu/16*0.0,16*0.0/
        data einstein_ll,set_ll/16*0.0,16*0.0/
        data A_l,A_l_x/4*0.0,16*0.0/
        data phi10_x,phi10_xx/4*0.0,16*0.0/
        data dxcar_dxqssph/16*0.0/

        data g0_ll_qssph/16*0.0/
        data g0_uu_qssph/16*0.0/
        data set_ll_qssph/16*0.0/
        data set_ul_qssph/16*0.0/
        data set_ul/16*0.0/



!----------------------------------------------------------------------


      ! Black hole mass
        M0=ief_bh_r0/2
      ! Minimum black hole mass. For M0 below this value, there is a naked singularity
        M0_min=((2*(1 + a_rot**2/L**2) + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2))*Sqrt(-1 + Sqrt((1 + a_rot**2/L**2)**2 
     &   + (12*a_rot**2)/L**2) - a_rot**2/L**2))/(3.*Sqrt(6.))

        if (a_rot.ge.L) then
         write (*,*) "sqrth1spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &   the rotation parameter a must be smaller than the AdS radius L"
          write (*,*) "a_rot,L=",a_rot,L
          stop
        end if

        if ((abs(M0).gt.10.0d0**(-10))
     &     .and.(M0.le.M0_min)) then
          write (*,*) "sqrth1spnormdensity_phi: ERROR in choice of Kerr-AdS initial parameters: 
     &      the black hole mass M0=2*r0 must be larger
     &      than the M0_min value"
          write (*,*) "M0,M0_min=",M0,M0_min
          stop
        end if


        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)


        ! (MAIN LOOP) loop through spacetime points x(i),y(j)
        do i=is,ie
         do j=js,je
          do k=ks,ke

              x0=x(i)
              y0=y(j)
              z0=z(k)
              rho0=sqrt(x0**2+y0**2+z0**2)
              theta0=acos(x0/rho0)
              if (z0.lt.0) then
               phi0=atan2(z0,y0)+2*PI
              else
               phi0=atan2(z0,y0)
              end if

              !determine radius at which we stop computing the norms. This is the apparent horizon (AH) radius.

              !Cartesian coordinates of grid point p=(i,j,k) w.r.t. the centre of the AH
              xp=x0-AH_xc(1)
              yp=y0-AH_xc(2)
              zp=z0-AH_xc(3)

              !The center of the AH xp=yp=zp=0 will be excluded or included depending on the excision mask chr(i,j,k), 
              !so no need to force inclusion or exclusion via the value of rhostop. 
              !In fact, xp=yp=zp=0 is also the centre of the grid if there is no AH (also if the AH is centred at the origin), 
              !so we should include it if we want to compute these norms in the case of no horizon.
              !For this reason, it's best to leave the inclusion or exclusion of that point to the excision mask chr(i,j,k).
              if ((abs(xp).lt.10.0d0**(-12)).and.
     &            (abs(yp).lt.10.0d0**(-12)).and.  
     &            (abs(zp).lt.10.0d0**(-12)) ) then
                rhostop=-1
              else
                !corresponding spherical coordinates of point p
                !thetap,phip is in general NOT on the AH_Ntheta,AH_Nphi grid.
                rhop=sqrt(xp**2+yp**2+zp**2)
                thetap=acos(xp/rhop)
                if (zp.lt.0) then
                    phip=atan2(zp,yp)+2*PI
                else
                    phip=atan2(zp,yp)
                end if

                ! i0 is the index of the point with chi value closest to chip on the AH_Nchi grid (indices going from 0 to AH_Ntheta-1),
                ! and larger than chip (except if thetap=Pi, then it is smaller)
                ! similarly for j0
                dahtheta=PI/(AH_Ntheta-1)
                dahphi=2*PI/(AH_Nphi-1)
                i0=floor(thetap/dahtheta+1)
                if (AH_Ntheta-1<i0) i0=AH_Ntheta-1
                j0=floor(phip/dahphi+1)
                if (AH_Nphi-1<j0) j0=AH_Nphi-1
                ft=((thetap-(i0-1)*dahtheta))/dahtheta
                fp=((phip-(j0-1)*dahphi))/dahphi

                !bi-linear interpolated value of AH radius at (thetap,phip)
                AH_Rp=  AH_R(i0,j0)*(1-ft)*(1-fp)+
     &                  AH_R(i0+1,j0)*(ft)*(1-fp)+
     &                  AH_R(i0,j0+1)*(1-ft)*(fp)+
     &                  AH_R(i0+1,j0+1)*(ft)*(fp)

                if (AH_Rp.lt.1) then 
                    rhostop=AH_Rp
                else
                    rhostop=-1
                end if
              end if



            if ((rho0.gt.rhostop).and.
     &         (chr(i,j,k).ne.ex)) then

    		  !compute metric and its derivatives 
              call tensor_init(
     &                gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                gb_tz_np1,gb_tz_n,gb_tz_nm1,
     &                gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                gb_xz_np1,gb_xz_n,gb_xz_nm1,
     &                gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                gb_yz_np1,gb_yz_n,gb_yz_nm1,
     &                gb_zz_np1,gb_zz_n,gb_zz_nm1,
     &                Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                Hb_z_np1,Hb_z_n,Hb_z_nm1,
     &                phi1_np1,phi1_n,phi1_nm1,
     &                g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                gads_ll,gads_uu,
     &                gads_ll_x,gads_uu_x,gads_ll_xx,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                A_l,A_l_x,Hads_l,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                einstein_ll,set_ll,
     &                phi10_x,phi10_xx,
     &                x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                ief_bh_r0,a_rot,kerrads_background,
     &                calc_der,calc_adv_quant)

             do a=1,4
               do b=1,4
                set_ul(a,b)=0.0d0
                do c=1,4
                  set_ul(a,b)=set_ul(a,b)
     &                    +g0_uu(a,c)*set_ll(c,b)
                end do
               end do
             end do

             lapse=1/sqrt(-g0_uu(1,1))


            !define generators of asymptotic time translations and asymptotic rotations around the x-axis: d/dt, d/dphi=-z(d/dy)+y(d/dz)
            k_u(1)=1
            k_u(2)=0
            k_u(3)=0
            k_u(4)=0

            m3_u(1)=0
            m3_u(2)=0
            m3_u(3)=-z0
            m3_u(4)=y0

            phiendensity=0.0d0
            phiangmomdensity=0.0d0
            do a=1,4
                phiendensity=phiendensity+
     &           (-lapse)*set_ul(1,a)*k_u(a)
                phiangmomdensity=phiangmomdensity+
     &           (-1.0d0)*(-lapse)*set_ul(1,a)*m3_u(a)
            end do

             detgamma=-g0_ll(2,4)**2*g0_ll(3,3)
     &          +2*g0_ll(2,3)*g0_ll(2,4)*g0_ll(3,4) 
     &          -g0_ll(2,2)*g0_ll(3,4)**2
     &          -g0_ll(2,3)**2*g0_ll(4,4)
     &          +g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)

            phiendensity=sqrt(detgamma)*phiendensity
            phiangmomdensity=sqrt(detgamma)*phiangmomdensity

            sqrtphiendensity(i,j,k)=sqrt(phiendensity)
            sqrtphiangmomdensity(i,j,k)=phiangmomdensity


     		else !i.e., the point is either excised or inside the apparent horizon

     			sqrtphiendensity(i,j,k)=0.0d0
             	sqrtphiangmomdensity(i,j,k)=0.0d0

     		end if


          end do
         end do
        end do


        return
        end

