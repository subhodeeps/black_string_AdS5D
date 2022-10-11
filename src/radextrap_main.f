c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for setting the mask chrbdy to a value different 
c from ex at the outermost grid points that will then be used 
c for radial extrapolation.  
c
c The extrapolation technique goes as follows:
c 1. define a range in which we look for grid points  
c to use for extrapolation. Consider the outermost  
c points in this range. 
c 2. Take one of such points, denoted by p1.  
c This will be the vertex of a cube (the other vertices  
c are other grid points). Consider the radial direction  
c connecting the point to the origin. This direction  
c intersects a face of the cube at a point p2. This will be  
c the second point used for extrapolation. The angular coords  
c of p2 are the same as p1, the radial coord must be obtained  
c and it depends on the quadrant of the 3-dimensional grid in  
c which p1 is.
c 3. To obtain the value of the function at p2, we use bilinear  
c extrapolation from the 4 vertices of the face of the cube (if  
c these are not accessible from the current process, we discard  
c the point p1 and move on to the next one).
c----------------------------------------------------------------------

        subroutine nexttobdypoints_radextrap(
     &                  chrbdy,
     &                  numbdypoints,
     &                  bdy_extrap_order,
     &                  currentres_ratio_Lhighres_Llowres,
     &                  half_steps_from_bdy_ext,
     &                  half_steps_from_bdy_int,
     &                  x,y,z,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz),chrbdy2(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer ip2a,jp2a,kp2a
        integer ip2b,jp2b,kp2b
        integer ip2c,jp2c,kp2c
        integer ip2d,jp2d,kp2d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1,rhop1,chip1,xip1
        real*8 maxxyzp1
        integer numbdypoints
        integer half_steps_from_bdy_ext
        integer half_steps_from_bdy_int
        real*8 currentres_ratio_Lhighres_Llowres
        integer bdy_extrap_order

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        numbdypoints=0

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

         do i=1,Nx
          do j=1,Ny
           do k=1,Nz
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)
           chip1=(1/PI)*acos(xp1/rhop1)
           if (zp1.lt.0) then
              xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
           else
              xip1=(1/(2*PI))*atan2(zp1,yp1)
           end if

!select only points that have radius smaller than (1-5*dx/2).
!For each resolution, points between rhobdy=1 and 1-dx/2 are excised, points between rhobdy=1 and 1-3*dx/2 use forward/backward stencils (so when we cannot expect convergence at these points because the stencils used are different for different resolutions), points between rhobdy=1 and 1-5*dx/2 use points that are set by forward/backward stencils, so we cannot expect convergence. If we want to check convergence at the grid points used for extrapolation at the boundary, we need to pick points that have rho<1-5*dx/2. However, we see that we get good convergence only if we use points with rho<1-9*dx/2. This should be understood in the future. The condition rhop1>currentres_ratio_Lhighres_Llowres ensures that the first suitable point for extrapolation is not too far from the AdS boundary.
            if ((rhop1.lt.(1-half_steps_from_bdy_ext*dx/2))
     &          .and.(rhop1.gt.
     &                (1-half_steps_from_bdy_int
     &                  *currentres_ratio_Lhighres_Llowres*dx/2))
     &          .and.(chr(i,j,k).ne.ex)
     &         ) then
              chrbdy(i,j,k) =ex-1.0d0
              chrbdy2(i,j,k)=ex-1.0d0
            else
              chrbdy(i,j,k) =ex
              chrbdy2(i,j,k)=ex
            end if


            if ((i.lt.is).or.(i.gt.ie).or.
     &          (j.lt.js).or.(j.gt.je).or.
     &          (k.lt.ks).or.(k.gt.ke)) then

               chrbdy(i,j,k) =ex
               chrbdy2(i,j,k)=ex
            end if


!! eliminate troublesome points
!!removing points with abs(x)=abs(y)=abs(z): even if the grid function quasiset_tracell converges at these points, the extrapolated boundary quantity quasisettrace, obtained using these points, is not converging. It would be interesting to understand why. For now, we just eliminate these points from the set of points used for extrapolation
!            if (
!     &          (abs(abs(xp1)-abs(yp1)).lt.10.0d0**(-10)).and.
!     &          (abs(abs(xp1)-abs(zp1)).lt.10.0d0**(-10))
!     &         ) then
!                  chrbdy(i,j,k)=ex
!                  chrbdy2(i,j,k)=ex
!            end if




! eliminate troublesome points
!y=z=0 (i.e. the x axis, where chi=0 or 1) are troublesome points where the xi coordinate is not defined  
!We will fill and impose regularity at these points in post-processing in Mathematica.
            if (
     &          (abs(yp1).lt.10.0d0**(-10)).and.
     &          (abs(zp1).lt.10.0d0**(-10)) 
     &         ) then
             chrbdy(i,j,k)=ex
             chrbdy2(i,j,k)=ex
            end if

           end do
          end do
         end do


         do i=is,ie
          do j=js,je
           do k=ks,ke

           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)
           chip1=(1/PI)*acos(xp1/rhop1)
           if (zp1.lt.0) then
              xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
           else
              xip1=(1/(2*PI))*atan2(zp1,yp1)
           end if


           if (chrbdy(i,j,k).ne.ex) then


            if (bdy_extrap_order.eq.1) then


! if some points are not suitable for radial extrapolation 
! (typically because either their neighbouring points needed for extrapolation are not available
! or because they are not the outermost points that are suitable for radial extrapolation)
! the following routine sets chrbdy(i,j,k)=ex at these points.
              call firstord_chrbdy_radextrap(
     &                  chrbdy,
     &                  chrbdy2,
     &                  is,ie,js,je,ks,ke,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  chr,ex,Nx,Ny,Nz)
           

            end if !closes condition on bdy_extrap_order.eq.1

            if (bdy_extrap_order.eq.2) then


! if some points are not suitable for radial extrapolation 
! (typically because either their neighbouring points needed for extrapolation are not available
! or because they are not the outermost points that are suitable for radial extrapolation)
! the following routine sets chrbdy(i,j,k)=ex at these points.
              call secondord_chrbdy_radextrap(
     &                  chrbdy,
     &                  chrbdy2,
     &                  is,ie,js,je,ks,ke,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  chr,ex,Nx,Ny,Nz)
           

            end if !closes condition on bdy_extrap_order.eq.2


          end if !closes condition on chrbdy(i,j,k).ne.ex





          if (chrbdy(i,j,k).ne.ex) then
            numbdypoints=numbdypoints+1
          end if


         end do
        end do
      end do

        return
        end


c------------------------------------------------------------------------------------------------------


c------------------------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the Cartesian coordinates of each AdS boundary point pbdy where we want to extrapolate the value of functions 
c AND the coordinates of the outermost point pout used for radial extrapolation for each of the boundary points
c-------------------------------------------------------------------------------------------------------------------------

        subroutine xyz_bdy_out_radextrap(
     &                  xpbdy,ypbdy,zpbdy,
     &                  xpout,ypout,zpout,
     &                  chrbdy,
     &                  numbdypoints,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz
        integer lind

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1
        real*8 rhop1,chip1,xip1
        real*8 rhopbdy,chipbdy,xipbdy
        real*8 maxxyzp1
        integer numbdypoints
        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)
        real*8 xpout(numbdypoints)
        real*8 ypout(numbdypoints)
        real*8 zpout(numbdypoints)

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

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

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke

            if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

              xp1=x(i)
              yp1=y(j)
              zp1=z(k)
              rhop1=sqrt(xp1**2+yp1**2+zp1**2)
              chip1=(1/PI)*acos(xp1/rhop1)
              if (zp1.lt.0) then
                xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
              else
                xip1=(1/(2*PI))*atan2(zp1,yp1)
              end if

              xpout(lind)=xp1
              ypout(lind)=yp1
              zpout(lind)=zp1

              rhopbdy=1
              chipbdy=chip1
              xipbdy=xip1

              xpbdy(lind)=rhopbdy*cos(PI*chipbdy)
              ypbdy(lind)=rhopbdy*sin(PI*chipbdy)
     &              *cos(2*PI*xipbdy)
              zpbdy(lind)=rhopbdy*sin(PI*chipbdy)
     &              *sin(2*PI*xipbdy)


            end if !closes condition on (chrbdy(i,j,k).ne.ex)
          end do
         end do
        end do

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the scalar field phi of the boundary CFT
c----------------------------------------------------------------------

        subroutine bdyphi_radextrap(bdyphi,
     &                  leadordcoeff_phi1,
     &                  xpbdy,ypbdy,zpbdy,
     &                  chrbdy,numbdypoints,
     &                  bdy_extrap_order,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        integer bdy_extrap_order
        real*8 chrbdy(Nx,Ny,Nz)
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        integer i,j,k,is,ie,js,je,ks,ke,lind
        integer ip2a,jp2a,kp2a
        integer ip2b,jp2b,kp2b
        integer ip2c,jp2c,kp2c
        integer ip2d,jp2d,kp2d
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0
        real*8 xp2a,xp2b,xp2c,xp2d
        real*8 yp2a,yp2b,yp2c,yp2d
        real*8 zp2a,zp2b,zp2c,zp2d

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 leadordcoeff_phi1(Nx,Ny,Nz)
        real*8 leadordcoeff_phi1_p1

        real*8 leadordcoeff_phi1_p2a
        real*8 leadordcoeff_phi1_p2b
        real*8 leadordcoeff_phi1_p2c
        real*8 leadordcoeff_phi1_p2d

        real*8 leadordcoeff_phi1_p2
        real*8 leadordcoeff_phi1_p3
        real*8 leadordcoeff_phi1_p4
        real*8 bdyphi(numbdypoints)

        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 xp1,yp1,zp1
        real*8 xp2,yp2,zp2
        real*8 xp3,yp3,zp3
        real*8 xp4,yp4,zp4
        real*8 xex,yex,zex,rhoex,chiex,xiex
        real*8 rhop1,chip1,xip1
        real*8 rhop2,chip2,xip2
        real*8 maxxyzp1

        real*8 bilinear_interp
        real*8 firstord_extrap
        real*8 secondord_extrap
        real*8 thirdord_extrap
        real*8 fourthord_extrap

        real*8 firstord_func_radextrap
        real*8 secondord_func_radextrap
!----------------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)


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

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke

           if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

              xp1=x(i)
              yp1=y(j)
              zp1=z(k)
              rhop1=sqrt(xp1**2+yp1**2+zp1**2)
              chip1=(1/PI)*acos(xp1/rhop1)
              if (zp1.lt.0) then
                xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
              else
                xip1=(1/(2*PI))*atan2(zp1,yp1)
              end if
  
              leadordcoeff_phi1_p1=leadordcoeff_phi1(i,j,k)

              xex=xpbdy(lind)
              yex=ypbdy(lind)
              zex=zpbdy(lind)
              rhoex=1
              chiex=(1/PI)*acos(xex/rhoex)
              if (zex.lt.0) then
                xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
              else
                xiex=(1/(2*PI))*atan2(zex,yex)
              end if



            if (bdy_extrap_order.eq.1) then

                      bdyphi(lind)=
     &                   firstord_func_radextrap(
     &                  leadordcoeff_phi1,
     &                  leadordcoeff_phi1_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)




!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.

            if (bdy_extrap_order.eq.2) then

                      bdyphi(lind)=
     &                   secondord_func_radextrap(
     &                  leadordcoeff_phi1,
     &                  leadordcoeff_phi1_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)




!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.2


           end if !closes condition on chrbdy(i,j,k).ne.ex




          end do
         end do
        end do


        return
        end
c--------------------------------------------------------------------------------------

c-------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D_polar  
c using a 1-rho expansion about rho=1 AT the boundary through extrapolation.
c The tensor components are given in spherical polar coordinates.
c-------------------------------------------------------------------------------------

        subroutine quasiset_radextrap(
     &                  quasiset_tt,quasiset_tchi,quasiset_txi,
     &                  quasiset_chichi,quasiset_chixi,
     &                  quasiset_xixi,
     &                  quasiset_trace,
     &                  quasiset_massdensity,
     &                  quasiset_angmomdensityx,
     &                  quasiset_angmomdensityy,
     &                  quasiset_angmomdensityz,
     &                  quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
     &                  quasiset_chichi_ll,quasiset_chixi_ll,
     &                  quasiset_xixi_ll,
     &                  quasiset_tracell,
     &                  quasiset_massdensityll,
     &                  quasiset_angmomdensityxll,
     &                  quasiset_angmomdensityyll,
     &                  quasiset_angmomdensityzll,
     &                  xpbdy,ypbdy,zpbdy,
     &                  chrbdy,numbdypoints,
     &                  bdy_extrap_order,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        integer bdy_extrap_order
        real*8 chrbdy(Nx,Ny,Nz)
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
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex
        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 quasiset_tt_ll(Nx,Ny,Nz),quasiset_tchi_ll(Nx,Ny,Nz)
        real*8 quasiset_txi_ll(Nx,Ny,Nz),quasiset_chichi_ll(Nx,Ny,Nz)
        real*8 quasiset_chixi_ll(Nx,Ny,Nz),quasiset_xixi_ll(Nx,Ny,Nz)
        real*8 quasiset_tracell(Nx,Ny,Nz)
        real*8 quasiset_massdensityll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityxll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityyll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityzll(Nx,Ny,Nz)

        real*8 quasiset_tt_p1,quasiset_tchi_p1
        real*8 quasiset_txi_p1,quasiset_chichi_p1
        real*8 quasiset_chixi_p1,quasiset_xixi_p1
        real*8 quasiset_trace_p1
        real*8 quasiset_massdensity_p1
        real*8 quasiset_angmomdensityx_p1
        real*8 quasiset_angmomdensityy_p1
        real*8 quasiset_angmomdensityz_p1

        real*8 quasiset_tt_p2a,quasiset_tchi_p2a
        real*8 quasiset_txi_p2a,quasiset_chichi_p2a
        real*8 quasiset_chixi_p2a,quasiset_xixi_p2a
        real*8 quasiset_trace_p2a
        real*8 quasiset_massdensity_p2a
        real*8 quasiset_angmomdensityx_p2a
        real*8 quasiset_angmomdensityy_p2a
        real*8 quasiset_angmomdensityz_p2a

        real*8 quasiset_tt_p2b,quasiset_tchi_p2b
        real*8 quasiset_txi_p2b,quasiset_chichi_p2b
        real*8 quasiset_chixi_p2b,quasiset_xixi_p2b
        real*8 quasiset_trace_p2b
        real*8 quasiset_massdensity_p2b
        real*8 quasiset_angmomdensityx_p2b
        real*8 quasiset_angmomdensityy_p2b
        real*8 quasiset_angmomdensityz_p2b

        real*8 quasiset_tt_p2c,quasiset_tchi_p2c
        real*8 quasiset_txi_p2c,quasiset_chichi_p2c
        real*8 quasiset_chixi_p2c,quasiset_xixi_p2c
        real*8 quasiset_trace_p2c
        real*8 quasiset_massdensity_p2c
        real*8 quasiset_angmomdensityx_p2c
        real*8 quasiset_angmomdensityy_p2c
        real*8 quasiset_angmomdensityz_p2c

        real*8 quasiset_tt_p2d,quasiset_tchi_p2d
        real*8 quasiset_txi_p2d,quasiset_chichi_p2d
        real*8 quasiset_chixi_p2d,quasiset_xixi_p2d
        real*8 quasiset_trace_p2d
        real*8 quasiset_massdensity_p2d
        real*8 quasiset_angmomdensityx_p2d
        real*8 quasiset_angmomdensityy_p2d
        real*8 quasiset_angmomdensityz_p2d

        real*8 quasiset_tt_p2,quasiset_tchi_p2
        real*8 quasiset_txi_p2,quasiset_chichi_p2
        real*8 quasiset_chixi_p2,quasiset_xixi_p2
        real*8 quasiset_trace_p2
        real*8 quasiset_massdensity_p2
        real*8 quasiset_angmomdensityx_p2
        real*8 quasiset_angmomdensityy_p2
        real*8 quasiset_angmomdensityz_p2

        real*8 quasiset_tt_p3,quasiset_tchi_p3
        real*8 quasiset_txi_p3,quasiset_chichi_p3
        real*8 quasiset_chixi_p3,quasiset_xixi_p3
        real*8 quasiset_trace_p3
        real*8 quasiset_massdensity_p3
        real*8 quasiset_angmomdensityx_p3
        real*8 quasiset_angmomdensityy_p3
        real*8 quasiset_angmomdensityz_p3

        real*8 quasiset_tt_p4,quasiset_tchi_p4
        real*8 quasiset_txi_p4,quasiset_chichi_p4
        real*8 quasiset_chixi_p4,quasiset_xixi_p4
        real*8 quasiset_trace_p4
        real*8 quasiset_massdensity_p4
        real*8 quasiset_angmomdensityx_p4
        real*8 quasiset_angmomdensityy_p4
        real*8 quasiset_angmomdensityz_p4


        real*8 gamma0sphbdy_uu_tt
        real*8 gamma0sphbdy_uu_tchi
        real*8 gamma0sphbdy_uu_txi
        real*8 gamma0sphbdy_uu_chichi
        real*8 gamma0sphbdy_uu_chixi
        real*8 gamma0sphbdy_uu_xixi

        real*8 quasiset_tt(numbdypoints),quasiset_tchi(numbdypoints)
        real*8 quasiset_txi(numbdypoints),quasiset_chichi(numbdypoints)
        real*8 quasiset_chixi(numbdypoints),quasiset_xixi(numbdypoints)
        real*8 quasiset_trace(numbdypoints)
        real*8 quasiset_massdensity(numbdypoints)
        real*8 quasiset_angmomdensityx(numbdypoints)
        real*8 quasiset_angmomdensityy(numbdypoints)
        real*8 quasiset_angmomdensityz(numbdypoints)

        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 rhoextrap
        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        integer bdy_Nchi,bdy_Nxi

!        real*8 chibdy(bdy_Nchi)
!        real*8 xibdy(bdy_Nxi)

        integer i,j,k,is,ie,js,je,ks,ke,lind,m,e,increase

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        integer ip2a,jp2a,kp2a
        integer ip2b,jp2b,kp2b
        integer ip2c,jp2c,kp2c
        integer ip2d,jp2d,kp2d

        real*8 x0,y0,z0,rho0,q
        real*8 xp2a,xp2b,xp2c,xp2d
        real*8 yp2a,yp2b,yp2c,yp2d
        real*8 zp2a,zp2b,zp2c,zp2d

        real*8 xp1,yp1,zp1
        real*8 xp2,yp2,zp2
        real*8 xp3,yp3,zp3
        real*8 xp4,yp4,zp4
        real*8 xex,yex,zex,rhoex,chiex,xiex
        real*8 rhop1,chip1,xip1
        real*8 rhop2,chip2,xip2
        real*8 maxxyzp1

        real*8 bilinear_interp
        real*8 firstord_extrap
        real*8 secondord_extrap
        real*8 thirdord_extrap
        real*8 fourthord_extrap

        real*8 firstord_func_radextrap
        real*8 secondord_func_radextrap

        real*8 dp1p2

        real*8 AdS_mass

!----------------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)


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

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke

           if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

              xp1=x(i)
              yp1=y(j)
              zp1=z(k)
              rhop1=sqrt(xp1**2+yp1**2+zp1**2)
              chip1=(1/PI)*acos(xp1/rhop1)
              if (zp1.lt.0) then
                xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
              else
                xip1=(1/(2*PI))*atan2(zp1,yp1)
              end if
  
              quasiset_tt_p1=
     &                        quasiset_tt_ll(i,j,k)
              quasiset_tchi_p1=
     &                        quasiset_tchi_ll(i,j,k)
              quasiset_txi_p1=
     &                        quasiset_txi_ll(i,j,k)
              quasiset_chichi_p1=
     &                        quasiset_chichi_ll(i,j,k)
              quasiset_chixi_p1=
     &                        quasiset_chixi_ll(i,j,k)
              quasiset_xixi_p1=
     &                        quasiset_xixi_ll(i,j,k)
              quasiset_trace_p1=
     &                        quasiset_tracell(i,j,k)
              quasiset_massdensity_p1=
     &                        quasiset_massdensityll(i,j,k)
              quasiset_angmomdensityx_p1=
     &                        quasiset_angmomdensityxll(i,j,k)
              quasiset_angmomdensityy_p1=
     &                        quasiset_angmomdensityyll(i,j,k)
              quasiset_angmomdensityz_p1=
     &                        quasiset_angmomdensityzll(i,j,k)

              xex=xpbdy(lind)
              yex=ypbdy(lind)
              zex=zpbdy(lind)
              rhoex=1
              chiex=(1/PI)*acos(xex/rhoex)
              if (zex.lt.0) then
                xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
              else
                xiex=(1/(2*PI))*atan2(zex,yex)
              end if

              !inverse of conformal metric on AdS boundary (needed for trace) at extrapolated point
              gamma0sphbdy_uu_tt=-1
              gamma0sphbdy_uu_tchi=0
              gamma0sphbdy_uu_txi=0
              gamma0sphbdy_uu_chichi=1/(PI**2)
              gamma0sphbdy_uu_chixi=0
              gamma0sphbdy_uu_xixi=1/((sin(PI*chiex))**2)/4/PI**2



            if (bdy_extrap_order.eq.1) then


                      quasiset_tt(lind)=
     &                   firstord_func_radextrap(
     &                  quasiset_tt_ll,
     &                  quasiset_tt_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_tchi(lind)=
     &                   firstord_func_radextrap(
     &                  quasiset_tchi_ll,
     &                  quasiset_tchi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_txi(lind)=
     &                   firstord_func_radextrap(
     &                  quasiset_txi_ll,
     &                  quasiset_txi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_chichi(lind)=
     &                   firstord_func_radextrap(
     &                  quasiset_chichi_ll,
     &                  quasiset_chichi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_chixi(lind)=
     &                   firstord_func_radextrap(
     &                  quasiset_chixi_ll,
     &                  quasiset_chixi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_xixi(lind)=
     &                   firstord_func_radextrap(
     &                  quasiset_xixi_ll,
     &                  quasiset_xixi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)


                      quasiset_trace(lind)=
     &                      firstord_func_radextrap(
     &                  quasiset_tracell,
     &                  quasiset_trace_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )


                      quasiset_massdensity(lind)=
     &                      firstord_func_radextrap(
     &                  quasiset_massdensityll,
     &                  quasiset_massdensity_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_angmomdensityx(lind)=
     &                      firstord_func_radextrap(
     &                  quasiset_angmomdensityxll,
     &                  quasiset_angmomdensityx_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_angmomdensityy(lind)=
     &                      firstord_func_radextrap(
     &                  quasiset_angmomdensityyll,
     &                  quasiset_angmomdensityy_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_angmomdensityz(lind)=
     &                      firstord_func_radextrap(
     &                  quasiset_angmomdensityzll,
     &                  quasiset_angmomdensityz_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)


!                  quasiset_tt(lind)=quasiset_tt_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "quasiset_tt(lind)=",quasiset_tt(lind)

            end if !closes condition on bdy_extrap_order.eq.1


            if (bdy_extrap_order.eq.2) then


                      quasiset_tt(lind)=
     &                   secondord_func_radextrap(
     &                  quasiset_tt_ll,
     &                  quasiset_tt_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_tchi(lind)=
     &                   secondord_func_radextrap(
     &                  quasiset_tchi_ll,
     &                  quasiset_tchi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_txi(lind)=
     &                   secondord_func_radextrap(
     &                  quasiset_txi_ll,
     &                  quasiset_txi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_chichi(lind)=
     &                   secondord_func_radextrap(
     &                  quasiset_chichi_ll,
     &                  quasiset_chichi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_chixi(lind)=
     &                   secondord_func_radextrap(
     &                  quasiset_chixi_ll,
     &                  quasiset_chixi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_xixi(lind)=
     &                   secondord_func_radextrap(
     &                  quasiset_xixi_ll,
     &                  quasiset_xixi_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)


                      quasiset_trace(lind)=
     &                      secondord_func_radextrap(
     &                  quasiset_tracell,
     &                  quasiset_trace_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)
!                       quasiset_trace(lind)=
!     &                  (
!     &                  gamma0sphbdy_uu_tt*quasiset_tt(lind)
!     &              +2*gamma0sphbdy_uu_tchi*quasiset_tchi(lind)
!     &              +2*gamma0sphbdy_uu_txi*quasiset_txi(lind)
!     &                +gamma0sphbdy_uu_chichi*quasiset_chichi(lind)
!     &               +2*gamma0sphbdy_uu_chixi*quasiset_chixi(lind)
!     &                 +gamma0sphbdy_uu_xixi*quasiset_xixi(lind)
!     &                         )


                      quasiset_massdensity(lind)=
     &                      secondord_func_radextrap(
     &                  quasiset_massdensityll,
     &                  quasiset_massdensity_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_angmomdensityx(lind)=
     &                      secondord_func_radextrap(
     &                  quasiset_angmomdensityxll,
     &                  quasiset_angmomdensityx_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_angmomdensityy(lind)=
     &                      secondord_func_radextrap(
     &                  quasiset_angmomdensityyll,
     &                  quasiset_angmomdensityy_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

                      quasiset_angmomdensityz(lind)=
     &                      secondord_func_radextrap(
     &                  quasiset_angmomdensityzll,
     &                  quasiset_angmomdensityz_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)


!                  quasiset_tt(lind)=quasiset_tt_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "quasiset_tt(lind)=",quasiset_tt(lind)

            end if !closes condition on bdy_extrap_order.eq.2



           end if !closes condition on chrbdy(i,j,k).ne.ex




          end do
         end do
        end do


        return
        end
c--------------------------------------------------------------------------------------


c----------------------------------------------------------------------
c bilinear interpolation on a plane using values 
c Tp2a at (xp2a,yp2a), Tp2b at (xp2b,yp2b), Tp2c at (xp2c,yp2c), Tp2d at (xp2d,yp2d)  
c to obtain the value at (xp0,yp0) 
c----------------------------------------------------------------------
        real*8 function bilinear_interp(Tp2a,Tp2b,Tp2c,Tp2d,
     &                                  xp2a,xp2b,xp2c,xp2d,
     &                                  yp2a,yp2b,yp2c,yp2d,
     &                                  xp0,yp0)
        implicit none
        real*8 Tp2a,Tp2b,Tp2c,Tp2d
        real*8 xp2a,xp2b,xp2c,xp2d
        real*8 yp2a,yp2b,yp2c,yp2d
        real*8 xp0,yp0

        !--------------------------------------------------------------

        bilinear_interp=(1/((xp2c-xp2a)*(yp2c-yp2a)))*
     &    (Tp2a*(xp2c-xp0)*(yp2c-yp0)+
     &     Tp2b*(xp2c-xp0)*(yp0-yp2a)+
     &     Tp2c*(xp0-xp2a)*(yp0-yp2a)+
     &     Tp2d*(xp0-xp2a)*(yp2c-yp0))

        return
        end
c--------------------------------------------------------------------------------------