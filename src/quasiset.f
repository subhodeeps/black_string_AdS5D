c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for setting the mask chrbdy to a value different from ex at grid points that are at large rho=sqrt(x**2+y**2+z**2) and next to excised points, so we're identifying the last not excised grid points near the boundary. We will use them to extrapolate the value of the quasi-local boundary stress-energy tensor at the boundary
c----------------------------------------------------------------------

        subroutine nexttobdypoints_freepts(
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


           if (chrbdy(i,j,k).ne.ex) then


! eliminate troublesome points
!removing points with abs(x)=abs(y)=abs(z): even if the grid function quasiset_tracell is converging at these points, the extrapolated boundary quantity quasisettrace, obtained using these points, is not converging. It would be interesting to understand why. For now, we just eliminate these points from the set of points used for extrapolation
            if (
     &          (abs(abs(xp1)-abs(yp1)).lt.10.0d0**(-10)).and.
     &          (abs(abs(xp1)-abs(zp1)).lt.10.0d0**(-10))
     &         ) then
                  chrbdy(i,j,k)=ex
            end if




! eliminate troublesome points
!y=z=0 (i.e. the x axis, where chi=0 or 1) are troublesome points where the xi coordinate is not defined  
!We will fill and impose regularity at these points in post-processing in Mathematica.
            if (
     &          (abs(yp1).lt.10.0d0**(-10)).and.
     &          (abs(zp1).lt.10.0d0**(-10)) 
     &         ) then
             chrbdy(i,j,k)=ex
            end if


!If we use derivatives to define near boundary quantities, we will only define them at points between is and ie (js and je, ks and ke). Therefore, for extrapolation, we can only select near boundary points whose neighbors used for extrapolation in the direction of the bulk along the axes (i.e. the direction of extrapolation) are within that range
!We also need to make sure that those neighbours are not excised.
!We also define near boundary quantities at points where y0 and z0 are not both 0. So we need to make sure that we select points such that neighbouring points used for extrapolation don't have such values of y0,z0. Notice, we've already imposed that y(j)=!=0 and z(k)=!=0 for the outermost point. Therefore, when extrapolation is along y or z, we need to impose the condition that the (y,z) coordinates of the other points used are different from 0.
!The condition (chrbdy2(i+1,j,k).ne.ex) makes sure that (i,j,k) is the outmost point satisfying the conditions of the previous for-loop, which sets chrbdy2 as well as chrbdy. In other words, if there's an outer point w.r.t. (i,j,k) that satisfies those conditions, then we don't want to use (i,j,k) for extrapolation, but we will use that other point. 


          maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

           if (bdy_extrap_order.eq.1) then
            if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
             if (xp1.gt.0) then 
              if ((i-1).lt.is) then !it ensures that the index i-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i-1,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i+1).le.Nx) then !ensures that the index i+1 is not out of bounds
                 if (chrbdy2(i+1,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((i+1).gt.ie) then !it ensures that the index i+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i+1,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i-1).ge.1) then !ensures that the index i-1 is not out of bounds
                 if (chrbdy2(i-1,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if


            else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
             if (yp1.gt.0) then
              if ((j-1).lt.js) then !it ensures that the index j-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation  at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-1,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j+1).le.Ny) then !ensures that the index j+1 is not out of bounds
                 if (chrbdy2(i,j+1,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((j+1).gt.je) then !it ensures that the index j+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+1,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j-1).ge.1) then !ensures that the index j-1 is not out of bounds
                 if (chrbdy2(i,j-1,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            else !i.e. when maxxyzp1.eq.abs(zp1)
             if (zp1.gt.0) then
              if ((k-1).lt.ks) then !it ensures that the index k-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-1).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k+1).le.Nz) then !ensures that the index k+1 is not out of bounds
                 if (chrbdy2(i,j,k+1).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((k+1).gt.ke) then !it ensures that the index k+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+1).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k-1).ge.1) then !ensures that the index k-1 is not out of bounds
                 if (chrbdy2(i,j,k-1).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if
            
            end if !closes condition on (abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))
           

           end if !closes condition on bdy_extrap_order.eq.1


           if (bdy_extrap_order.eq.2) then
            if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
             if (xp1.gt.0) then
              if ((i-2).lt.is) then !it ensures that the index i-2,i-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i-1,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i-2,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i+1).le.Nx) then !ensures that the index i+1 is not out of bounds
                 if (chrbdy2(i+1,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((i+2).gt.ie) then !it ensures that the index i+2,i+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i+1,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i+2,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i-1).ge.1) then !ensures that the index i-1 is not out of bounds
                 if (chrbdy2(i-1,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if


            else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
             if (yp1.gt.0) then
              if ((j-2).lt.js) then !it ensures that the index j-2,j-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-1,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-2,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-2))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j+1).le.Ny) then !ensures that the index j+1 is not out of bounds
                 if (chrbdy2(i,j+1,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((j+2).gt.je) then !it ensures that the index j+2,j+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+1,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+2,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+2))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j-1).ge.1) then !ensures that the index j-1 is not out of bounds
                 if (chrbdy2(i,j-1,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            else !i.e. when maxxyzp1.eq.abs(zp1)
             if (zp1.gt.0) then
              if ((k-2).lt.ks) then !it ensures that the index k-2,k-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-1).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-2).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-2))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k+1).le.Nz) then !ensures that the index k+1 is not out of bounds
                 if (chrbdy2(i,j,k+1).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((k+2).gt.ke) then !it ensures that the index k+2,k+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+1).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+2).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+2))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k-1).ge.1) then !ensures that the index k-1 is not out of bounds
                 if (chrbdy2(i,j,k-1).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            end if !closes condition on (abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))


           end if !closes condition on bdy_extrap_order.eq.2



          if (bdy_extrap_order.eq.3) then
            if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
             if (xp1.gt.0) then
              if ((i-3).lt.is) then !it ensures that the index i-3,i-2,i-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i-1,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i-2,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i-3,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex     
              else if ((i+1).le.Nx) then !ensures that the index i+1 is not out of bounds
                 if (chrbdy2(i+1,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((i+3).gt.ie) then !it ensures that the index i+2,i+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i+1,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i+2,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i+3,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i-1).ge.1) then !ensures that the index i-1 is not out of bounds
                 if (chrbdy2(i-1,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if


            else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
             if (yp1.gt.0) then
              if ((j-3).lt.js) then !it ensures that the index j-2,j-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-1,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-2,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-3,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-2))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-3))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j+1).le.Ny) then !ensures that the index j+1 is not out of bounds
                 if (chrbdy2(i,j+1,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((j+3).gt.je) then !it ensures that the index j+2,j+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+1,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+2,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+3,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+2))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+3))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j-1).ge.1) then !ensures that the index j-1 is not out of bounds
                 if (chrbdy2(i,j-1,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            else !i.e. when maxxyzp1.eq.abs(zp1)
             if (zp1.gt.0) then
              if ((k-3).lt.ks) then !it ensures that the index k-2,k-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-1).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-2).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-3).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-2))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-3))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k+1).le.Nz) then !ensures that the index k+1 is not out of bounds
                 if (chrbdy2(i,j,k+1).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((k+3).gt.ke) then !it ensures that the index k+2,k+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+1).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+2).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+3).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+1))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+2))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+3))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k-1).ge.1) then !ensures that the index k-1 is not out of bounds
                 if (chrbdy2(i,j,k-1).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            end if !closes condition on (abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))


           end if !closes condition on bdy_extrap_order.eq.3





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

c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for setting the mask chrbdy to a value different from ex at grid points that are at large rho=sqrt(x**2+y**2+z**2), are not excised points, have coordinates in a range of fixed values for all resolutions (in order to be able to show convergence of boundary extrapolated quantities), and such that their neighbours used for extrapolation are also not excised.
c----------------------------------------------------------------------

        subroutine nexttobdypoints_fixedpts(
     &                  chrbdy,
     &                  numbdypoints,
     &                  bdy_extrap_order,
     &                  ind_distance_fixedpts,
     &                  currentres_ratio_Lhighres_Llowres,
     &                  half_steps_from_bdy_ext,
     &                  half_steps_from_bdy_int,
     &                  num_fixed_coords,
     &                  fixed_coords,
     &                  x,y,z,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz),chrbdy2(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer m

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1,rhop1,chip1,xip1
        real*8 maxxyzp1
        integer numbdypoints
        integer half_steps_from_bdy_ext
        integer half_steps_from_bdy_int
        integer bdy_extrap_order

        real*8 currentres_ratio_Lhighres_Llowres
        integer num_fixed_coords
        integer ind_distance_fixedpts
        real*8 fixed_coords(num_fixed_coords)

        integer isxp1fixcoord
        integer isyp1fixcoord
        integer iszp1fixcoord

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------


        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

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

!select only points that have radius smaller than (1-5*currentres_ratio_Lhighres_Llowres*dx/2).
!for the lowest degree resolution, currentres_ratio_Lhighres_Llowres=1. For each resolution, points between rhobdy=1 and 1-dx/2 are excised, points between rhobdy=1 and 1-3*dx/2 use forward/backward stencils (so when we cannot expect convergence at these points because the stencils used are different for different resolutions), points between rhobdy=1 and 1-5*dx/2 use points that are set by forward/backward stencils, so we cannot expect convergence. If we want to check convergence at the grid points used for extrapolation at the boundary, we need to pick points that have rho<1-5*dx/2. However, we see that we get good convergence only if we use points with rho<1-9*dx/2. This should be understood in the future. The condition rhop1>currentres_ratio_Lhighres_Llowres ensures that the first suitable point for extrapolation is not too far from the AdS boundary.

            if ((rhop1.lt.
     &          (1-half_steps_from_bdy_ext
     &             *currentres_ratio_Lhighres_Llowres*dx/2))
     &          .and.(rhop1.gt.
     &                (1-half_steps_from_bdy_int
     &                 *currentres_ratio_Lhighres_Llowres*dx/2))
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

           if (chrbdy(i,j,k).ne.ex) then

!select only points whose coordinates are in the group of fixed_coordinates for all resolutions
             isxp1fixcoord=0
             isyp1fixcoord=0
             iszp1fixcoord=0
             do m=1,num_fixed_coords
              if ((abs(xp1-fixed_coords(m)).lt.10.0d0**(-10)).and.
     &            (isxp1fixcoord.ne.1)) then
                   isxp1fixcoord=1
              end if
              if ((abs(yp1-fixed_coords(m)).lt.10.0d0**(-10)).and.
     &            (isyp1fixcoord.ne.1)) then
                   isyp1fixcoord=1
              end if
              if ((abs(zp1-fixed_coords(m)).lt.10.0d0**(-10)).and.
     &            (iszp1fixcoord.ne.1)) then
                   iszp1fixcoord=1
              end if
             end do

             if ((isxp1fixcoord.ne.1).or.
     &           (isyp1fixcoord.ne.1).or.
     &           (iszp1fixcoord.ne.1)) then
                   chrbdy(i,j,k)=ex
                   chrbdy2(i,j,k)=ex
             end if


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

           if (chrbdy(i,j,k).ne.ex) then


! eliminate troublesome points
!removing points with abs(x)=abs(y)=abs(z): even if the grid function quasiset_tracell is converging at these points, the extrapolated boundary quantity quasisettrace, obtained using these points, is not converging. It would be interesting to understand why. For now, we just eliminate these points from the set of points used for extrapolation
            if (
     &          (abs(abs(xp1)-abs(yp1)).lt.10.0d0**(-10)).and.
     &          (abs(abs(xp1)-abs(zp1)).lt.10.0d0**(-10))
     &         ) then
                  chrbdy(i,j,k)=ex
            end if



! eliminate troublesome points
!y=z=0 (i.e. the x axis, where chi=0 or 1) are troublesome points where the xi coordinate is not defined  
!We will fill and impose regularity at these points in post-processing in Mathematica.
            if (
     &          (abs(yp1).lt.10.0d0**(-10)).and.
     &          (abs(zp1).lt.10.0d0**(-10)) 
     &         ) then
             chrbdy(i,j,k)=ex
            end if


!NOTICE: for example, in the case of extrapolation along x, x>0, if the closest point to the AdS boundary that we use is (i,j,k), the second fixed (for all resolutions) point that we want to use is (i-ind_distance_fixedpts,j,k)
!If we use derivatives to define near boundary quantities, we will only define them at points between is and ie (js and je, ks and ke). Therefore, for extrapolation, we can only select near boundary points whose neighbors used for extrapolation in the direction of the bulk along the axes (i.e. the direction of extrapolation) are within that range
!We also define near boundary quantities at points where y0 and z0 are not both 0. So we need to make sure that we select points such that neighbouring points used for extrapolation don't have such values of y0,z0. Notice, we've already imposed that y(j)=!=0 and z(k)=!=0 for the outermost point. Therefore, when extrapolation is along y or z, we need to impose the condition that the (y,z) coordinates of the other points used are different from 0.
!The condition (chrbdy2(i+ind_distance_fixedpts,j,k).ne.ex) makes sure that (i,j,k) is the outmost point satisfying the conditions of the previous for-loop, which sets chrbdy2 as well as chrbdy. In other words, if there's an outer point w.r.t. (i,j,k) that satisfies those conditions, then we don't want to use (i,j,k) for extrapolation, but we will use that other point.

           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))
           if (bdy_extrap_order.eq.1) then
            if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
             if (xp1.gt.0) then
              if ((i-ind_distance_fixedpts).lt.is) then !it ensures that the index i-ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i-ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i+ind_distance_fixedpts).le.Nx) then !ensures that the index i+ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i+ind_distance_fixedpts,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((i+ind_distance_fixedpts).gt.ie) then !it ensures that the index i+ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i+ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i-ind_distance_fixedpts).ge.1) then !ensures that the index i-ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i-ind_distance_fixedpts,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
             if (yp1.gt.0) then
              if ((j-ind_distance_fixedpts).lt.js) then !it ensures that the index j-ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation  at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j+ind_distance_fixedpts).le.Ny) then !ensures that the index j+ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i,j+ind_distance_fixedpts,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((j+ind_distance_fixedpts).gt.je) then !it ensures that the index j+ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j-ind_distance_fixedpts).ge.1) then !ensures that the index j-ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i,j-ind_distance_fixedpts,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            else !i.e. when maxxyzp1.eq.abs(zp1)
             if (zp1.gt.0) then
              if ((k-ind_distance_fixedpts).lt.ks) then !it ensures that the index k-ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k+ind_distance_fixedpts).le.Nz) then !ensures that the index k+ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i,j,k+ind_distance_fixedpts).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((k+ind_distance_fixedpts).gt.ke) then !it ensures that the index k+ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k-ind_distance_fixedpts).ge.1) then !ensures that the index k-ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i,j,k-ind_distance_fixedpts).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            end if !closes condition on (abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))

           end if !closes condition on bdy_extrap_order.eq.1

           

           if (bdy_extrap_order.eq.2) then
            if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
             if (xp1.gt.0) then
              if ((i-2*ind_distance_fixedpts).lt.is) then !it ensures that the index i-2*ind_distance_fixedpts,i-ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i-ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i-2*ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i+ind_distance_fixedpts).le.Nx) then !ensures that the index i+ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i+ind_distance_fixedpts,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((i+2*ind_distance_fixedpts).gt.ie) then !it ensures that the index i+2*ind_distance_fixedpts,i+ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i+ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i+2*ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i-ind_distance_fixedpts).ge.1) then !ensures that the index i-ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i-ind_distance_fixedpts,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
             if (yp1.gt.0) then
              if ((j-2*ind_distance_fixedpts).lt.js) then !it ensures that the index j-2*ind_distance_fixedpts,j-ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-2*ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-2*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j+ind_distance_fixedpts).le.Ny) then !ensures that the index j+ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i,j+ind_distance_fixedpts,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((j+2*ind_distance_fixedpts).gt.je) then !it ensures that the index j+2*ind_distance_fixedpts,j+ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+2*ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+2*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j-ind_distance_fixedpts).ge.1) then !ensures that the index j-ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i,j-ind_distance_fixedpts,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            else !i.e. when maxxyzp1.eq.abs(zp1)
             if (zp1.gt.0) then
              if ((k-2*ind_distance_fixedpts).lt.ks) then !it ensures that the index k-2*ind_distance_fixedpts,k-ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-2*ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-2*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k+ind_distance_fixedpts).le.Nz) then !ensures that the index k+ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i,j,k+ind_distance_fixedpts).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((k+2*ind_distance_fixedpts).gt.ke) then !it ensures that the index k+2*ind_distance_fixedpts,k+ind_distance_fixedpts is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+2*ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+2*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k-ind_distance_fixedpts).ge.1) then !ensures that the index k-ind_distance_fixedpts is not out of bounds
                 if (chrbdy2(i,j,k-ind_distance_fixedpts).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            end if !closes condition on (abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))

           end if !closes condition on bdy_extrap_order.eq.2



          if (bdy_extrap_order.eq.3) then
            if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
             if (xp1.gt.0) then
              if ((i-3*ind_distance_fixedpts).lt.is) then !it ensures that the index i-3,i-2,i-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i-ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i-2*ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i-3*ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex     
              else if ((i+ind_distance_fixedpts).le.Nx) then !ensures that the index i+1 is not out of bounds
                 if (chrbdy2(i+ind_distance_fixedpts,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((i+3*ind_distance_fixedpts).gt.ie) then !it ensures that the index i+2,i+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i+ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i+2*ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i+3*ind_distance_fixedpts,j,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if ((i-ind_distance_fixedpts).ge.1) then !ensures that the index i-1 is not out of bounds
                 if (chrbdy2(i-ind_distance_fixedpts,j,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if


            else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
             if (yp1.gt.0) then
              if ((j-3*ind_distance_fixedpts).lt.js) then !it ensures that the index j-2,j-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-2*ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j-3*ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-2*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j-3*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j+ind_distance_fixedpts).le.Ny) then !ensures that the index j+1 is not out of bounds
                 if (chrbdy2(i,j+ind_distance_fixedpts,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((j+3*ind_distance_fixedpts).gt.je) then !it ensures that the index j+2,j+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+2*ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j+3*ind_distance_fixedpts,k).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+2*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(z(k)).lt.10.0d0**(-10)).and.
     &             (abs(y(j+3*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((j-ind_distance_fixedpts).ge.1) then !ensures that the index j-1 is not out of bounds
                 if (chrbdy2(i,j-ind_distance_fixedpts,k).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            else !i.e. when maxxyzp1.eq.abs(zp1)
             if (zp1.gt.0) then
              if ((k-3*ind_distance_fixedpts).lt.ks) then !it ensures that the index k-2,k-1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-2*ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k-3*ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-2*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k-3*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k+ind_distance_fixedpts).le.Nz) then !ensures that the index k+1 is not out of bounds
                 if (chrbdy2(i,j,k+ind_distance_fixedpts).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             else
              if ((k+3*ind_distance_fixedpts).gt.ke) then !it ensures that the index k+2,k+1 is not out of bounds and that near-boundary quantities are defined at other points used for extrapolation
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+2*ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (chr(i,j,k+3*ind_distance_fixedpts).eq.ex) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+2*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if (
     &             (abs(y(j)).lt.10.0d0**(-10)).and.
     &             (abs(z(k+3*ind_distance_fixedpts))
     &                 .lt.10.0d0**(-10))
     &            ) then
                   chrbdy(i,j,k)=ex
              else if ((k-ind_distance_fixedpts).ge.1) then !ensures that the index k-1 is not out of bounds
                 if (chrbdy2(i,j,k-ind_distance_fixedpts).ne.ex) then
                   chrbdy(i,j,k)=ex
                 end if
              end if
             end if

            end if !closes condition on (abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))


           end if !closes condition on bdy_extrap_order.eq.3



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
c routine for computing the Cartesian coordinates of each boundary point where we want to extrapolate the value of functions AND the coordinates the outermost point used for extrapolation at each boundary point
c-------------------------------------------------------------------------------------------------------------------------

        subroutine xyz_extrap_outermost(
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
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

            if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

              xpout(lind)=xp1
              ypout(lind)=yp1
              zpout(lind)=zp1

              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
               if (xp1.gt.0) then
                  xpbdy(lind)=sqrt(1-yp1**2-zp1**2)
                  ypbdy(lind)=yp1
                  zpbdy(lind)=zp1
              else
                  xpbdy(lind)=-sqrt(1-yp1**2-zp1**2)
                  ypbdy(lind)=yp1
                  zpbdy(lind)=zp1
              end if
             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  ypbdy(lind)=sqrt(1-xp1**2-zp1**2)
                  xpbdy(lind)=xp1
                  zpbdy(lind)=zp1
              else
                  ypbdy(lind)=-sqrt(1-xp1**2-zp1**2)
                  xpbdy(lind)=xp1
                  zpbdy(lind)=zp1
              end if
             else
                 if (zp1.gt.0) then
                  zpbdy(lind)=sqrt(1-yp1**2-xp1**2)
                  ypbdy(lind)=yp1
                  xpbdy(lind)=xp1
              else
                  zpbdy(lind)=-sqrt(1-yp1**2-xp1**2)
                  ypbdy(lind)=yp1
                  xpbdy(lind)=xp1
               end if
              end if

!               xpbdy(lind)=xp1   !TEST
!               ypbdy(lind)=yp1   !TEST
!               zpbdy(lind)=zp1   !TEST

            end if
          end do
         end do
        end do

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the leading order coefficient of the scalar field phi1 in the near bdy expansion
c in powers of q=1-rho about  q=0.
c----------------------------------------------------------------------

        subroutine calc_leadordcoeff_phi1(leadordcoeff_phi1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 leadordcoeff_phi1(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        logical no_derivatives
        data no_derivatives/.false./

        real*8 dphi1_drho

        real*8 test1(Nx,Ny,Nz)
        real*8 dtest1_drho
!----------------------------------------------------------------------

        if (no_derivatives) then
         do i=1,Nx
          do j=1,Ny
           do k=1,Nz
             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0

            if ((chr(i,j,k).ne.ex)
     &          .and.(
     &                (abs(y0).ge.10.0d0**(-10)).or.
     &                (abs(z0).ge.10.0d0**(-10))
     &               )
     &         ) then
               leadordcoeff_phi1(i,j,k)=phi1_n(i,j,k)/q
            else
               leadordcoeff_phi1(i,j,k)=0
            end if
           end do
          end do
         end do

        else !i.e. if .not.no_derivatives
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

         do i=is,ie
          do j=js,je
           do k=ks,ke
             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0

            if ((chr(i,j,k).ne.ex)
     &          .and.(
!the xi coordinate is not defined at y0=z0=0
     &                (abs(y0).ge.10.0d0**(-10)).or.
     &                (abs(z0).ge.10.0d0**(-10))
     &               )
     &         ) then

               call df_drho(phi1_n,dphi1_drho,
     &          x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
               leadordcoeff_phi1(i,j,k)=
     &                  -dphi1_drho
            else
               leadordcoeff_phi1(i,j,k)=0
            end if

           end do
          end do
         end do

        end if

        return
        end
c-----------------------------------------------------------------------------------


c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the scalar field phi of the boundary CFT
c----------------------------------------------------------------------

        subroutine extrap_bdyphi_freepts(bdyphi,
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
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 leadordcoeff_phi1(Nx,Ny,Nz)
        real*8 leadordcoeff_phi1_p1
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
        real*8 maxxyzp1

        real*8 firstord_extrap
        real*8 secondord_extrap
        real*8 thirdord_extrap
        real*8 fourthord_extrap
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
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)

           leadordcoeff_phi1_p1=leadordcoeff_phi1(i,j,k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

           if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

                   xex=xpbdy(lind)
                   yex=ypbdy(lind)
                   zex=zpbdy(lind)
                   rhoex=sqrt(xex**2+yex**2+zex**2)
                   chiex=(1/PI)*acos(xex/rhoex)
                if (zex.lt.0) then
                   xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
                else
                   xiex=(1/(2*PI))*atan2(zex,yex)
                end if

             if (bdy_extrap_order.eq.1) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then
                  xp2=x(i-1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i-1,j,k)
              else
                  xp2=x(i+1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i+1,j,k)
               end if
                  bdyphi(lind)=
     &                 firstord_extrap(leadordcoeff_phi1_p1
     &                  ,leadordcoeff_phi1_p2,xp1,xp2,xex)
             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j-1,k)
              else
                  yp2=y(j+1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j+1,k)
              end if
                  bdyphi(lind)=
     &                 firstord_extrap(leadordcoeff_phi1_p1
     &                  ,leadordcoeff_phi1_p2,yp1,yp2,yex)
             else
                 if (zp1.gt.0) then
                  zp2=z(k-1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j,k-1)
                 else
                  zp2=z(k+1)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j,k+1)
               end if
                  bdyphi(lind)=
     &                 firstord_extrap(leadordcoeff_phi1_p1
     &                  ,leadordcoeff_phi1_p2,zp1,zp2,zex)
              end if


!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.1



             if (bdy_extrap_order.eq.2) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then
                  xp2=x(i-1)
                  xp3=x(i-2)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i-1,j,k)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i-2,j,k)
              else
                  xp2=x(i+1)
                  xp3=x(i+2)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i+1,j,k)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i+2,j,k)
               end if
                  bdyphi(lind)=
     &                 secondord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   xp1,xp2,xp3,xex)
             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-1)
                  yp3=y(j-2)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j-1,k)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i,j-2,k)
              else
                  yp2=y(j+1)
                  yp3=y(j+2)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j+1,k)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i,j+2,k)
              end if
                  bdyphi(lind)=
     &                 secondord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   yp1,yp2,yp3,yex)
             else
                 if (zp1.gt.0) then
                  zp2=z(k-1)
                  zp3=z(k-2)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j,k-1)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i,j,k-2)
                 else
                  zp2=z(k+1)
                  zp3=z(k+2)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j,k+1)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i,j,k+2)
               end if
                  bdyphi(lind)=
     &                 secondord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   zp1,zp2,zp3,zex)
              end if


!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.2




             if (bdy_extrap_order.eq.3) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then
                  xp2=x(i-1)
                  xp3=x(i-2)
                  xp4=x(i-3)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i-1,j,k)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i-2,j,k)
                  leadordcoeff_phi1_p4=leadordcoeff_phi1(i-3,j,k)
              else
                  xp2=x(i+1)
                  xp3=x(i+2)
                  xp4=x(i+3)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i+1,j,k)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i+2,j,k)
                  leadordcoeff_phi1_p4=leadordcoeff_phi1(i+3,j,k)
               end if
                  bdyphi(lind)=
     &                 thirdord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   leadordcoeff_phi1_p4,
     &                   xp1,xp2,xp3,xp4,xex)
             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-1)
                  yp3=y(j-2)
                  yp4=y(j-3)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j-1,k)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i,j-2,k)
                  leadordcoeff_phi1_p4=leadordcoeff_phi1(i,j-3,k)
              else
                  yp2=y(j+1)
                  yp3=y(j+2)
                  yp4=y(j+3)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j+1,k)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i,j+2,k)
                  leadordcoeff_phi1_p4=leadordcoeff_phi1(i,j+3,k)
              end if
                  bdyphi(lind)=
     &                 thirdord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   leadordcoeff_phi1_p4,
     &                   yp1,yp2,yp3,yp4,yex)
             else
                 if (zp1.gt.0) then
                  zp2=z(k-1)
                  zp3=z(k-2)
                  zp4=z(k-3)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j,k-1)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i,j,k-2)
                  leadordcoeff_phi1_p4=leadordcoeff_phi1(i,j,k-3)
                 else
                  zp2=z(k+1)
                  zp3=z(k+2)
                  zp4=z(k+3)
                  leadordcoeff_phi1_p2=leadordcoeff_phi1(i,j,k+1)
                  leadordcoeff_phi1_p3=leadordcoeff_phi1(i,j,k+2)
                  leadordcoeff_phi1_p4=leadordcoeff_phi1(i,j,k+3)
               end if
                  bdyphi(lind)=
     &                 thirdord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   leadordcoeff_phi1_p4,
     &                   zp1,zp2,zp3,zp4,zex)
              end if


!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.3






           end if !closes condition on chrbdy(i,j,k).ne.ex




          end do
         end do
        end do


        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the scalar field phi of the boundary CFT, extrapolation with fixed points for all resolutions
c----------------------------------------------------------------------

        subroutine extrap_bdyphi_fixedpts(bdyphi,
     &                  leadordcoeff_phi1,
     &                  xpbdy,ypbdy,zpbdy,
     &                  chrbdy,numbdypoints,
     &                  bdy_extrap_order,
     &                  ind_distance_fixedpts,
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
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 leadordcoeff_phi1(Nx,Ny,Nz)
        real*8 leadordcoeff_phi1_p1
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
        real*8 maxxyzp1

        real*8 firstord_extrap
        real*8 secondord_extrap
        real*8 thirdord_extrap
        real*8 fourthord_extrap

        real*8 ratio_Lhighres_Llowres
        integer resolution_degree
        integer max_resolution_degree
        integer reduction_factor
        integer ind_distance_fixedpts
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
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)

           leadordcoeff_phi1_p1=leadordcoeff_phi1(i,j,k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

           if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

                   xex=xpbdy(lind)
                   yex=ypbdy(lind)
                   zex=zpbdy(lind)
                   rhoex=sqrt(xex**2+yex**2+zex**2)
                   chiex=(1/PI)*acos(xex/rhoex)
                if (zex.lt.0) then
                   xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
                else
                   xiex=(1/(2*PI))*atan2(zex,yex)
                end if

             if (bdy_extrap_order.eq.1) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then
                  xp2=x(i-ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i-ind_distance_fixedpts,j,k)
              else
                  xp2=x(i+ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i+ind_distance_fixedpts,j,k)
               end if
                  bdyphi(lind)=
     &                 firstord_extrap(leadordcoeff_phi1_p1
     &                  ,leadordcoeff_phi1_p2,xp1,xp2,xex)
             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j-ind_distance_fixedpts,k)
              else
                  yp2=y(j+ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j+ind_distance_fixedpts,k)
              end if
                  bdyphi(lind)=
     &                 firstord_extrap(leadordcoeff_phi1_p1
     &                  ,leadordcoeff_phi1_p2,yp1,yp2,yex)
             else
                 if (zp1.gt.0) then
                  zp2=z(k-ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j,k-ind_distance_fixedpts)
                 else
                  zp2=z(k+ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j,k+ind_distance_fixedpts)
               end if
                  bdyphi(lind)=
     &                 firstord_extrap(leadordcoeff_phi1_p1
     &                  ,leadordcoeff_phi1_p2,zp1,zp2,zex)
              end if


!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.1

             if (bdy_extrap_order.eq.2) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then
                  xp2=x(i-ind_distance_fixedpts)
                  xp3=x(i-2*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i-ind_distance_fixedpts,j,k)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i-2*ind_distance_fixedpts,j,k)
              else
                  xp2=x(i+ind_distance_fixedpts)
                  xp3=x(i+2*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i+ind_distance_fixedpts,j,k)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i+2*ind_distance_fixedpts,j,k)
               end if
                  bdyphi(lind)=
     &                 secondord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   xp1,xp2,xp3,xex)
             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-ind_distance_fixedpts)
                  yp3=y(j-2*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j-ind_distance_fixedpts,k)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i,j-2*ind_distance_fixedpts,k)
              else
                  yp2=y(j+ind_distance_fixedpts)
                  yp3=y(j+2*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j+ind_distance_fixedpts,k)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i,j+2*ind_distance_fixedpts,k)
              end if
                  bdyphi(lind)=
     &                 secondord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   yp1,yp2,yp3,yex)
             else
                 if (zp1.gt.0) then
                  zp2=z(k-ind_distance_fixedpts)
                  zp3=z(k-2*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j,k-ind_distance_fixedpts)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i,j,k-2*ind_distance_fixedpts)
                 else
                  zp2=z(k+ind_distance_fixedpts)
                  zp3=z(k+2*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j,k+ind_distance_fixedpts)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i,j,k+2*ind_distance_fixedpts)
               end if
                  bdyphi(lind)=
     &                 secondord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   zp1,zp2,zp3,zex)
              end if


!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.2



            if (bdy_extrap_order.eq.3) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then
                  xp2=x(i-1*ind_distance_fixedpts)
                  xp3=x(i-2*ind_distance_fixedpts)
                  xp4=x(i-3*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i-1*ind_distance_fixedpts,j,k)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i-2*ind_distance_fixedpts,j,k)
                  leadordcoeff_phi1_p4=
     &                leadordcoeff_phi1(i-3*ind_distance_fixedpts,j,k)
              else
                  xp2=x(i+1*ind_distance_fixedpts)
                  xp3=x(i+2*ind_distance_fixedpts)
                  xp4=x(i+3*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i+1*ind_distance_fixedpts,j,k)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i+2*ind_distance_fixedpts,j,k)
                  leadordcoeff_phi1_p4=
     &                leadordcoeff_phi1(i+3*ind_distance_fixedpts,j,k)
               end if
                  bdyphi(lind)=
     &                 thirdord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   leadordcoeff_phi1_p4,
     &                   xp1,xp2,xp3,xp4,xex)
             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-1*ind_distance_fixedpts)
                  yp3=y(j-2*ind_distance_fixedpts)
                  yp4=y(j-3*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j-1*ind_distance_fixedpts,k)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i,j-2*ind_distance_fixedpts,k)
                  leadordcoeff_phi1_p4=
     &                leadordcoeff_phi1(i,j-3*ind_distance_fixedpts,k)
              else
                  yp2=y(j+1*ind_distance_fixedpts)
                  yp3=y(j+2*ind_distance_fixedpts)
                  yp4=y(j+3*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j+1*ind_distance_fixedpts,k)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i,j+2*ind_distance_fixedpts,k)
                  leadordcoeff_phi1_p4=
     &                leadordcoeff_phi1(i,j+3*ind_distance_fixedpts,k)
              end if
                  bdyphi(lind)=
     &                 thirdord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   leadordcoeff_phi1_p4,
     &                   yp1,yp2,yp3,yp4,yex)
             else
                 if (zp1.gt.0) then
                  zp2=z(k-1*ind_distance_fixedpts)
                  zp3=z(k-2*ind_distance_fixedpts)
                  zp4=z(k-3*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j,k-1*ind_distance_fixedpts)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i,j,k-2*ind_distance_fixedpts)
                  leadordcoeff_phi1_p4=
     &                leadordcoeff_phi1(i,j,k-3*ind_distance_fixedpts)
                 else
                  zp2=z(k+1*ind_distance_fixedpts)
                  zp3=z(k+2*ind_distance_fixedpts)
                  zp4=z(k+3*ind_distance_fixedpts)
                  leadordcoeff_phi1_p2=
     &                leadordcoeff_phi1(i,j,k+1*ind_distance_fixedpts)
                  leadordcoeff_phi1_p3=
     &                leadordcoeff_phi1(i,j,k+2*ind_distance_fixedpts)
                  leadordcoeff_phi1_p4=
     &                leadordcoeff_phi1(i,j,k+3*ind_distance_fixedpts)
               end if
                  bdyphi(lind)=
     &                 thirdord_extrap(
     &                   leadordcoeff_phi1_p1,
     &                   leadordcoeff_phi1_p2,
     &                   leadordcoeff_phi1_p3,
     &                   leadordcoeff_phi1_p4,
     &                   zp1,zp2,zp3,zp4,zex)
              end if


!                  bdyphi(lind)=leadordcoeff_phi1_p1  !TEST
!              write(*,*) "lind-1,xp1,yp1,zp1",lind-1,xp1,yp1,zp1
!             write(*,*) "bdyphi(lind)=",bdyphi(lind)

            end if !closes condition on bdy_extrap_order.eq.3





           end if !closes condition on chrbdy(i,j,k).ne.ex
          end do
         end do
        end do


        return
        end
c--------------------------------------------------------------------------------------


        


c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D
c using a 1-rho expansion about rho=1 at points near the boundary. 
c The tensor components are given in quasi-spherical polar coordinates.
c The calculation is done as follows:
c 1. computing the deviations from pure AdS (these are not the evolution variables gb's if we evolve
c around a background metric that is not AdS)
c 2. using the formulae for the boundary stress-tensor obtained for a perturbation of pure AdS.
c----------------------------------------------------------------------

        subroutine calc_quasiset_ll(
     &                  quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
     &                  quasiset_chichi_ll,quasiset_chixi_ll,
     &                  quasiset_xixi_ll,
     &                  quasiset_tracell,
     &                  quasiset_massdensityll,
     &                  quasiset_angmomdensityxll,
     &                  quasiset_angmomdensityyll,
     &                  quasiset_angmomdensityzll,
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
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width,
     &                  ief_bh_r0,a_rot,kerrads_background)

!----------------------------------------------------------------------


        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
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
        real*8 phi1_np1(Nx,Ny,Nz),phi1_n(Nx,Ny,Nz),phi1_nm1(Nx,Ny,Nz)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 quasiset_tt_ll(Nx,Ny,Nz),quasiset_tchi_ll(Nx,Ny,Nz)
        real*8 quasiset_txi_ll(Nx,Ny,Nz),quasiset_chichi_ll(Nx,Ny,Nz)
        real*8 quasiset_chixi_ll(Nx,Ny,Nz),quasiset_xixi_ll(Nx,Ny,Nz)
        real*8 quasiset_tracell(Nx,Ny,Nz)
        real*8 quasiset_massdensityll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityxll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityyll(Nx,Ny,Nz)
        real*8 quasiset_angmomdensityzll(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)



        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t, gb_tt_x, gb_tt_y
        real*8 gb_tt_z
        real*8 gb_tt_tt,gb_tt_tx,gb_tt_ty
        real*8 gb_tt_tz
        real*8 gb_tt_xx,gb_tt_xy
        real*8 gb_tt_xz
        real*8 gb_tt_yy
        real*8 gb_tt_yz
        real*8 gb_tt_zz


        real*8 gb_tx_t, gb_tx_x, gb_tx_y
        real*8 gb_tx_z
        real*8 gb_tx_tt,gb_tx_tx,gb_tx_ty
        real*8 gb_tx_tz
        real*8 gb_tx_xx,gb_tx_xy
        real*8 gb_tx_xz
        real*8 gb_tx_yy
        real*8 gb_tx_yz
        real*8 gb_tx_zz

        real*8 gb_ty_t, gb_ty_x, gb_ty_y
        real*8 gb_ty_z
        real*8 gb_ty_tt,gb_ty_tx,gb_ty_ty
        real*8 gb_ty_tz
        real*8 gb_ty_xx,gb_ty_xy
        real*8 gb_ty_xz
        real*8 gb_ty_yy
        real*8 gb_ty_yz
        real*8 gb_ty_zz

        real*8 gb_tz_t, gb_tz_x, gb_tz_y
        real*8 gb_tz_z
        real*8 gb_tz_tt,gb_tz_tx,gb_tz_ty
        real*8 gb_tz_tz
        real*8 gb_tz_xx,gb_tz_xy
        real*8 gb_tz_xz
        real*8 gb_tz_yy
        real*8 gb_tz_yz
        real*8 gb_tz_zz


        real*8 gb_xx_t, gb_xx_x, gb_xx_y
        real*8 gb_xx_z
        real*8 gb_xx_tt,gb_xx_tx,gb_xx_ty
        real*8 gb_xx_tz
        real*8 gb_xx_xx,gb_xx_xy
        real*8 gb_xx_xz
        real*8 gb_xx_yy
        real*8 gb_xx_yz
        real*8 gb_xx_zz

        real*8 gb_xy_t, gb_xy_x, gb_xy_y
        real*8 gb_xy_z
        real*8 gb_xy_tt,gb_xy_tx,gb_xy_ty
        real*8 gb_xy_tz
        real*8 gb_xy_xx,gb_xy_xy
        real*8 gb_xy_xz
        real*8 gb_xy_yy
        real*8 gb_xy_yz
        real*8 gb_xy_zz

        real*8 gb_xz_t, gb_xz_x, gb_xz_y
        real*8 gb_xz_z
        real*8 gb_xz_tt,gb_xz_tx,gb_xz_ty
        real*8 gb_xz_tz
        real*8 gb_xz_xx,gb_xz_xy
        real*8 gb_xz_xz
        real*8 gb_xz_yy
        real*8 gb_xz_yz
        real*8 gb_xz_zz

        real*8 gb_yy_t, gb_yy_x, gb_yy_y
        real*8 gb_yy_z
        real*8 gb_yy_tt,gb_yy_tx,gb_yy_ty
        real*8 gb_yy_tz
        real*8 gb_yy_xx,gb_yy_xy
        real*8 gb_yy_xz
        real*8 gb_yy_yy
        real*8 gb_yy_yz
        real*8 gb_yy_zz

        real*8 gb_yz_t, gb_yz_x, gb_yz_y
        real*8 gb_yz_z
        real*8 gb_yz_tt,gb_yz_tx,gb_yz_ty
        real*8 gb_yz_tz
        real*8 gb_yz_xx,gb_yz_xy
        real*8 gb_yz_xz
        real*8 gb_yz_yy
        real*8 gb_yz_yz
        real*8 gb_yz_zz

        real*8 gb_zz_t, gb_zz_x, gb_zz_y
        real*8 gb_zz_z
        real*8 gb_zz_tt,gb_zz_tx,gb_zz_ty
        real*8 gb_zz_tz
        real*8 gb_zz_xx,gb_zz_xy
        real*8 gb_zz_xz
        real*8 gb_zz_yy
        real*8 gb_zz_yz
        real*8 gb_zz_zz


        real*8 dtdt,dtdrho,dtdchi,dtdxi
        real*8 dxdt,dxdrho,dxdchi,dxdxi
        real*8 dydt,dydrho,dydchi,dydxi
        real*8 dzdt,dzdrho,dzdchi,dzdxi

        real*8 g0_tt_ads0,g0_xx_ads0
        real*8 g0_tx_ads0,g0_ty_ads0,g0_tz_ads0
        real*8 g0_xy_ads0,g0_yy_ads0,g0_zz_ads0
        real*8 g0_xz_ads0,g0_yz_ads0

        real*8  ief_bh_r0,a_rot,M0,M0_min
        integer kerrads_background

        logical calc_der,calc_adv_quant
        data calc_der/.false./
        data calc_adv_quant/.false./

        real*8 gads_ll(4,4),gads_uu(4,4)
        real*8 gads_ll_x(4,4,4),gads_uu_x(4,4,4),gads_ll_xx(4,4,4,4)
        real*8 gads_ll_qssph(4,4),gads_uu_qssph(4,4)
        real*8 gads_ll_qssph_x(4,4,4),gads_uu_qssph_x(4,4,4)
        real*8 gads_ll_qssph_xx(4,4,4,4)
        real*8 gammaads_ull(4,4,4)
        real*8 boxadsx_u(4)
        real*8 Hads_l(4)
        real*8 phi1ads, phi1ads_x(4),phi1ads_xx(4,4)
        real*8 gammaads_ull_x(4,4,4,4)
        real*8 riemannads_ulll(4,4,4,4)
        real*8 ricciads_ll(4,4),ricciads_lu(4,4),ricciads
        real*8 einsteinads_ll(4,4),setads_ll(4,4)

        real*8 dxcar_dxqssph(4,4)
        real*8 hcar_n(4,4),g0car_n(4,4)
        real*8 hqssph_n(4,4),gbqssph_n(4,4),g0qssph_n(4,4)
        real*8 gamma0qssph_ll(4,4),gamma0qssph_uu(4,4)
        real*8 gamma0qssphbdy_ll(3,3),gamma0qssphbdy_uu(3,3)
        real*8 detgamma3
        real*8 detgamma0qssphbdy

        logical no_derivatives
        data no_derivatives/.false./
 
        real*8 gbqssph_tt_n(Nx,Ny,Nz),gbqssph_trho_n(Nx,Ny,Nz)
        real*8 gbqssph_tchi_n(Nx,Ny,Nz),gbqssph_txi_n(Nx,Ny,Nz)
        real*8 gbqssph_rhorho_n(Nx,Ny,Nz),gbqssph_rhochi_n(Nx,Ny,Nz)
        real*8 gbqssph_rhoxi_n(Nx,Ny,Nz)
        real*8 gbqssph_chichi_n(Nx,Ny,Nz),gbqssph_chixi_n(Nx,Ny,Nz)
        real*8 gbqssph_xixi_n(Nx,Ny,Nz)
        real*8 dgbqssph_tt_drho_n,dgbqssph_trho_drho_n
        real*8 dgbqssph_tchi_drho_n,dgbqssph_txi_drho_n
        real*8 dgbqssph_rhorho_drho_n
        real*8 dgbqssph_rhochi_drho_n
        real*8 dgbqssph_rhoxi_drho_n
        real*8 dgbqssph_chichi_drho_n
        real*8 dgbqssph_chixi_drho_n
        real*8 dgbqssph_xixi_drho_n

        real*8 gbqssph_tt_n_x,gbqssph_tt_n_y
        real*8 gb_tt_n_x,gb_tt_n_y
        real*8 gb_xx_n_x,gb_xx_n_y

        real*8 dergb_tt_x_n,dergbqssph_tt_x_n
        real*8 dergb_tt_y_n,dergbqssph_tt_y_n
        real*8 dergb_tt_z_n,dergbqssph_tt_z_n
        real*8 dergb_xx_x_n
        real*8 dergb_xx_y_n
        real*8 dergb_xx_z_n
        real*8 dergb_tt_x_np1,dergb_tt_y_np1,dergb_tt_z_np1

!----------------------------------------------------------------------

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

! define quantities necessary to define the quasi local stress energy tensor over the whole grid
       do i=1,Nx
        do j=1,Ny
         do k=1,Nz

             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0
             chi0=(1/PI)*acos(x0/rho0)
             if (z0.lt.0) then
                xi0=(1/(2*PI))*(atan2(z0,y0)+2*PI)
             else
                xi0=(1/(2*PI))*atan2(z0,y0)
             end if
!calculate regularized metric components in quasi-spherical coordinates in terms of regularized metric components in Cartesian coordinates
! we use the following coordinate transformation (notice that the angles are rescaled w.r.t. the usual quasi-spherical coordinates): x=rho*cos(PI*chi),y=rho*sin(PI*chi)*cos(2*PI*xi),z=rho*sin(PI*chi)*sin(2*PI*xi)

          if ( (chr(i,j,k).ne.ex)
     &          .and.(
!the xi coordinate is not defined at y0=z0=0
     &                (abs(y0).ge.10.0d0**(-10)).or.
     &                (abs(z0).ge.10.0d0**(-10))
     &               )
     &         ) then

!transformation matrix

             dtdt=1
             dtdrho=0
             dtdchi=0
             dtdxi=0
             dxdt=0
             dxdrho=cos(PI*chi0)
             dxdchi=-PI*rho0*sin(PI*chi0)
             dxdxi=0
             dydt=0
             dydrho=sin(PI*chi0)*cos(2*PI*xi0)
             dydchi=PI*rho0*cos(PI*chi0)*cos(2*PI*xi0)
             dydxi=-2*PI*rho0*sin(PI*chi0)*sin(2*PI*xi0)
             dzdt=0
             dzdrho=sin(PI*chi0)*sin(2*PI*xi0)
             dzdchi=PI*rho0*cos(PI*chi0)*sin(2*PI*xi0)
             dzdxi=2*PI*rho0*sin(PI*chi0)*cos(2*PI*xi0)

             dxcar_dxqssph(1,1)=dtdt
             dxcar_dxqssph(1,2)=dtdrho
             dxcar_dxqssph(1,3)=dtdchi
             dxcar_dxqssph(1,4)=dtdxi
             dxcar_dxqssph(2,1)=dxdt
             dxcar_dxqssph(2,2)=dxdrho
             dxcar_dxqssph(2,3)=dxdchi
             dxcar_dxqssph(2,4)=dxdxi
             dxcar_dxqssph(3,1)=dydt
             dxcar_dxqssph(3,2)=dydrho
             dxcar_dxqssph(3,3)=dydchi
             dxcar_dxqssph(3,4)=dydxi
             dxcar_dxqssph(4,1)=dzdt
             dxcar_dxqssph(4,2)=dzdrho
             dxcar_dxqssph(4,3)=dzdchi
             dxcar_dxqssph(4,4)=dzdxi

        !metric components of pure AdS in Cartesian coordinates
        g0_tt_ads0 =-(4*rho0**2+L**2*(-1+rho0**2)**2)
     &               /L**2/(-1+rho0**2)**2
        g0_tx_ads0 =0
        g0_ty_ads0 =0
        g0_tz_ads0 =0
        g0_xx_ads0 =(8*(-1+L**2)*(x0**2-y0**2-z0**2)
     &              +8*rho0**2+4*L**2*(1+rho0**4))
     &              /(-1+rho0**2)**2/(4*rho0**2+L**2*(-1+rho0**2)**2)
        g0_xy_ads0 =(16 *(-1 + L**2) *x0* y0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
        g0_xz_ads0 =(16 *(-1 + L**2) *x0* z0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
        g0_yy_ads0 =(4*(4*(x0**2+z0**2)+L**2*(x0**4+(1+y0**2)**2
     &              +2*(-1+y0**2)*z0**2+z0**4
     &              +2*x0**2*(-1+y0**2+z0**2))))
     &              /(L**2*(-1+rho0**2)**4+4*(-1+rho0**2)**2*(rho0**2))
        g0_yz_ads0 =(16 *(-1 + L**2) *y0* z0)
     &              /((-1 + rho0**2)**2
     &               *(4 *rho0**2 +L**2 *(-1 +rho0**2)**2))
        g0_zz_ads0=(4*(4*(x0**2+y0**2)+L**2*((-1+x0**2+y0**2)**2
     &              +2*(1+x0**2+y0**2)*z0**2+z0**4)))
     &              /(L**2*(-1+rho0**2)**4
     &              +4*(-1+rho0**2)**2*(rho0**2))


    !compute background metric
    ! NOTE: even if the background metric is not pure AdS, we still denote it by gads_ll,Hads_l,etc.
        if ((kerrads_background.eq.0).or.
     &      ((kerrads_background.eq.1).and.
     &        (ief_bh_r0.lt.10.0d0**(-10)).and.
     &        (a_rot.lt.10.0d0**(-10))) ) then
            call ads_derivs_cartcoords(
     &                  gads_ll,gads_uu,gads_ll_x,
     &                  gads_uu_x,gads_ll_xx,
     &                  Hads_l,
     &                  gammaads_ull,
     &                  phi1ads,
     &                  phi1ads_x,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                  calc_der,calc_adv_quant)
        else if ((kerrads_background.eq.1).and.
     &           (a_rot.lt.10.0d0**(-10)) ) then !the background metric is Schw-AdS in Kerr-Schild coords
            call schwads_derivs_kerrschildcoords(
     &                  gads_ll,gads_uu,gads_ll_x,
     &                  gads_uu_x,gads_ll_xx,
     &                  Hads_l,
     &                  gammaads_ull,
     &                  phi1ads,
     &                  phi1ads_x,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                  ief_bh_r0,
     &                  calc_der,calc_adv_quant)
        else if ((kerrads_background.eq.1).and.
     &           (a_rot.gt.10.0d0**(-10)) )  then !the background metric is Kerr-AdS in Kerr-Schild coords
            call kerrads_derivs_kerrschildcoords(
     &                  gads_ll,gads_uu,gads_ll_x,
     &                  gads_uu_x,gads_ll_xx,
     &                  Hads_l,
     &                  gammaads_ull,
     &                  phi1ads,
     &                  phi1ads_x,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                  ief_bh_r0,a_rot,
     &                  calc_der,calc_adv_quant)
        end if



    !compute background metric and its derivatives 
    ! NOTE: even if the background metric is not pure AdS, we still denote by g0_ll_ads0.

      !deviation from PURE AdS in Cartesian coordinates
      !This is the deviation from pure AdS regardless of the background metric we are using in the evolution
        hcar_n(1,1)=(gads_ll(1,1)-g0_tt_ads0)+gb_tt_n(i,j,k)
        hcar_n(1,2)=(gads_ll(1,2)-g0_tx_ads0)+gb_tx_n(i,j,k)
        hcar_n(1,3)=(gads_ll(1,3)-g0_ty_ads0)+gb_ty_n(i,j,k)
        hcar_n(1,4)=(gads_ll(1,4)-g0_tz_ads0)+gb_tz_n(i,j,k)
        hcar_n(2,2)=(gads_ll(2,2)-g0_xx_ads0)+gb_xx_n(i,j,k)
        hcar_n(2,3)=(gads_ll(2,3)-g0_xy_ads0)+gb_xy_n(i,j,k)
        hcar_n(2,4)=(gads_ll(2,4)-g0_xz_ads0)+gb_xz_n(i,j,k)
        hcar_n(3,3)=(gads_ll(3,3)-g0_yy_ads0)+gb_yy_n(i,j,k)
        hcar_n(3,4)=(gads_ll(3,4)-g0_yz_ads0)+gb_yz_n(i,j,k)
        hcar_n(4,4)=(gads_ll(4,4)-g0_zz_ads0)+gb_zz_n(i,j,k)

       do a=1,3
          do b=a+1,4
            hcar_n(b,a)=hcar_n(a,b)
          end do
        end do

!       !full metric in Cartesian coordinates
!        g0car_n(1,1)=g0_tt_ads0+hcar_n(1,1)
!        g0car_n(1,2)=g0_tx_ads0+hcar_n(1,2)
!        g0car_n(1,3)=g0_ty_ads0+hcar_n(1,3)
!        g0car_n(1,4)=g0_tz_ads0+hcar_n(1,4)
!        g0car_n(2,2)=g0_xx_ads0+hcar_n(2,2)
!        g0car_n(2,3)=g0_xy_ads0+hcar_n(2,3)
!        g0car_n(2,4)=g0_xz_ads0+hcar_n(2,4)
!        g0car_n(3,3)=g0_yy_ads0+hcar_n(3,3)
!        g0car_n(3,4)=g0_yz_ads0+hcar_n(3,4)
!        g0car_n(4,4)=g0_zz_ads0+hcar_n(4,4)
!
!       do a=1,3
!          do b=a+1,4
!            g0car_n(b,a)=g0car_n(a,b)
!          end do
!        end do


       !deviation from pure AdS in (rescaled) quasi-spherical coordinates
        do a=1,4
          do b=1,4
           hqssph_n(a,b)=0.0d0
           do c=1,4
            do d=1,4
             hqssph_n(a,b)=hqssph_n(a,b)
     &             +dxcar_dxqssph(c,a)*dxcar_dxqssph(d,b)*hcar_n(c,d)
            end do
           end do
          end do
        end do


        !regularised metric components in (rescaled) quasi-spherical coordinates for perturbation of pure AdS
         gbqssph_n(1,1)=hqssph_n(1,1)
         gbqssph_n(1,2)=hqssph_n(1,2)
         gbqssph_n(1,3)=hqssph_n(1,3)
         gbqssph_n(1,4)=hqssph_n(1,4)
         gbqssph_n(2,2)=hqssph_n(2,2)
         gbqssph_n(2,3)=hqssph_n(2,3)
         gbqssph_n(2,4)=hqssph_n(2,4)
         gbqssph_n(3,3)=hqssph_n(3,3)
         gbqssph_n(3,4)=hqssph_n(3,4)
         gbqssph_n(4,4)=hqssph_n(4,4)

         gbqssph_tt_n(i,j,k)    =gbqssph_n(1,1)
         gbqssph_trho_n(i,j,k)  =gbqssph_n(1,2)
         gbqssph_tchi_n(i,j,k)  =gbqssph_n(1,3)
         gbqssph_txi_n(i,j,k)   =gbqssph_n(1,4)
         gbqssph_rhorho_n(i,j,k)=gbqssph_n(2,2)
         gbqssph_rhochi_n(i,j,k)=gbqssph_n(2,3)
         gbqssph_rhoxi_n(i,j,k) =gbqssph_n(2,4)
         gbqssph_chichi_n(i,j,k)=gbqssph_n(3,3)
         gbqssph_chixi_n(i,j,k) =gbqssph_n(3,4)
         gbqssph_xixi_n(i,j,k)  =gbqssph_n(4,4)

! calculate the AdS-subtracted bulk gravity quasilocal stress-energy tensor,
! identified with the bdy CFT stress-energy tensor one-point function
!these are the coefficients of the lowest order terms (i.e. those contributing to the AdS mass) in the expansion of the non-zero components of the quasi-local stress-energy tensor

!initialise
              quasiset_tt_ll(i,j,k)=0
              quasiset_tchi_ll(i,j,k)=0
              quasiset_txi_ll(i,j,k)=0
              quasiset_chichi_ll(i,j,k)=0
              quasiset_chixi_ll(i,j,k)=0
              quasiset_xixi_ll(i,j,k)=0
              quasiset_tracell(i,j,k)=0
              quasiset_massdensityll(i,j,k)=0
              quasiset_angmomdensityxll(i,j,k)=0
              quasiset_angmomdensityyll(i,j,k)=0
              quasiset_angmomdensityzll(i,j,k)=0

            if (no_derivatives) then
!             if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
!
               quasiset_tt_ll(i,j,k)= 
     &              (12*(gbqssph_n(3,3)/q)
     &       + 8*PI**2*(gbqssph_n(2,2)/q)
     &       +  (3*(gbqssph_n(4,4)/q)/(sin(PI*chi0))**2)
     &       )/(64*PI**3)


               quasiset_tchi_ll(i,j,k)  = 
     &              (3*(gbqssph_n(1,3)/q))/(16*PI)


               quasiset_txi_ll(i,j,k)   = 
     &               (3*(gbqssph_n(1,4)/q))/(16*PI)

               quasiset_chichi_ll(i,j,k)=(3.0d0/16.0d0)*PI
     &                                   *(gbqssph_n(1,1)/q)
     &                                  -(1.0d0/8.0d0)*PI
     &                                   *(gbqssph_n(2,2)/q)
     &                    -(3*(gbqssph_n(4,4)/q)/(sin(PI*chi0))**2)
     &                                   /(64*PI)


               quasiset_chixi_ll(i,j,k) = 
     &              (3*(gbqssph_n(3,4)/q))/(16*PI)


               quasiset_xixi_ll(i,j,k)  =((sin(PI*chi0))**2*(-3
     &                                   *(gbqssph_n(3,3)/q)
     &                                  +PI**2*(3*(gbqssph_n(1,1)/q)
     &                                  -2*(gbqssph_n(2,2)/q)))
     &                                  )/(4*PI)

      !trace of quasi local stress-tensor in terms of regularised metric components (this should be the same as the one above, within numerical error)
               quasiset_tracell(i,j,k)=
     &           (3/(32*PI**3))
     &           *(4*PI**2*(gbqssph_n(1,1)/q)
     &             -4*(gbqssph_n(3,3)/q)
     &             -4*PI**2*(gbqssph_n(2,2)/q)
     &             -(gbqssph_n(4,4)/q)/((sin(PI*chi0))**2) )

      !multiply by the square root of the determinant on the boundary sphere, 2*PI^2 * sin(PI*chi0), 
      ! and integrate in chi0 from 0 to 1 and xi0 from 0 to 1 to obtain the total mass in AdS
               quasiset_massdensityll(i,j,k)=
     &                          (
     &                           12*(gbqssph_n(3,3)/q)
     &                          +8*PI**2*(gbqssph_n(2,2)/q)
     &                          +3*(gbqssph_n(4,4)/q)
     &                           /((sin(PI*chi0))**2)
     &                          )
     &                          /(64*PI**3)
               quasiset_angmomdensityxll(i,j,k)=
     &                          -3*(gbqssph_n(1,4)/q)
     &                         /(32*PI**2)
               quasiset_angmomdensityyll(i,j,k)=
     &                         3*(
     &                          2*sin(2*PI*xi0)*(gbqssph_n(1,3)/q)
     &                         +cos(2*PI*xi0)*
     &                          (cos(PI*chi0)/sin(PI*chi0))*
     &                             (gbqssph_n(1,4)/q)
     &                         )
     &                         /(32*PI**2)
               quasiset_angmomdensityzll(i,j,k)=
     &                         3*(
     &                          -2*cos(2*PI*xi0)*(gbqssph_n(1,3)/q)
     &                         +(cos(PI*chi0)/sin(PI*chi0))*
     &                           sin(2*PI*xi0)*(gbqssph_n(1,4)/q)
     &                         )
     &                         /(32*PI**2)


            end if

          else  !excised points or points with y0=z0=0

              gbqssph_tt_n(i,j,k)    =0
              gbqssph_trho_n(i,j,k)  =0
              gbqssph_tchi_n(i,j,k)  =0
              gbqssph_txi_n(i,j,k)   =0
              gbqssph_rhorho_n(i,j,k)=0
              gbqssph_rhochi_n(i,j,k)=0
              gbqssph_rhoxi_n(i,j,k) =0
              gbqssph_chichi_n(i,j,k)=0
              gbqssph_chixi_n(i,j,k) =0
              gbqssph_xixi_n(i,j,k)  =0


              quasiset_tt_ll(i,j,k)=0
              quasiset_tchi_ll(i,j,k)=0
              quasiset_txi_ll(i,j,k)=0
              quasiset_chichi_ll(i,j,k)=0
              quasiset_chixi_ll(i,j,k)=0
              quasiset_xixi_ll(i,j,k)=0
              quasiset_tracell(i,j,k)=0
              quasiset_massdensityll(i,j,k)=0

          end if


         end do
        end do
       end do


      if (.not.no_derivatives) then !we need to restric the range of definition to points owned only by the current process (not any other process) if we want to use derivatives
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1


!!!ghost_width not needed as long as we don't use derivatives
        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        do i=is,ie
         do j=js,je
          do k=ks,ke

             x0=x(i)
             y0=y(j)
             z0=z(k)
             rho0=sqrt(x0**2+y0**2+z0**2)
             q=1-rho0
             chi0=(1/PI)*acos(x0/rho0)
             if (z0.lt.0) then
                xi0=(1/(2*PI))*(atan2(z0,y0)+2*PI)
             else
                xi0=(1/(2*PI))*atan2(z0,y0)
             end if


!                ! calculate the AdS-subtracted bulk gravity quasilocal stress-energy tensor,
!                ! identified with the bdy CFT stress-energy tensor one-point function
!                !these are the coefficients of the lowest order terms (i.e. those contributing to the AdS mass) in the expansion of the non-zero components of the quasi-local stress-energy tensor
            if ( (chr(i,j,k).ne.ex)
     &          .and.(
!the xi coordinate is not defined at y0=z0=0
     &                (abs(y0).ge.10.0d0**(-10)).or.
     &                (abs(z0).ge.10.0d0**(-10))
     &               )
     &         ) then

              call df_drho(gbqssph_tt_n,    dgbqssph_tt_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_trho_n,  dgbqssph_trho_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_tchi_n,  dgbqssph_tchi_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_txi_n,   dgbqssph_txi_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_rhorho_n,dgbqssph_rhorho_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_rhochi_n,dgbqssph_rhochi_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_rhoxi_n, dgbqssph_rhoxi_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_chichi_n,dgbqssph_chichi_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_chixi_n, dgbqssph_chixi_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)
              call df_drho(gbqssph_xixi_n,  dgbqssph_xixi_drho_n,
     &              x,y,z,i,j,k,chr,ex,Nx,Ny,Nz)

               quasiset_tt_ll(i,j,k)=
     &                         (12*(-dgbqssph_chichi_drho_n)
     &                         + 8*PI**2*(-dgbqssph_rhorho_drho_n)
     &                         +  (3*(-dgbqssph_xixi_drho_n)
     &                               /(sin(PI*chi0))**2)
     &                         )/(64*PI**3)


               quasiset_tchi_ll(i,j,k)  =
     &                        (3*(-dgbqssph_tchi_drho_n))
     &                                       /(16*PI)


               quasiset_txi_ll(i,j,k)   =
     &                        (3*(-dgbqssph_txi_drho_n))
     &                                       /(16*PI)


               quasiset_chichi_ll(i,j,k)=
     &                                   (3.0d0/16.0d0)*PI
     &                                   *(-dgbqssph_tt_drho_n)
     &                                  -(1.0d0/8.0d0)*PI
     &                                   *(-dgbqssph_rhorho_drho_n)
     &                           -(3*(-dgbqssph_xixi_drho_n)
     &                               /(sin(PI*chi0))**2)
     &                                   /(64*PI)


               quasiset_chixi_ll(i,j,k) =
     &                                   (3*(-dgbqssph_chixi_drho_n))
     &                                      /(16*PI)



               quasiset_xixi_ll(i,j,k)=
     &                                  ((sin(PI*chi0))**2*(-3
     &                                   *(-dgbqssph_chichi_drho_n)
     &                                  +PI**2*(3*(-dgbqssph_tt_drho_n)
     &                                  -2*(-dgbqssph_rhorho_drho_n)))
     &                                  )/(4*PI)


!        write(*,*) "COMPUTING TRACE DIRECTLY FROM DERIVATIVES 
!    &               OF METRIC COMPONENTS"
!
             quasiset_tracell(i,j,k)=
     &        (3/(32*PI**3))
     &           *(4*PI**2*(-dgbqssph_tt_drho_n)
     &             -4*(-dgbqssph_chichi_drho_n)
     &             -4*PI**2*(-dgbqssph_rhorho_drho_n)
     &             -(-dgbqssph_xixi_drho_n)/((sin(PI*chi0))**2) )

      !multiply by the square root of the determinant on the boundary sphere, 2*PI^2 * sin(PI*chi0), 
      ! and integrate in chi0 from 0 to 1 and xi0 from 0 to 1 to obtain the total mass in AdS
               quasiset_massdensityll(i,j,k)=
     &                        (
     &                         12*(-dgbqssph_chichi_drho_n)
     &                        +8*PI**2
     &                            *(-dgbqssph_rhorho_drho_n)
     &                        +3*(-dgbqssph_xixi_drho_n)
     &                         /((sin(PI*chi0))**2)
     &                        )
     &                        /(64*PI**3)

              quasiset_angmomdensityxll(i,j,k)=
     &                        -3*(-dgbqssph_txi_drho_n)
     &                       /(32*PI**2)
              quasiset_angmomdensityyll(i,j,k)=
     &                        3*(
     &                         2*sin(2*PI*xi0)*(-dgbqssph_tchi_drho_n)
     &                        +cos(2*PI*xi0)*
     &                          (cos(PI*chi0)/sin(PI*chi0))*
     &                             (-dgbqssph_txi_drho_n)
     &                        )
     &                        /(32*PI**2)
              quasiset_angmomdensityzll(i,j,k)=
     &                        3*(
     &                         -2*cos(2*PI*xi0)*(-dgbqssph_tchi_drho_n)
     &                        +(cos(PI*chi0)/sin(PI*chi0))*
     &                          sin(2*PI*xi0)*(-dgbqssph_txi_drho_n)
     &                        )
     &                        /(32*PI**2)

           else !excised points or points with y0=z0=0

              quasiset_tt_ll(i,j,k)=0
              quasiset_tchi_ll(i,j,k)=0
              quasiset_txi_ll(i,j,k)=0
              quasiset_chichi_ll(i,j,k)=0
              quasiset_chixi_ll(i,j,k)=0
              quasiset_xixi_ll(i,j,k)=0
              quasiset_tracell(i,j,k)=0
              quasiset_massdensityll(i,j,k)=0
              quasiset_angmomdensityxll(i,j,k)=0
              quasiset_angmomdensityyll(i,j,k)=0
              quasiset_angmomdensityzll(i,j,k)=0
            end if

          end do
         end do
        end do

       end if

       return
       end
c--------------------------------------------------------------------------------------



c-------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D_polar  
c using a 1-rho expansion about rho=1 AT the boundary through extrapolation.
c The tensor components are given in quasi-spherical polar coordinates.
c-------------------------------------------------------------------------------------

        subroutine extrap_quasiset_freepts(
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


        real*8 gamma0qssphbdy_uu_tt
        real*8 gamma0qssphbdy_uu_tchi
        real*8 gamma0qssphbdy_uu_txi
        real*8 gamma0qssphbdy_uu_chichi
        real*8 gamma0qssphbdy_uu_chixi
        real*8 gamma0qssphbdy_uu_xixi

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

        real*8 x0,y0,z0,rho0,q

        real*8 xp1,yp1,zp1,rhop1
        real*8 xp2,yp2,zp2
        real*8 xp3,yp3,zp3
        real*8 xp4,yp4,zp4
        real*8 xex,yex,zex,rhoex,chiex,xiex
        real*8 maxxyzp1

        real*8 firstord_extrap
        real*8 secondord_extrap
        real*8 thirdord_extrap
        real*8 fourthord_extrap

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
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)

           quasiset_tt_p1=quasiset_tt_ll(i,j,k)
           quasiset_tchi_p1=quasiset_tchi_ll(i,j,k)
           quasiset_txi_p1=quasiset_txi_ll(i,j,k)
           quasiset_chichi_p1=quasiset_chichi_ll(i,j,k)
           quasiset_chixi_p1=quasiset_chixi_ll(i,j,k)
           quasiset_xixi_p1=quasiset_xixi_ll(i,j,k)
           quasiset_trace_p1=quasiset_tracell(i,j,k)
           quasiset_massdensity_p1=quasiset_massdensityll(i,j,k)
           quasiset_angmomdensityx_p1=quasiset_angmomdensityxll(i,j,k)
           quasiset_angmomdensityy_p1=quasiset_angmomdensityyll(i,j,k)
           quasiset_angmomdensityz_p1=quasiset_angmomdensityzll(i,j,k)
           maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

           if (chrbdy(i,j,k).ne.ex) then

              lind=lind+1

                  xex=xpbdy(lind)
                  yex=ypbdy(lind)
                  zex=zpbdy(lind)
                  rhoex=sqrt(xex**2+yex**2+zex**2)
                  chiex=(1/PI)*acos(xex/rhoex)
                  if (zex.lt.0) then
                    xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
                  else
                    xiex=(1/(2*PI))*atan2(zex,yex)
                  end if

                  !inverse of conformal metric on AdS boundary (needed for trace) at extrapolated point
                  gamma0qssphbdy_uu_tt=-1
                  gamma0qssphbdy_uu_tchi=0
                  gamma0qssphbdy_uu_txi=0
                  gamma0qssphbdy_uu_chichi=1/(PI**2)
                  gamma0qssphbdy_uu_chixi=0
                  gamma0qssphbdy_uu_xixi=1/((sin(PI*chiex))**2)/4/PI**2

             if (bdy_extrap_order.eq.1) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then

                  xp2=x(i-1)
                  quasiset_tt_p2=quasiset_tt_ll(i-1,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i-1,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i-1,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i-1,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i-1,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i-1,j,k)
                  quasiset_trace_p2=quasiset_tracell(i-1,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i-1,j,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i-1,j,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i-1,j,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i-1,j,k)


              else
                  xp2=x(i+1)
                  quasiset_tt_p2=quasiset_tt_ll(i+1,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i+1,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i+1,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i+1,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i+1,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i+1,j,k)
                  quasiset_trace_p2=quasiset_tracell(i+1,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i+1,j,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i+1,j,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i+1,j,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i+1,j,k)


               end if

                 quasiset_tt(lind)=
     &                 firstord_extrap(quasiset_tt_p1
     &                  ,quasiset_tt_p2,xp1,xp2,xex)
                 quasiset_tchi(lind)=
     &                 firstord_extrap(quasiset_tchi_p1
     &                  ,quasiset_tchi_p2,xp1,xp2,xex)
                 quasiset_txi(lind)=
     &                 firstord_extrap(quasiset_txi_p1
     &                  ,quasiset_txi_p2,xp1,xp2,xex)
                 quasiset_chichi(lind)=
     &                 firstord_extrap(quasiset_chichi_p1
     &                  ,quasiset_chichi_p2,xp1,xp2,xex)
                 quasiset_chixi(lind)=
     &                 firstord_extrap(quasiset_chixi_p1
     &                  ,quasiset_chixi_p2,xp1,xp2,xex)
                 quasiset_xixi(lind)=
     &                 firstord_extrap(quasiset_xixi_p1
     &                  ,quasiset_xixi_p2,xp1,xp2,xex)
                 quasiset_trace(lind)=
     &                 firstord_extrap(quasiset_trace_p1
     &                  ,quasiset_trace_p2,xp1,xp2,xex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 firstord_extrap(quasiset_massdensity_p1
     &                  ,quasiset_massdensity_p2,xp1,xp2,xex)
                 quasiset_angmomdensityx(lind)=
     &                 firstord_extrap(quasiset_angmomdensityx_p1
     &                  ,quasiset_angmomdensityx_p2,xp1,xp2,xex)
                 quasiset_angmomdensityy(lind)=
     &                 firstord_extrap(quasiset_angmomdensityy_p1
     &                  ,quasiset_angmomdensityy_p2,xp1,xp2,xex)
                 quasiset_angmomdensityz(lind)=
     &                 firstord_extrap(quasiset_angmomdensityz_p1
     &                  ,quasiset_angmomdensityz_p2,xp1,xp2,xex)

             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j-1,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j-1,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j-1,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j-1,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j-1,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j-1,k)
                  quasiset_trace_p2=quasiset_tracell(i,j-1,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j-1,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j-1,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j-1,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j-1,k)


 
              else
                  yp2=y(j+1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j+1,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j+1,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j+1,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j+1,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j+1,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j+1,k)
                  quasiset_trace_p2=quasiset_tracell(i,j+1,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j+1,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j+1,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j+1,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j+1,k)


              end if             

                 quasiset_tt(lind)=
     &                 firstord_extrap(quasiset_tt_p1
     &                  ,quasiset_tt_p2,yp1,yp2,yex)
                 quasiset_tchi(lind)=
     &                 firstord_extrap(quasiset_tchi_p1
     &                  ,quasiset_tchi_p2,yp1,yp2,yex)
                 quasiset_txi(lind)=
     &                 firstord_extrap(quasiset_txi_p1
     &                  ,quasiset_txi_p2,yp1,yp2,yex)
                 quasiset_chichi(lind)=
     &                 firstord_extrap(quasiset_chichi_p1
     &                  ,quasiset_chichi_p2,yp1,yp2,yex)
                 quasiset_chixi(lind)=
     &                 firstord_extrap(quasiset_chixi_p1
     &                  ,quasiset_chixi_p2,yp1,yp2,yex)
                 quasiset_xixi(lind)=
     &                 firstord_extrap(quasiset_xixi_p1
     &                  ,quasiset_xixi_p2,yp1,yp2,yex)
                 quasiset_trace(lind)=
     &                 firstord_extrap(quasiset_trace_p1
     &                  ,quasiset_trace_p2,yp1,yp2,yex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 firstord_extrap(quasiset_massdensity_p1
     &                  ,quasiset_massdensity_p2,yp1,yp2,yex)
                 quasiset_angmomdensityx(lind)=
     &                 firstord_extrap(quasiset_angmomdensityx_p1
     &                  ,quasiset_angmomdensityx_p2,yp1,yp2,yex)
                 quasiset_angmomdensityy(lind)=
     &                 firstord_extrap(quasiset_angmomdensityy_p1
     &                  ,quasiset_angmomdensityy_p2,yp1,yp2,yex)
                 quasiset_angmomdensityz(lind)=
     &                 firstord_extrap(quasiset_angmomdensityz_p1
     &                  ,quasiset_angmomdensityz_p2,yp1,yp2,yex)
             else
                 if (zp1.gt.0) then
                  zp2=z(k-1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k-1)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k-1)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k-1)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k-1)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k-1)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k-1)
                  quasiset_trace_p2=quasiset_tracell(i,j,k-1)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k-1)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j,k-1)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j,k-1)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j,k-1)


              else
                  zp2=z(k+1)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k+1)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k+1)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k+1)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k+1)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k+1)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k+1)
                  quasiset_trace_p2=quasiset_tracell(i,j,k+1)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k+1)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j,k+1)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j,k+1)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j,k+1)


               end if

                 quasiset_tt(lind)=
     &                 firstord_extrap(quasiset_tt_p1
     &                  ,quasiset_tt_p2,zp1,zp2,zex)
                 quasiset_tchi(lind)=
     &                 firstord_extrap(quasiset_tchi_p1
     &                  ,quasiset_tchi_p2,zp1,zp2,zex)
                 quasiset_txi(lind)=
     &                 firstord_extrap(quasiset_txi_p1
     &                  ,quasiset_txi_p2,zp1,zp2,zex)
                 quasiset_chichi(lind)=
     &                 firstord_extrap(quasiset_chichi_p1
     &                  ,quasiset_chichi_p2,zp1,zp2,zex)
                 quasiset_chixi(lind)=
     &                 firstord_extrap(quasiset_chixi_p1
     &                  ,quasiset_chixi_p2,zp1,zp2,zex)
                 quasiset_xixi(lind)=
     &                 firstord_extrap(quasiset_xixi_p1
     &                  ,quasiset_xixi_p2,zp1,zp2,zex)
                 quasiset_trace(lind)=
     &                 firstord_extrap(quasiset_trace_p1
     &                  ,quasiset_trace_p2,zp1,zp2,zex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 firstord_extrap(quasiset_massdensity_p1
     &                  ,quasiset_massdensity_p2,zp1,zp2,zex)
                  quasiset_angmomdensityx(lind)=
     &                 firstord_extrap(quasiset_angmomdensityx_p1
     &                  ,quasiset_angmomdensityx_p2,zp1,zp2,zex)
                  quasiset_angmomdensityy(lind)=
     &                 firstord_extrap(quasiset_angmomdensityy_p1
     &                  ,quasiset_angmomdensityy_p2,zp1,zp2,zex)
                  quasiset_angmomdensityz(lind)=
     &                 firstord_extrap(quasiset_angmomdensityz_p1
     &                  ,quasiset_angmomdensityz_p2,zp1,zp2,zex)
              end if

!               quasiset_massdensity(lind)=sin(PI*chiex)*cos(2*PI*xiex)  !TEST

            end if !closes condition on bdy_extrap_order.eq.1




             if (bdy_extrap_order.eq.2) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then

                  xp2=x(i-1)
                  xp3=x(i-2)
                  quasiset_tt_p2=quasiset_tt_ll(i-1,j,k)
                  quasiset_tt_p3=quasiset_tt_ll(i-2,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i-1,j,k)
                  quasiset_tchi_p3=quasiset_tchi_ll(i-2,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i-1,j,k)
                  quasiset_txi_p3=quasiset_txi_ll(i-2,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i-1,j,k)
                  quasiset_chichi_p3=quasiset_chichi_ll(i-2,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i-1,j,k)
                  quasiset_chixi_p3=quasiset_chixi_ll(i-2,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i-1,j,k)
                  quasiset_xixi_p3=quasiset_xixi_ll(i-2,j,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i-1,j,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i-2,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i-1,j,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i-2,j,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i-1,j,k)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i-2,j,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i-1,j,k)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i-2,j,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i-1,j,k)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i-2,j,k)


              else
                  xp2=x(i+1)
                  xp3=x(i+2)
                  quasiset_tt_p2=quasiset_tt_ll(i+1,j,k)
                  quasiset_tt_p3=quasiset_tt_ll(i+2,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i+1,j,k)
                  quasiset_tchi_p3=quasiset_tchi_ll(i+2,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i+1,j,k)
                  quasiset_txi_p3=quasiset_txi_ll(i+2,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i+1,j,k)
                  quasiset_chichi_p3=quasiset_chichi_ll(i+2,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i+1,j,k)
                  quasiset_chixi_p3=quasiset_chixi_ll(i+2,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i+1,j,k)
                  quasiset_xixi_p3=quasiset_xixi_ll(i+2,j,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i+1,j,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i+2,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i+1,j,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i+2,j,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i+1,j,k)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i+2,j,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i+1,j,k)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i+2,j,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i+1,j,k)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i+2,j,k)


               end if

                 quasiset_tt(lind)=
     &                 secondord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_tchi(lind)=
     &                 secondord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_txi(lind)=
     &                 secondord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_chichi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_chixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_xixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_trace(lind)=
     &                 secondord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   xp1,xp2,xp3,xex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 secondord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_angmomdensityx(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_angmomdensityy(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_angmomdensityz(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   xp1,xp2,xp3,xex)


             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-1)
                  yp3=y(j-2)
                  quasiset_tt_p2=quasiset_tt_ll(i,j-1,k)
                  quasiset_tt_p3=quasiset_tt_ll(i,j-2,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j-1,k)
                  quasiset_tchi_p3=quasiset_tchi_ll(i,j-2,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j-1,k)
                  quasiset_txi_p3=quasiset_txi_ll(i,j-2,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j-1,k)
                  quasiset_chichi_p3=quasiset_chichi_ll(i,j-2,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j-1,k)
                  quasiset_chixi_p3=quasiset_chixi_ll(i,j-2,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j-1,k)
                  quasiset_xixi_p3=quasiset_xixi_ll(i,j-2,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j-1,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j-2,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j-1,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j-2,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j-1,k)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i,j-2,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j-1,k)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i,j-2,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j-1,k)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i,j-2,k)



              else
                  yp2=y(j+1)
                  yp3=y(j+2)
                  quasiset_tt_p2=quasiset_tt_ll(i,j+1,k)
                  quasiset_tt_p3=quasiset_tt_ll(i,j+2,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j+1,k)
                  quasiset_tchi_p3=quasiset_tchi_ll(i,j+2,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j+1,k)
                  quasiset_txi_p3=quasiset_txi_ll(i,j+2,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j+1,k)
                  quasiset_chichi_p3=quasiset_chichi_ll(i,j+2,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j+1,k)
                  quasiset_chixi_p3=quasiset_chixi_ll(i,j+2,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j+1,k)
                  quasiset_xixi_p3=quasiset_xixi_ll(i,j+2,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j+1,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j+2,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j+1,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j+2,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j+1,k)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i,j+2,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j+1,k)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i,j+2,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j+1,k)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i,j+2,k)


              end if
                 quasiset_tt(lind)=
     &                 secondord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_tchi(lind)=
     &                 secondord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_txi(lind)=
     &                 secondord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_chichi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_chixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_xixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_trace(lind)=
     &                 secondord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   yp1,yp2,yp3,yex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 secondord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_angmomdensityx(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_angmomdensityy(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_angmomdensityz(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   yp1,yp2,yp3,yex)
             else
                 if (zp1.gt.0) then
                  zp2=z(k-1)
                  zp3=z(k-2)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k-1)
                  quasiset_tt_p3=quasiset_tt_ll(i,j,k-2)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k-1)
                  quasiset_tchi_p3=quasiset_tchi_ll(i,j,k-2)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k-1)
                  quasiset_txi_p3=quasiset_txi_ll(i,j,k-2)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k-1)
                  quasiset_chichi_p3=quasiset_chichi_ll(i,j,k-2)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k-1)
                  quasiset_chixi_p3=quasiset_chixi_ll(i,j,k-2)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k-1)
                  quasiset_xixi_p3=quasiset_xixi_ll(i,j,k-2)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j,k-1)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j,k-2)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k-1)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j,k-2)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j,k-1)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i,j,k-2)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j,k-1)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i,j,k-2)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j,k-1)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i,j,k-2)

              else
                  zp2=z(k+1)
                  zp3=z(k+2)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k+1)
                  quasiset_tt_p3=quasiset_tt_ll(i,j,k+2)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k+1)
                  quasiset_tchi_p3=quasiset_tchi_ll(i,j,k+2)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k+1)
                  quasiset_txi_p3=quasiset_txi_ll(i,j,k+2)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k+1)
                  quasiset_chichi_p3=quasiset_chichi_ll(i,j,k+2)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k+1)
                  quasiset_chixi_p3=quasiset_chixi_ll(i,j,k+2)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k+1)
                  quasiset_xixi_p3=quasiset_xixi_ll(i,j,k+2)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j,k+1)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j,k+2)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k+1)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j,k+2)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j,k+1)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i,j,k+2)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j,k+1)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i,j,k+2)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j,k+1)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i,j,k+2)

               end if
                 quasiset_tt(lind)=
     &                 secondord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_tchi(lind)=
     &                 secondord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_txi(lind)=
     &                 secondord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_chichi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_chixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_xixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_trace(lind)=
     &                 secondord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   zp1,zp2,zp3,zex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 secondord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_angmomdensityx(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_angmomdensityy(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_angmomdensityz(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   zp1,zp2,zp3,zex)
              end if

!               quasiset_massdensity(lind)=sin(PI*chiex)*cos(2*PI*xiex)  !TEST

            end if !closes condition on bdy_extrap_order.eq.2




             if (bdy_extrap_order.eq.3) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then

                  xp2=x(i-1)
                  xp3=x(i-2)
                  xp4=x(i-3)
                  quasiset_tt_p2=quasiset_tt_ll(i-1,j,k)
                  quasiset_tt_p3=quasiset_tt_ll(i-2,j,k)
                  quasiset_tt_p4=quasiset_tt_ll(i-3,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i-1,j,k)
                  quasiset_tchi_p3=quasiset_tchi_ll(i-2,j,k)
                  quasiset_tchi_p4=quasiset_tchi_ll(i-3,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i-1,j,k)
                  quasiset_txi_p3=quasiset_txi_ll(i-2,j,k)
                  quasiset_txi_p4=quasiset_txi_ll(i-3,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i-1,j,k)
                  quasiset_chichi_p3=quasiset_chichi_ll(i-2,j,k)
                  quasiset_chichi_p4=quasiset_chichi_ll(i-3,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i-1,j,k)
                  quasiset_chixi_p3=quasiset_chixi_ll(i-2,j,k)
                  quasiset_chixi_p4=quasiset_chixi_ll(i-3,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i-1,j,k)
                  quasiset_xixi_p3=quasiset_xixi_ll(i-2,j,k)
                  quasiset_xixi_p4=quasiset_xixi_ll(i-3,j,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i-1,j,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i-2,j,k)
                  quasiset_trace_p4=
     &             quasiset_tracell(i-3,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i-1,j,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i-2,j,k)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i-3,j,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i-1,j,k)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i-2,j,k)
                  quasiset_angmomdensityx_p4=
     &             quasiset_angmomdensityxll(i-3,j,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i-1,j,k)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i-2,j,k)
                  quasiset_angmomdensityy_p4=
     &             quasiset_angmomdensityyll(i-3,j,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i-1,j,k)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i-2,j,k)
                  quasiset_angmomdensityz_p4=
     &             quasiset_angmomdensityzll(i-3,j,k)


              else
                  xp2=x(i+1)
                  xp3=x(i+2)
                  xp4=x(i+3)
                  quasiset_tt_p2=quasiset_tt_ll(i+1,j,k)
                  quasiset_tt_p3=quasiset_tt_ll(i+2,j,k)
                  quasiset_tt_p4=quasiset_tt_ll(i+3,j,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i+1,j,k)
                  quasiset_tchi_p3=quasiset_tchi_ll(i+2,j,k)
                  quasiset_tchi_p4=quasiset_tchi_ll(i+3,j,k)
                  quasiset_txi_p2=quasiset_txi_ll(i+1,j,k)
                  quasiset_txi_p3=quasiset_txi_ll(i+2,j,k)
                  quasiset_txi_p4=quasiset_txi_ll(i+3,j,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i+1,j,k)
                  quasiset_chichi_p3=quasiset_chichi_ll(i+2,j,k)
                  quasiset_chichi_p4=quasiset_chichi_ll(i+3,j,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i+1,j,k)
                  quasiset_chixi_p3=quasiset_chixi_ll(i+2,j,k)
                  quasiset_chixi_p4=quasiset_chixi_ll(i+3,j,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i+1,j,k)
                  quasiset_xixi_p3=quasiset_xixi_ll(i+2,j,k)
                  quasiset_xixi_p4=quasiset_xixi_ll(i+3,j,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i+1,j,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i+2,j,k)
                  quasiset_trace_p4=
     &             quasiset_tracell(i+3,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i+1,j,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i+2,j,k)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i+3,j,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i+1,j,k)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i+2,j,k)
                  quasiset_angmomdensityx_p4=
     &             quasiset_angmomdensityxll(i+3,j,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i+1,j,k)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i+2,j,k)
                  quasiset_angmomdensityy_p4=
     &             quasiset_angmomdensityyll(i+3,j,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i+1,j,k)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i+2,j,k)
                  quasiset_angmomdensityz_p4=
     &             quasiset_angmomdensityzll(i+3,j,k)


               end if

                 quasiset_tt(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   quasiset_tt_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_tchi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   quasiset_tchi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_txi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   quasiset_txi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_chichi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   quasiset_chichi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_chixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   quasiset_chixi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_xixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   quasiset_xixi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_trace(lind)=
     &                 thirdord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   quasiset_trace_p4,
     &                   xp1,xp2,xp3,xp4,xex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 thirdord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   quasiset_massdensity_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_angmomdensityx(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   quasiset_angmomdensityx_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_angmomdensityy(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   quasiset_angmomdensityy_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_angmomdensityz(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   quasiset_angmomdensityz_p4,
     &                   xp1,xp2,xp3,xp4,xex)


             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
                if (yp1.gt.0) then

                  yp2=y(j-1)
                  yp3=y(j-2)
                  yp4=y(j-3)
                  quasiset_tt_p2=quasiset_tt_ll(i,j-1,k)
                  quasiset_tt_p3=quasiset_tt_ll(i,j-2,k)
                  quasiset_tt_p4=quasiset_tt_ll(i,j-3,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j-1,k)
                  quasiset_tchi_p3=quasiset_tchi_ll(i,j-2,k)
                  quasiset_tchi_p4=quasiset_tchi_ll(i,j-3,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j-1,k)
                  quasiset_txi_p3=quasiset_txi_ll(i,j-2,k)
                  quasiset_txi_p4=quasiset_txi_ll(i,j-3,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j-1,k)
                  quasiset_chichi_p3=quasiset_chichi_ll(i,j-2,k)
                  quasiset_chichi_p4=quasiset_chichi_ll(i,j-3,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j-1,k)
                  quasiset_chixi_p3=quasiset_chixi_ll(i,j-2,k)
                  quasiset_chixi_p4=quasiset_chixi_ll(i,j-3,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j-1,k)
                  quasiset_xixi_p3=quasiset_xixi_ll(i,j-2,k)
                  quasiset_xixi_p4=quasiset_xixi_ll(i,j-3,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j-1,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j-2,k)
                  quasiset_trace_p4=
     &             quasiset_tracell(i,j-3,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j-1,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j-2,k)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i,j-3,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j-1,k)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i,j-2,k)
                  quasiset_angmomdensityx_p4=
     &             quasiset_angmomdensityxll(i,j-3,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j-1,k)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i,j-2,k)
                  quasiset_angmomdensityy_p4=
     &             quasiset_angmomdensityyll(i,j-3,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j-1,k)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i,j-2,k)
                  quasiset_angmomdensityz_p4=
     &             quasiset_angmomdensityzll(i,j-3,k)


              else
                  yp2=y(j+1)
                  yp3=y(j+2)
                  yp4=y(j+3)
                  quasiset_tt_p2=quasiset_tt_ll(i,j+1,k)
                  quasiset_tt_p3=quasiset_tt_ll(i,j+2,k)
                  quasiset_tt_p4=quasiset_tt_ll(i,j+3,k)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j+1,k)
                  quasiset_tchi_p3=quasiset_tchi_ll(i,j+2,k)
                  quasiset_tchi_p4=quasiset_tchi_ll(i,j+3,k)
                  quasiset_txi_p2=quasiset_txi_ll(i,j+1,k)
                  quasiset_txi_p3=quasiset_txi_ll(i,j+2,k)
                  quasiset_txi_p4=quasiset_txi_ll(i,j+3,k)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j+1,k)
                  quasiset_chichi_p3=quasiset_chichi_ll(i,j+2,k)
                  quasiset_chichi_p4=quasiset_chichi_ll(i,j+3,k)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j+1,k)
                  quasiset_chixi_p3=quasiset_chixi_ll(i,j+2,k)
                  quasiset_chixi_p4=quasiset_chixi_ll(i,j+3,k)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j+1,k)
                  quasiset_xixi_p3=quasiset_xixi_ll(i,j+2,k)
                  quasiset_xixi_p4=quasiset_xixi_ll(i,j+3,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j+1,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j+2,k)
                  quasiset_trace_p4=
     &             quasiset_tracell(i,j+3,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j+1,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j+2,k)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i,j+3,k)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j+1,k)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i,j+2,k)
                  quasiset_angmomdensityx_p4=
     &             quasiset_angmomdensityxll(i,j+3,k)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j+1,k)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i,j+2,k)
                  quasiset_angmomdensityy_p4=
     &             quasiset_angmomdensityyll(i,j+3,k)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j+1,k)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i,j+2,k)
                  quasiset_angmomdensityz_p4=
     &             quasiset_angmomdensityzll(i,j+3,k)



               end if

                 quasiset_tt(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   quasiset_tt_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_tchi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   quasiset_tchi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_txi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   quasiset_txi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_chichi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   quasiset_chichi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_chixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   quasiset_chixi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_xixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   quasiset_xixi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_trace(lind)=
     &                 thirdord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   quasiset_trace_p4,
     &                   yp1,yp2,yp3,yp4,yex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 thirdord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   quasiset_massdensity_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_angmomdensityx(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   quasiset_angmomdensityx_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_angmomdensityy(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   quasiset_angmomdensityy_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_angmomdensityz(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   quasiset_angmomdensityz_p4,
     &                   yp1,yp2,yp3,yp4,yex)
             else
                if (zp1.gt.0) then

                  zp2=z(k-1)
                  zp3=z(k-2)
                  zp4=z(k-3)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k-1)
                  quasiset_tt_p3=quasiset_tt_ll(i,j,k-2)
                  quasiset_tt_p4=quasiset_tt_ll(i,j,k-3)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k-1)
                  quasiset_tchi_p3=quasiset_tchi_ll(i,j,k-2)
                  quasiset_tchi_p4=quasiset_tchi_ll(i,j,k-3)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k-1)
                  quasiset_txi_p3=quasiset_txi_ll(i,j,k-2)
                  quasiset_txi_p4=quasiset_txi_ll(i,j,k-3)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k-1)
                  quasiset_chichi_p3=quasiset_chichi_ll(i,j,k-2)
                  quasiset_chichi_p4=quasiset_chichi_ll(i,j,k-3)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k-1)
                  quasiset_chixi_p3=quasiset_chixi_ll(i,j,k-2)
                  quasiset_chixi_p4=quasiset_chixi_ll(i,j,k-3)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k-1)
                  quasiset_xixi_p3=quasiset_xixi_ll(i,j,k-2)
                  quasiset_xixi_p4=quasiset_xixi_ll(i,j,k-3)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j,k-1)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j,k-2)
                  quasiset_trace_p4=
     &             quasiset_tracell(i,j,k-3)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k-1)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j,k-2)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i,j,k-3)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j,k-1)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i,j,k-2)
                  quasiset_angmomdensityx_p4=
     &             quasiset_angmomdensityxll(i,j,k-3)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j,k-1)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i,j,k-2)
                  quasiset_angmomdensityy_p4=
     &             quasiset_angmomdensityyll(i,j,k-3)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j,k-1)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i,j,k-2)
                  quasiset_angmomdensityz_p4=
     &             quasiset_angmomdensityzll(i,j,k-3)


              else
                  zp2=z(k+1)
                  zp3=z(k+2)
                  zp4=z(k+3)
                  quasiset_tt_p2=quasiset_tt_ll(i,j,k+1)
                  quasiset_tt_p3=quasiset_tt_ll(i,j,k+2)
                  quasiset_tt_p4=quasiset_tt_ll(i,j,k+3)
                  quasiset_tchi_p2=quasiset_tchi_ll(i,j,k+1)
                  quasiset_tchi_p3=quasiset_tchi_ll(i,j,k+2)
                  quasiset_tchi_p4=quasiset_tchi_ll(i,j,k+3)
                  quasiset_txi_p2=quasiset_txi_ll(i,j,k+1)
                  quasiset_txi_p3=quasiset_txi_ll(i,j,k+2)
                  quasiset_txi_p4=quasiset_txi_ll(i,j,k+3)
                  quasiset_chichi_p2=quasiset_chichi_ll(i,j,k+1)
                  quasiset_chichi_p3=quasiset_chichi_ll(i,j,k+2)
                  quasiset_chichi_p4=quasiset_chichi_ll(i,j,k+3)
                  quasiset_chixi_p2=quasiset_chixi_ll(i,j,k+1)
                  quasiset_chixi_p3=quasiset_chixi_ll(i,j,k+2)
                  quasiset_chixi_p4=quasiset_chixi_ll(i,j,k+3)
                  quasiset_xixi_p2=quasiset_xixi_ll(i,j,k+1)
                  quasiset_xixi_p3=quasiset_xixi_ll(i,j,k+2)
                  quasiset_xixi_p4=quasiset_xixi_ll(i,j,k+3)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j,k+1)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j,k+2)
                  quasiset_trace_p4=
     &             quasiset_tracell(i,j,k+3)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k+1)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j,k+2)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i,j,k+3)
                  quasiset_angmomdensityx_p2=
     &             quasiset_angmomdensityxll(i,j,k+1)
                  quasiset_angmomdensityx_p3=
     &             quasiset_angmomdensityxll(i,j,k+2)
                  quasiset_angmomdensityx_p4=
     &             quasiset_angmomdensityxll(i,j,k+3)
                  quasiset_angmomdensityy_p2=
     &             quasiset_angmomdensityyll(i,j,k+1)
                  quasiset_angmomdensityy_p3=
     &             quasiset_angmomdensityyll(i,j,k+2)
                  quasiset_angmomdensityy_p4=
     &             quasiset_angmomdensityyll(i,j,k+3)
                  quasiset_angmomdensityz_p2=
     &             quasiset_angmomdensityzll(i,j,k+1)
                  quasiset_angmomdensityz_p3=
     &             quasiset_angmomdensityzll(i,j,k+2)
                  quasiset_angmomdensityz_p4=
     &             quasiset_angmomdensityzll(i,j,k+3)


               end if

                 quasiset_tt(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   quasiset_tt_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_tchi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   quasiset_tchi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_txi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   quasiset_txi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_chichi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   quasiset_chichi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_chixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   quasiset_chixi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_xixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   quasiset_xixi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_trace(lind)=
     &                 thirdord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   quasiset_trace_p4,
     &                   zp1,zp2,zp3,zp4,zex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 thirdord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   quasiset_massdensity_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_angmomdensityx(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   quasiset_angmomdensityx_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_angmomdensityy(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   quasiset_angmomdensityy_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_angmomdensityz(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   quasiset_angmomdensityz_p4,
     &                   zp1,zp2,zp3,zp4,zex)
              end if

!               quasiset_massdensity(lind)=sin(PI*chiex)*cos(2*PI*xiex)  !TEST

            end if !closes condition on bdy_extrap_order.eq.3




           end if !closes condition on chrbdy(i,j,k).ne.ex
          end do
         end do
        end do

!        do lind=1,numbdypoints
!         write(*,*) "lind,xpbdy(lind),ypbdy(lind),zpbdy(lind)="
!     &              ,lind,xpbdy(lind),ypbdy(lind),zpbdy(lind)
!         write(*,*) "quasiset_trace(lind)=",quasiset_trace(lind)
!        end do

        return
        end
c--------------------------------------------------------------------------------------

c-------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the asymptotic quasilocal stress-energy of AdS4D_polar  
c using a 1-rho expansion about rho=1 AT the boundary through extrapolation.
c The tensor components are given in quasi-spherical polar coordinates.i
c Extrapolation at fixed points for all resolutions
c-------------------------------------------------------------------------------------

        subroutine extrap_quasiset_fixedpts(
     &          quasiset_tt,quasiset_tchi,quasiset_txi,
     &          quasiset_chichi,quasiset_chixi,
     &          quasiset_xixi,
     &          quasiset_trace,
     &          quasiset_massdensity,
     &          quasiset_angmomdensityx,
     &          quasiset_angmomdensityy,
     &          quasiset_angmomdensityz,
     &          quasiset_tt_ll,quasiset_tchi_ll,quasiset_txi_ll,
     &          quasiset_chichi_ll,quasiset_chixi_ll,
     &          quasiset_xixi_ll,
     &          quasiset_tracell,
     &          quasiset_massdensityll,
     &          quasiset_angmomdensityxll,
     &          quasiset_angmomdensityyll,
     &          quasiset_angmomdensityzll,
     &          xpbdy,ypbdy,zpbdy,
     &          chrbdy,numbdypoints,
     &          bdy_extrap_order,
     &          ind_distance_fixedpts,
     &          x,y,z,dt,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none
 
        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        integer numbdypoints
        integer bdy_extrap_order
        real*8 chrbdy(Nx,Ny,Nz)
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


        real*8 gamma0qssphbdy_uu_tt
        real*8 gamma0qssphbdy_uu_tchi
        real*8 gamma0qssphbdy_uu_txi
        real*8 gamma0qssphbdy_uu_chichi
        real*8 gamma0qssphbdy_uu_chixi
        real*8 gamma0qssphbdy_uu_xixi

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

!       real*8 chibdy(bdy_Nchi)
!       real*8 xibdy(bdy_Nxi)

       integer i,j,k,is,ie,js,je,ks,ke,lind,m,e,increase

       integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
       integer ix,jy,kz

       real*8 x0,y0,z0,rho0,q

       real*8 xp1,yp1,zp1,rhop1
       real*8 xp2,yp2,zp2
       real*8 xp3,yp3,zp3
       real*8 xp4,yp4,zp4
       real*8 xex,yex,zex,rhoex,chiex,xiex
       real*8 maxxyzp1

       real*8 firstord_extrap
       real*8 secondord_extrap
       real*8 thirdord_extrap
       real*8 fourthord_extrap

       real*8 dp1p2

       real*8 AdS_mass

       real*8 ratio_Lhighres_Llowres
       integer resolution_degree
       integer max_resolution_degree
       real*8 reduction_factor
       integer ind_distance_fixedpts

!----------------------------------------------------------------------

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

       lind=0
       do i=is,ie
        do j=js,je
         do k=ks,ke
          xp1=x(i)
          yp1=y(j)
          zp1=z(k)
          rhop1=sqrt(xp1**2+yp1**2+zp1**2)

          quasiset_tt_p1=quasiset_tt_ll(i,j,k)
          quasiset_tchi_p1=quasiset_tchi_ll(i,j,k)
          quasiset_txi_p1=quasiset_txi_ll(i,j,k)
          quasiset_chichi_p1=quasiset_chichi_ll(i,j,k)
          quasiset_chixi_p1=quasiset_chixi_ll(i,j,k)
          quasiset_xixi_p1=quasiset_xixi_ll(i,j,k)
          quasiset_trace_p1=quasiset_tracell(i,j,k)
          quasiset_massdensity_p1=quasiset_massdensityll(i,j,k)
          quasiset_angmomdensityx_p1=quasiset_angmomdensityxll(i,j,k)
          quasiset_angmomdensityy_p1=quasiset_angmomdensityyll(i,j,k)
          quasiset_angmomdensityz_p1=quasiset_angmomdensityzll(i,j,k)
          maxxyzp1=max(abs(xp1),abs(yp1),abs(zp1))

          if (chrbdy(i,j,k).ne.ex) then

             lind=lind+1

                 xex=xpbdy(lind)
                 yex=ypbdy(lind)
                 zex=zpbdy(lind)
                 rhoex=sqrt(xex**2+yex**2+zex**2)
                 chiex=(1/PI)*acos(xex/rhoex)
                 if (zex.lt.0) then
                   xiex=(1/(2*PI))*(atan2(zex,yex)+2*PI)
                 else
                   xiex=(1/(2*PI))*atan2(zex,yex)
                 end if

                 !inverse of conformal metric on AdS boundary (needed for trace) at extrapolated point
                 gamma0qssphbdy_uu_tt=-1
                 gamma0qssphbdy_uu_tchi=0
                 gamma0qssphbdy_uu_txi=0
                 gamma0qssphbdy_uu_chichi=1/(PI**2)
                 gamma0qssphbdy_uu_chixi=0
                 gamma0qssphbdy_uu_xixi=1/((sin(PI*chiex))**2)/4/PI**2

           if (bdy_extrap_order.eq.1) then
             if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
               if (xp1.gt.0) then

                 xp2=x(i-ind_distance_fixedpts)
                 quasiset_tt_p2=
     &               quasiset_tt_ll(i-ind_distance_fixedpts,j,k)
                 quasiset_tchi_p2=
     &               quasiset_tchi_ll(i-ind_distance_fixedpts,j,k)
                 quasiset_txi_p2=
     &               quasiset_txi_ll(i-ind_distance_fixedpts,j,k)
                 quasiset_chichi_p2=
     &               quasiset_chichi_ll(i-ind_distance_fixedpts,j,k)
                 quasiset_chixi_p2=
     &               quasiset_chixi_ll(i-ind_distance_fixedpts,j,k)
                 quasiset_xixi_p2=
     &               quasiset_xixi_ll(i-ind_distance_fixedpts,j,k)
                 quasiset_trace_p2=
     &               quasiset_tracell(i-ind_distance_fixedpts,j,k)
                 quasiset_massdensity_p2=
     &               quasiset_massdensityll(i-ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityx_p2=
     &      quasiset_angmomdensityxll(i-ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityy_p2=
     &      quasiset_angmomdensityyll(i-ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityz_p2=
     &      quasiset_angmomdensityzll(i-ind_distance_fixedpts,j,k)



             else
                 xp2=x(i+ind_distance_fixedpts)
                 quasiset_tt_p2=
     &               quasiset_tt_ll(i+ind_distance_fixedpts,j,k)
                 quasiset_tchi_p2=
     &               quasiset_tchi_ll(i+ind_distance_fixedpts,j,k)
                 quasiset_txi_p2=
     &               quasiset_txi_ll(i+ind_distance_fixedpts,j,k)
                 quasiset_chichi_p2=
     &               quasiset_chichi_ll(i+ind_distance_fixedpts,j,k)
                 quasiset_chixi_p2=
     &               quasiset_chixi_ll(i+ind_distance_fixedpts,j,k)
                 quasiset_xixi_p2=
     &               quasiset_xixi_ll(i+ind_distance_fixedpts,j,k)
                 quasiset_trace_p2=
     &               quasiset_tracell(i+ind_distance_fixedpts,j,k)
                 quasiset_massdensity_p2=
     &               quasiset_massdensityll(i+ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityx_p2=
     &      quasiset_angmomdensityxll(i+ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityy_p2=
     &      quasiset_angmomdensityyll(i+ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityz_p2=
     &      quasiset_angmomdensityzll(i+ind_distance_fixedpts,j,k)


              end if

                 quasiset_tt(lind)=
     &                 firstord_extrap(quasiset_tt_p1
     &                  ,quasiset_tt_p2,xp1,xp2,xex)
                 quasiset_tchi(lind)=
     &                 firstord_extrap(quasiset_tchi_p1
     &                  ,quasiset_tchi_p2,xp1,xp2,xex)
                 quasiset_txi(lind)=
     &                 firstord_extrap(quasiset_txi_p1
     &                  ,quasiset_txi_p2,xp1,xp2,xex)
                 quasiset_chichi(lind)=
     &                 firstord_extrap(quasiset_chichi_p1
     &                  ,quasiset_chichi_p2,xp1,xp2,xex)
                 quasiset_chixi(lind)=
     &                 firstord_extrap(quasiset_chixi_p1
     &                  ,quasiset_chixi_p2,xp1,xp2,xex)
                 quasiset_xixi(lind)=
     &                 firstord_extrap(quasiset_xixi_p1
     &                  ,quasiset_xixi_p2,xp1,xp2,xex)
                 quasiset_trace(lind)=
     &                 firstord_extrap(quasiset_trace_p1
     &                  ,quasiset_trace_p2,xp1,xp2,xex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 firstord_extrap(quasiset_massdensity_p1
     &                  ,quasiset_massdensity_p2,xp1,xp2,xex)
                 quasiset_angmomdensityx(lind)=
     &                 firstord_extrap(quasiset_angmomdensityx_p1
     &                  ,quasiset_angmomdensityx_p2,xp1,xp2,xex)
                 quasiset_angmomdensityy(lind)=
     &                 firstord_extrap(quasiset_angmomdensityy_p1
     &                  ,quasiset_angmomdensityy_p2,xp1,xp2,xex)
                 quasiset_angmomdensityz(lind)=
     &                 firstord_extrap(quasiset_angmomdensityz_p1
     &                  ,quasiset_angmomdensityz_p2,xp1,xp2,xex)
 
                  
            else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
             if (yp1.gt.0) then
                 yp2=y(j-ind_distance_fixedpts)
                 quasiset_tt_p2=
     &               quasiset_tt_ll(i,j-ind_distance_fixedpts,k)
                 quasiset_tchi_p2=
     &               quasiset_tchi_ll(i,j-ind_distance_fixedpts,k)
                 quasiset_txi_p2=
     &               quasiset_txi_ll(i,j-ind_distance_fixedpts,k)
                 quasiset_chichi_p2=
     &               quasiset_chichi_ll(i,j-ind_distance_fixedpts,k)
                 quasiset_chixi_p2=
     &               quasiset_chixi_ll(i,j-ind_distance_fixedpts,k)
                 quasiset_xixi_p2=
     &               quasiset_xixi_ll(i,j-ind_distance_fixedpts,k)
                 quasiset_trace_p2=
     &               quasiset_tracell(i,j-ind_distance_fixedpts,k)
                 quasiset_massdensity_p2=
     &               quasiset_massdensityll(i,j-ind_distance_fixedpts,k)
                 quasiset_angmomdensityx_p2=
     &      quasiset_angmomdensityxll(i,j-ind_distance_fixedpts,k)
                 quasiset_angmomdensityy_p2=
     &      quasiset_angmomdensityyll(i,j-ind_distance_fixedpts,k)
                 quasiset_angmomdensityz_p2=
     &      quasiset_angmomdensityzll(i,j-ind_distance_fixedpts,k)


             else
                 yp2=y(j+ind_distance_fixedpts)
                 quasiset_tt_p2=
     &               quasiset_tt_ll(i,j+ind_distance_fixedpts,k)
                 quasiset_tchi_p2=
     &               quasiset_tchi_ll(i,j+ind_distance_fixedpts,k)
                 quasiset_txi_p2=
     &               quasiset_txi_ll(i,j+ind_distance_fixedpts,k)
                 quasiset_chichi_p2=
     &               quasiset_chichi_ll(i,j+ind_distance_fixedpts,k)
                 quasiset_chixi_p2=
     &               quasiset_chixi_ll(i,j+ind_distance_fixedpts,k)
                 quasiset_xixi_p2=
     &               quasiset_xixi_ll(i,j+ind_distance_fixedpts,k)
                 quasiset_trace_p2=
     &               quasiset_tracell(i,j+ind_distance_fixedpts,k)
                 quasiset_massdensity_p2=
     &               quasiset_massdensityll(i,j+ind_distance_fixedpts,k)
                 quasiset_angmomdensityx_p2=
     &      quasiset_angmomdensityxll(i,j+ind_distance_fixedpts,k)
                 quasiset_angmomdensityy_p2=
     &      quasiset_angmomdensityyll(i,j+ind_distance_fixedpts,k)
                 quasiset_angmomdensityz_p2=
     &      quasiset_angmomdensityzll(i,j+ind_distance_fixedpts,k)




             end if

                 quasiset_tt(lind)=
     &                 firstord_extrap(quasiset_tt_p1
     &                  ,quasiset_tt_p2,yp1,yp2,yex)
                 quasiset_tchi(lind)=
     &                 firstord_extrap(quasiset_tchi_p1
     &                  ,quasiset_tchi_p2,yp1,yp2,yex)
                 quasiset_txi(lind)=
     &                 firstord_extrap(quasiset_txi_p1
     &                  ,quasiset_txi_p2,yp1,yp2,yex)
                 quasiset_chichi(lind)=
     &                 firstord_extrap(quasiset_chichi_p1
     &                  ,quasiset_chichi_p2,yp1,yp2,yex)
                 quasiset_chixi(lind)=
     &                 firstord_extrap(quasiset_chixi_p1
     &                  ,quasiset_chixi_p2,yp1,yp2,yex)
                 quasiset_xixi(lind)=
     &                 firstord_extrap(quasiset_xixi_p1
     &                  ,quasiset_xixi_p2,yp1,yp2,yex)
                 quasiset_trace(lind)=
     &                 firstord_extrap(quasiset_trace_p1
     &                  ,quasiset_trace_p2,yp1,yp2,yex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 firstord_extrap(quasiset_massdensity_p1
     &                  ,quasiset_massdensity_p2,yp1,yp2,yex)
                 quasiset_angmomdensityx(lind)=
     &                 firstord_extrap(quasiset_angmomdensityx_p1
     &                  ,quasiset_angmomdensityx_p2,yp1,yp2,yex)
                 quasiset_angmomdensityy(lind)=
     &                 firstord_extrap(quasiset_angmomdensityy_p1
     &                  ,quasiset_angmomdensityy_p2,yp1,yp2,yex)
                 quasiset_angmomdensityz(lind)=
     &                 firstord_extrap(quasiset_angmomdensityz_p1
     &                  ,quasiset_angmomdensityz_p2,yp1,yp2,yex)
            else
                if (zp1.gt.0) then
                 zp2=z(k-ind_distance_fixedpts)
                 quasiset_tt_p2=
     &               quasiset_tt_ll(i,j,k-ind_distance_fixedpts)
                 quasiset_tchi_p2=
     &               quasiset_tchi_ll(i,j,k-ind_distance_fixedpts)
                 quasiset_txi_p2=
     &               quasiset_txi_ll(i,j,k-ind_distance_fixedpts)
                 quasiset_chichi_p2=
     &               quasiset_chichi_ll(i,j,k-ind_distance_fixedpts)
                 quasiset_chixi_p2=
     &               quasiset_chixi_ll(i,j,k-ind_distance_fixedpts)
                 quasiset_xixi_p2=
     &               quasiset_xixi_ll(i,j,k-ind_distance_fixedpts)
                 quasiset_trace_p2=
     &               quasiset_tracell(i,j,k-ind_distance_fixedpts)
                 quasiset_massdensity_p2=
     &               quasiset_massdensityll(i,j,k-ind_distance_fixedpts)
                 quasiset_angmomdensityx_p2=
     &      quasiset_angmomdensityxll(i,j,k-ind_distance_fixedpts)
                 quasiset_angmomdensityy_p2=
     &      quasiset_angmomdensityyll(i,j,k-ind_distance_fixedpts)
                 quasiset_angmomdensityz_p2=
     &      quasiset_angmomdensityzll(i,j,k-ind_distance_fixedpts)


             else
                 zp2=z(k+ind_distance_fixedpts)
                 quasiset_tt_p2=
     &               quasiset_tt_ll(i,j,k+ind_distance_fixedpts)
                 quasiset_tchi_p2=
     &               quasiset_tchi_ll(i,j,k+ind_distance_fixedpts)
                 quasiset_txi_p2=
     &               quasiset_txi_ll(i,j,k+ind_distance_fixedpts)
                 quasiset_chichi_p2=
     &               quasiset_chichi_ll(i,j,k+ind_distance_fixedpts)
                 quasiset_chixi_p2=
     &               quasiset_chixi_ll(i,j,k+ind_distance_fixedpts)
                 quasiset_xixi_p2=
     &               quasiset_xixi_ll(i,j,k+ind_distance_fixedpts)
                 quasiset_trace_p2=
     &               quasiset_tracell(i,j,k+ind_distance_fixedpts)
                 quasiset_massdensity_p2=
     &               quasiset_massdensityll(i,j,k+ind_distance_fixedpts)
                 quasiset_angmomdensityx_p2=
     &      quasiset_angmomdensityxll(i,j,k+ind_distance_fixedpts)
                 quasiset_angmomdensityy_p2=
     &      quasiset_angmomdensityyll(i,j,k+ind_distance_fixedpts)
                 quasiset_angmomdensityz_p2=
     &      quasiset_angmomdensityzll(i,j,k+ind_distance_fixedpts)



              end if

                 quasiset_tt(lind)=
     &                 firstord_extrap(quasiset_tt_p1
     &                  ,quasiset_tt_p2,zp1,zp2,zex)
                 quasiset_tchi(lind)=
     &                 firstord_extrap(quasiset_tchi_p1
     &                  ,quasiset_tchi_p2,zp1,zp2,zex)
                 quasiset_txi(lind)=
     &                 firstord_extrap(quasiset_txi_p1
     &                  ,quasiset_txi_p2,zp1,zp2,zex)
                 quasiset_chichi(lind)=
     &                 firstord_extrap(quasiset_chichi_p1
     &                  ,quasiset_chichi_p2,zp1,zp2,zex)
                 quasiset_chixi(lind)=
     &                 firstord_extrap(quasiset_chixi_p1
     &                  ,quasiset_chixi_p2,zp1,zp2,zex)
                 quasiset_xixi(lind)=
     &                 firstord_extrap(quasiset_xixi_p1
     &                  ,quasiset_xixi_p2,zp1,zp2,zex)
                 quasiset_trace(lind)=
     &                 firstord_extrap(quasiset_trace_p1
     &                  ,quasiset_trace_p2,zp1,zp2,zex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 firstord_extrap(quasiset_massdensity_p1
     &                  ,quasiset_massdensity_p2,zp1,zp2,zex)
                 quasiset_angmomdensityx(lind)=
     &                 firstord_extrap(quasiset_angmomdensityx_p1
     &                  ,quasiset_angmomdensityx_p2,zp1,zp2,zex)
                 quasiset_angmomdensityy(lind)=
     &                 firstord_extrap(quasiset_angmomdensityy_p1
     &                  ,quasiset_angmomdensityy_p2,zp1,zp2,zex)
                 quasiset_angmomdensityz(lind)=
     &                 firstord_extrap(quasiset_angmomdensityz_p1
     &                  ,quasiset_angmomdensityz_p2,zp1,zp2,zex)
             end if

!              quasiset_massdensity(lind)=sin(PI*chiex)*cos(2*PI*xiex)  !TEST

            end if !closes condition on bdy_extrap_order.eq.1


             if (bdy_extrap_order.eq.2) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then

                  xp2=x(i-ind_distance_fixedpts)
                  xp3=x(i-2*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &             quasiset_tt_ll(i-ind_distance_fixedpts,j,k)
                  quasiset_tt_p3=
     &             quasiset_tt_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_tchi_p2=
     &             quasiset_tchi_ll(i-ind_distance_fixedpts,j,k)
                  quasiset_tchi_p3=
     &             quasiset_tchi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_txi_p2=
     &             quasiset_txi_ll(i-ind_distance_fixedpts,j,k)
                  quasiset_txi_p3=
     &             quasiset_txi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_chichi_p2=
     &             quasiset_chichi_ll(i-ind_distance_fixedpts,j,k)
                  quasiset_chichi_p3=
     &             quasiset_chichi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_chixi_p2=
     &             quasiset_chixi_ll(i-ind_distance_fixedpts,j,k)
                  quasiset_chixi_p3=
     &             quasiset_chixi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_xixi_p2=
     &             quasiset_xixi_ll(i-ind_distance_fixedpts,j,k)
                  quasiset_xixi_p3=
     &             quasiset_xixi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i-ind_distance_fixedpts,j,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i-2*ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i-ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i-ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i-ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i-ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i-2*ind_distance_fixedpts,j,k)


              else
                  xp2=x(i+ind_distance_fixedpts)
                  xp3=x(i+2*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &             quasiset_tt_ll(i+ind_distance_fixedpts,j,k)
                  quasiset_tt_p3=
     &             quasiset_tt_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_tchi_p2=
     &             quasiset_tchi_ll(i+ind_distance_fixedpts,j,k)
                  quasiset_tchi_p3=
     &             quasiset_tchi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_txi_p2=
     &             quasiset_txi_ll(i+ind_distance_fixedpts,j,k)
                  quasiset_txi_p3=
     &             quasiset_txi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_chichi_p2=
     &             quasiset_chichi_ll(i+ind_distance_fixedpts,j,k)
                  quasiset_chichi_p3=
     &             quasiset_chichi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_chixi_p2=
     &             quasiset_chixi_ll(i+ind_distance_fixedpts,j,k)
                  quasiset_chixi_p3=
     &             quasiset_chixi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_xixi_p2=
     &             quasiset_xixi_ll(i+ind_distance_fixedpts,j,k)
                  quasiset_xixi_p3=
     &             quasiset_xixi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i+ind_distance_fixedpts,j,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i+2*ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i+ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i+ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i+ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i+ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i+2*ind_distance_fixedpts,j,k)


               end if

                 quasiset_tt(lind)=
     &                 secondord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_tchi(lind)=
     &                 secondord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_txi(lind)=
     &                 secondord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_chichi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_chixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_xixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_trace(lind)=
     &                 secondord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   xp1,xp2,xp3,xex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 secondord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_angmomdensityx(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_angmomdensityy(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   xp1,xp2,xp3,xex)
                 quasiset_angmomdensityz(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   xp1,xp2,xp3,xex)


             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
              if (yp1.gt.0) then
                  yp2=y(j-ind_distance_fixedpts)
                  yp3=y(j-2*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &             quasiset_tt_ll(i,j-ind_distance_fixedpts,k)
                  quasiset_tt_p3=
     &             quasiset_tt_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_tchi_p2=
     &             quasiset_tchi_ll(i,j-ind_distance_fixedpts,k)
                  quasiset_tchi_p3=
     &             quasiset_tchi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_txi_p2=
     &             quasiset_txi_ll(i,j-ind_distance_fixedpts,k)
                  quasiset_txi_p3=
     &             quasiset_txi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_chichi_p2=
     &             quasiset_chichi_ll(i,j-ind_distance_fixedpts,k)
                  quasiset_chichi_p3=
     &             quasiset_chichi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_chixi_p2=
     &             quasiset_chixi_ll(i,j-ind_distance_fixedpts,k)
                  quasiset_chixi_p3=
     &             quasiset_chixi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_xixi_p2=
     &             quasiset_xixi_ll(i,j-ind_distance_fixedpts,k)
                  quasiset_xixi_p3=
     &             quasiset_xixi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j-ind_distance_fixedpts,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j-2*ind_distance_fixedpts,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j-ind_distance_fixedpts,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i,j-ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i,j-ind_distance_fixedpts,k)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i,j-ind_distance_fixedpts,k)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i,j-2*ind_distance_fixedpts,k)



              else
                  yp2=y(j+ind_distance_fixedpts)
                  yp3=y(j+2*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &             quasiset_tt_ll(i,j+ind_distance_fixedpts,k)
                  quasiset_tt_p3=
     &             quasiset_tt_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_tchi_p2=
     &             quasiset_tchi_ll(i,j+ind_distance_fixedpts,k)
                  quasiset_tchi_p3=
     &             quasiset_tchi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_txi_p2=
     &             quasiset_txi_ll(i,j+ind_distance_fixedpts,k)
                  quasiset_txi_p3=
     &             quasiset_txi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_chichi_p2=
     &             quasiset_chichi_ll(i,j+ind_distance_fixedpts,k)
                  quasiset_chichi_p3=
     &             quasiset_chichi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_chixi_p2=
     &             quasiset_chixi_ll(i,j+ind_distance_fixedpts,k)
                  quasiset_chixi_p3=
     &             quasiset_chixi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_xixi_p2=
     &             quasiset_xixi_ll(i,j+ind_distance_fixedpts,k)
                  quasiset_xixi_p3=
     &             quasiset_xixi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j+ind_distance_fixedpts,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j+2*ind_distance_fixedpts,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j+ind_distance_fixedpts,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i,j+ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i,j+ind_distance_fixedpts,k)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i,j+ind_distance_fixedpts,k)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i,j+2*ind_distance_fixedpts,k)

              end if
                 quasiset_tt(lind)=
     &                 secondord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_tchi(lind)=
     &                 secondord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_txi(lind)=
     &                 secondord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_chichi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_chixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_xixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_trace(lind)=
     &                 secondord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   yp1,yp2,yp3,yex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 secondord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_angmomdensityx(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_angmomdensityy(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   yp1,yp2,yp3,yex)
                 quasiset_angmomdensityz(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   yp1,yp2,yp3,yex)

             else
                 if (zp1.gt.0) then
                  zp2=z(k-ind_distance_fixedpts)
                  zp3=z(k-2*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &             quasiset_tt_ll(i,j,k-ind_distance_fixedpts)
                  quasiset_tt_p3=
     &             quasiset_tt_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_tchi_p2=
     &             quasiset_tchi_ll(i,j,k-ind_distance_fixedpts)
                  quasiset_tchi_p3=
     &             quasiset_tchi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_txi_p2=
     &             quasiset_txi_ll(i,j,k-ind_distance_fixedpts)
                  quasiset_txi_p3=
     &             quasiset_txi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_chichi_p2=
     &             quasiset_chichi_ll(i,j,k-ind_distance_fixedpts)
                  quasiset_chichi_p3=
     &             quasiset_chichi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_chixi_p2=
     &             quasiset_chixi_ll(i,j,k-ind_distance_fixedpts)
                  quasiset_chixi_p3=
     &             quasiset_chixi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_xixi_p2=
     &             quasiset_xixi_ll(i,j,k-ind_distance_fixedpts)
                  quasiset_xixi_p3=
     &             quasiset_xixi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j,k-ind_distance_fixedpts)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j,k-2*ind_distance_fixedpts)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k-ind_distance_fixedpts)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i,j,k-ind_distance_fixedpts)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i,j,k-ind_distance_fixedpts)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i,j,k-ind_distance_fixedpts)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i,j,k-2*ind_distance_fixedpts)

              else
                  zp2=z(k+ind_distance_fixedpts)
                  zp3=z(k+2*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &             quasiset_tt_ll(i,j,k+ind_distance_fixedpts)
                  quasiset_tt_p3=
     &             quasiset_tt_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_tchi_p2=
     &             quasiset_tchi_ll(i,j,k+ind_distance_fixedpts)
                  quasiset_tchi_p3=
     &             quasiset_tchi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_txi_p2=
     &             quasiset_txi_ll(i,j,k+ind_distance_fixedpts)
                  quasiset_txi_p3=
     &             quasiset_txi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_chichi_p2=
     &             quasiset_chichi_ll(i,j,k+ind_distance_fixedpts)
                  quasiset_chichi_p3=
     &             quasiset_chichi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_chixi_p2=
     &             quasiset_chixi_ll(i,j,k+ind_distance_fixedpts)
                  quasiset_chixi_p3=
     &             quasiset_chixi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_xixi_p2=
     &             quasiset_xixi_ll(i,j,k+ind_distance_fixedpts)
                  quasiset_xixi_p3=
     &             quasiset_xixi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j,k+ind_distance_fixedpts)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j,k+2*ind_distance_fixedpts)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k+ind_distance_fixedpts)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i,j,k+ind_distance_fixedpts)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i,j,k+ind_distance_fixedpts)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i,j,k+ind_distance_fixedpts)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i,j,k+2*ind_distance_fixedpts)

               end if
                 quasiset_tt(lind)=
     &                 secondord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_tchi(lind)=
     &                 secondord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_txi(lind)=
     &                 secondord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_chichi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_chixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_xixi(lind)=
     &                 secondord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_trace(lind)=
     &                 secondord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   zp1,zp2,zp3,zex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 secondord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_angmomdensityx(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_angmomdensityy(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   zp1,zp2,zp3,zex)
                 quasiset_angmomdensityz(lind)=
     &                 secondord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   zp1,zp2,zp3,zex)
              end if

!               quasiset_massdensity(lind)=sin(PI*chiex)*cos(2*PI*xiex)  !TEST

            end if !closes condition on bdy_extrap_order.eq.2


             if (bdy_extrap_order.eq.3) then
              if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
                if (xp1.gt.0) then

                  xp2=x(i-1*ind_distance_fixedpts)
                  xp3=x(i-2*ind_distance_fixedpts)
                  xp4=x(i-3*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &                quasiset_tt_ll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_tt_p3=
     &                quasiset_tt_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_tt_p4=
     &                quasiset_tt_ll(i-3*ind_distance_fixedpts,j,k)
                  quasiset_tchi_p2=
     &                quasiset_tchi_ll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_tchi_p3=
     &                quasiset_tchi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_tchi_p4=
     &                quasiset_tchi_ll(i-3*ind_distance_fixedpts,j,k)
                  quasiset_txi_p2=
     &                quasiset_txi_ll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_txi_p3=
     &                quasiset_txi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_txi_p4=
     &                quasiset_txi_ll(i-3*ind_distance_fixedpts,j,k)
                  quasiset_chichi_p2=
     &                quasiset_chichi_ll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_chichi_p3=
     &                quasiset_chichi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_chichi_p4=
     &                quasiset_chichi_ll(i-3*ind_distance_fixedpts,j,k)
                  quasiset_chixi_p2=
     &                quasiset_chixi_ll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_chixi_p3=
     &                quasiset_chixi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_chixi_p4=
     &                quasiset_chixi_ll(i-3*ind_distance_fixedpts,j,k)
                  quasiset_xixi_p2=
     &                quasiset_xixi_ll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_xixi_p3=
     &                quasiset_xixi_ll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_xixi_p4=
     &                quasiset_xixi_ll(i-3*ind_distance_fixedpts,j,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i-1*ind_distance_fixedpts,j,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i-2*ind_distance_fixedpts,j,k)
                  quasiset_trace_p4=
     &             quasiset_tracell(i-3*ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i-3*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p4=
     &       quasiset_angmomdensityxll(i-3*ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityy_p4=
     &       quasiset_angmomdensityyll(i-3*ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i-1*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i-2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityz_p4=
     &       quasiset_angmomdensityzll(i-3*ind_distance_fixedpts,j,k)


              else
                  xp2=x(i+1*ind_distance_fixedpts)
                  xp3=x(i+2*ind_distance_fixedpts)
                  xp4=x(i+3*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &                quasiset_tt_ll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_tt_p3=
     &                quasiset_tt_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_tt_p4=
     &                quasiset_tt_ll(i+3*ind_distance_fixedpts,j,k)
                  quasiset_tchi_p2=
     &                quasiset_tchi_ll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_tchi_p3=
     &                quasiset_tchi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_tchi_p4=
     &                quasiset_tchi_ll(i+3*ind_distance_fixedpts,j,k)
                  quasiset_txi_p2=
     &                quasiset_txi_ll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_txi_p3=
     &                quasiset_txi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_txi_p4=
     &                quasiset_txi_ll(i+3*ind_distance_fixedpts,j,k)
                  quasiset_chichi_p2=
     &                quasiset_chichi_ll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_chichi_p3=
     &                quasiset_chichi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_chichi_p4=
     &                quasiset_chichi_ll(i+3*ind_distance_fixedpts,j,k)
                  quasiset_chixi_p2=
     &                quasiset_chixi_ll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_chixi_p3=
     &                quasiset_chixi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_chixi_p4=
     &                quasiset_chixi_ll(i+3*ind_distance_fixedpts,j,k)
                  quasiset_xixi_p2=
     &                quasiset_xixi_ll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_xixi_p3=
     &                quasiset_xixi_ll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_xixi_p4=
     &                quasiset_xixi_ll(i+3*ind_distance_fixedpts,j,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i+1*ind_distance_fixedpts,j,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i+2*ind_distance_fixedpts,j,k)
                  quasiset_trace_p4=
     &             quasiset_tracell(i+3*ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i+3*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityx_p4=
     &       quasiset_angmomdensityxll(i+3*ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityy_p4=
     &       quasiset_angmomdensityyll(i+3*ind_distance_fixedpts,j,k)
                 quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i+1*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i+2*ind_distance_fixedpts,j,k)
                  quasiset_angmomdensityz_p4=
     &       quasiset_angmomdensityzll(i+3*ind_distance_fixedpts,j,k)


               end if

                 quasiset_tt(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   quasiset_tt_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_tchi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   quasiset_tchi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_txi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   quasiset_txi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_chichi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   quasiset_chichi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_chixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   quasiset_chixi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_xixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   quasiset_xixi_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_trace(lind)=
     &                 thirdord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   quasiset_trace_p4,
     &                   xp1,xp2,xp3,xp4,xex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 thirdord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   quasiset_massdensity_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_angmomdensityx(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   quasiset_angmomdensityx_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_angmomdensityy(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   quasiset_angmomdensityy_p4,
     &                   xp1,xp2,xp3,xp4,xex)
                 quasiset_angmomdensityz(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   quasiset_angmomdensityz_p4,
     &                   xp1,xp2,xp3,xp4,xex)


             else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
                if (yp1.gt.0) then

                  yp2=y(j-1*ind_distance_fixedpts)
                  yp3=y(j-2*ind_distance_fixedpts)
                  yp4=y(j-3*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &                quasiset_tt_ll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_tt_p3=
     &                quasiset_tt_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_tt_p4=
     &                quasiset_tt_ll(i,j-3*ind_distance_fixedpts,k)
                  quasiset_tchi_p2=
     &                quasiset_tchi_ll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_tchi_p3=
     &                quasiset_tchi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_tchi_p4=
     &                quasiset_tchi_ll(i,j-3*ind_distance_fixedpts,k)
                  quasiset_txi_p2=
     &                quasiset_txi_ll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_txi_p3=
     &                quasiset_txi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_txi_p4=
     &                quasiset_txi_ll(i,j-3*ind_distance_fixedpts,k)
                  quasiset_chichi_p2=
     &                quasiset_chichi_ll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_chichi_p3=
     &                quasiset_chichi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_chichi_p4=
     &                quasiset_chichi_ll(i,j-3*ind_distance_fixedpts,k)
                  quasiset_chixi_p2=
     &                quasiset_chixi_ll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_chixi_p3=
     &                quasiset_chixi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_chixi_p4=
     &                quasiset_chixi_ll(i,j-3*ind_distance_fixedpts,k)
                  quasiset_xixi_p2=
     &                quasiset_xixi_ll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_xixi_p3=
     &                quasiset_xixi_ll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_xixi_p4=
     &                quasiset_xixi_ll(i,j-3*ind_distance_fixedpts,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j-1*ind_distance_fixedpts,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j-2*ind_distance_fixedpts,k)
                  quasiset_trace_p4=
     &             quasiset_tracell(i,j-3*ind_distance_fixedpts,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i,j-3*ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p4=
     &       quasiset_angmomdensityxll(i,j-3*ind_distance_fixedpts,k)
                 quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityy_p4=
     &       quasiset_angmomdensityyll(i,j-3*ind_distance_fixedpts,k)
                 quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i,j-1*ind_distance_fixedpts,k)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i,j-2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityz_p4=
     &       quasiset_angmomdensityzll(i,j-3*ind_distance_fixedpts,k)


              else
                  yp2=y(j+1*ind_distance_fixedpts)
                  yp3=y(j+2*ind_distance_fixedpts)
                  yp4=y(j+3*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &                quasiset_tt_ll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_tt_p3=
     &                quasiset_tt_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_tt_p4=
     &                quasiset_tt_ll(i,j+3*ind_distance_fixedpts,k)
                  quasiset_tchi_p2=
     &                quasiset_tchi_ll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_tchi_p3=
     &                quasiset_tchi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_tchi_p4=
     &                quasiset_tchi_ll(i,j+3*ind_distance_fixedpts,k)
                  quasiset_txi_p2=
     &                quasiset_txi_ll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_txi_p3=
     &                quasiset_txi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_txi_p4=
     &                quasiset_txi_ll(i,j+3*ind_distance_fixedpts,k)
                  quasiset_chichi_p2=
     &                quasiset_chichi_ll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_chichi_p3=
     &                quasiset_chichi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_chichi_p4=
     &                quasiset_chichi_ll(i,j+3*ind_distance_fixedpts,k)
                  quasiset_chixi_p2=
     &                quasiset_chixi_ll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_chixi_p3=
     &                quasiset_chixi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_chixi_p4=
     &                quasiset_chixi_ll(i,j+3*ind_distance_fixedpts,k)
                  quasiset_xixi_p2=
     &                quasiset_xixi_ll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_xixi_p3=
     &                quasiset_xixi_ll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_xixi_p4=
     &                quasiset_xixi_ll(i,j+3*ind_distance_fixedpts,k)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j+1*ind_distance_fixedpts,k)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j+2*ind_distance_fixedpts,k)
                  quasiset_trace_p4=
     &             quasiset_tracell(i,j+3*ind_distance_fixedpts,k)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i,j+3*ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityx_p4=
     &       quasiset_angmomdensityxll(i,j+3*ind_distance_fixedpts,k)
                 quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityy_p4=
     &       quasiset_angmomdensityyll(i,j+3*ind_distance_fixedpts,k)
                 quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i,j+1*ind_distance_fixedpts,k)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i,j+2*ind_distance_fixedpts,k)
                  quasiset_angmomdensityz_p4=
     &       quasiset_angmomdensityzll(i,j+3*ind_distance_fixedpts,k)


               end if

                 quasiset_tt(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   quasiset_tt_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_tchi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   quasiset_tchi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_txi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   quasiset_txi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_chichi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   quasiset_chichi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_chixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   quasiset_chixi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_xixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   quasiset_xixi_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_trace(lind)=
     &                 thirdord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   quasiset_trace_p4,
     &                   yp1,yp2,yp3,yp4,yex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 thirdord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   quasiset_massdensity_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_angmomdensityx(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   quasiset_angmomdensityx_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_angmomdensityy(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   quasiset_angmomdensityy_p4,
     &                   yp1,yp2,yp3,yp4,yex)
                 quasiset_angmomdensityz(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityz_p1,
     &                   quasiset_angmomdensityz_p2,
     &                   quasiset_angmomdensityz_p3,
     &                   quasiset_angmomdensityz_p4,
     &                   yp1,yp2,yp3,yp4,yex)
             else
                if (zp1.gt.0) then

                  zp2=z(k-1*ind_distance_fixedpts)
                  zp3=z(k-2*ind_distance_fixedpts)
                  zp4=z(k-3*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &                quasiset_tt_ll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_tt_p3=
     &                quasiset_tt_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_tt_p4=
     &                quasiset_tt_ll(i,j,k-3*ind_distance_fixedpts)
                  quasiset_tchi_p2=
     &                quasiset_tchi_ll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_tchi_p3=
     &                quasiset_tchi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_tchi_p4=
     &                quasiset_tchi_ll(i,j,k-3*ind_distance_fixedpts)
                  quasiset_txi_p2=
     &                quasiset_txi_ll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_txi_p3=
     &                quasiset_txi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_txi_p4=
     &                quasiset_txi_ll(i,j,k-3*ind_distance_fixedpts)
                  quasiset_chichi_p2=
     &                quasiset_chichi_ll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_chichi_p3=
     &                quasiset_chichi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_chichi_p4=
     &                quasiset_chichi_ll(i,j,k-3*ind_distance_fixedpts)
                  quasiset_chixi_p2=
     &                quasiset_chixi_ll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_chixi_p3=
     &                quasiset_chixi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_chixi_p4=
     &                quasiset_chixi_ll(i,j,k-3*ind_distance_fixedpts)
                  quasiset_xixi_p2=
     &                quasiset_xixi_ll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_xixi_p3=
     &                quasiset_xixi_ll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_xixi_p4=
     &                quasiset_xixi_ll(i,j,k-3*ind_distance_fixedpts)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j,k-1*ind_distance_fixedpts)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j,k-2*ind_distance_fixedpts)
                  quasiset_trace_p4=
     &             quasiset_tracell(i,j,k-3*ind_distance_fixedpts)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i,j,k-3*ind_distance_fixedpts)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_angmomdensityx_p4=
     &       quasiset_angmomdensityxll(i,j,k-3*ind_distance_fixedpts)
                 quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_angmomdensityy_p4=
     &       quasiset_angmomdensityyll(i,j,k-3*ind_distance_fixedpts)
                 quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i,j,k-1*ind_distance_fixedpts)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i,j,k-2*ind_distance_fixedpts)
                  quasiset_angmomdensityz_p4=
     &       quasiset_angmomdensityzll(i,j,k-3*ind_distance_fixedpts)


              else
                  zp2=z(k+1*ind_distance_fixedpts)
                  zp3=z(k+2*ind_distance_fixedpts)
                  zp4=z(k+3*ind_distance_fixedpts)
                  quasiset_tt_p2=
     &                quasiset_tt_ll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_tt_p3=
     &                quasiset_tt_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_tt_p4=
     &                quasiset_tt_ll(i,j,k+3*ind_distance_fixedpts)
                  quasiset_tchi_p2=
     &                quasiset_tchi_ll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_tchi_p3=
     &                quasiset_tchi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_tchi_p4=
     &                quasiset_tchi_ll(i,j,k+3*ind_distance_fixedpts)
                  quasiset_txi_p2=
     &                quasiset_txi_ll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_txi_p3=
     &                quasiset_txi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_txi_p4=
     &                quasiset_txi_ll(i,j,k+3*ind_distance_fixedpts)
                  quasiset_chichi_p2=
     &                quasiset_chichi_ll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_chichi_p3=
     &                quasiset_chichi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_chichi_p4=
     &                quasiset_chichi_ll(i,j,k+3*ind_distance_fixedpts)
                  quasiset_chixi_p2=
     &                quasiset_chixi_ll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_chixi_p3=
     &                quasiset_chixi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_chixi_p4=
     &                quasiset_chixi_ll(i,j,k+3*ind_distance_fixedpts)
                  quasiset_xixi_p2=
     &                quasiset_xixi_ll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_xixi_p3=
     &                quasiset_xixi_ll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_xixi_p4=
     &                quasiset_xixi_ll(i,j,k+3*ind_distance_fixedpts)
                  quasiset_trace_p2=
     &             quasiset_tracell(i,j,k+1*ind_distance_fixedpts)
                  quasiset_trace_p3=
     &             quasiset_tracell(i,j,k+2*ind_distance_fixedpts)
                  quasiset_trace_p4=
     &             quasiset_tracell(i,j,k+3*ind_distance_fixedpts)
                  quasiset_massdensity_p2=
     &             quasiset_massdensityll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_massdensity_p3=
     &             quasiset_massdensityll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_massdensity_p4=
     &             quasiset_massdensityll(i,j,k+3*ind_distance_fixedpts)
                  quasiset_angmomdensityx_p2=
     &       quasiset_angmomdensityxll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_angmomdensityx_p3=
     &       quasiset_angmomdensityxll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_angmomdensityx_p4=
     &       quasiset_angmomdensityxll(i,j,k+3*ind_distance_fixedpts)
                 quasiset_angmomdensityy_p2=
     &       quasiset_angmomdensityyll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_angmomdensityy_p3=
     &       quasiset_angmomdensityyll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_angmomdensityy_p4=
     &       quasiset_angmomdensityyll(i,j,k+3*ind_distance_fixedpts)
                 quasiset_angmomdensityz_p2=
     &       quasiset_angmomdensityzll(i,j,k+1*ind_distance_fixedpts)
                  quasiset_angmomdensityz_p3=
     &       quasiset_angmomdensityzll(i,j,k+2*ind_distance_fixedpts)
                  quasiset_angmomdensityz_p4=
     &       quasiset_angmomdensityzll(i,j,k+3*ind_distance_fixedpts)


               end if

                 quasiset_tt(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tt_p1,
     &                   quasiset_tt_p2,
     &                   quasiset_tt_p3,
     &                   quasiset_tt_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_tchi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_tchi_p1,
     &                   quasiset_tchi_p2,
     &                   quasiset_tchi_p3,
     &                   quasiset_tchi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_txi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_txi_p1,
     &                   quasiset_txi_p2,
     &                   quasiset_txi_p3,
     &                   quasiset_txi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_chichi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chichi_p1,
     &                   quasiset_chichi_p2,
     &                   quasiset_chichi_p3,
     &                   quasiset_chichi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_chixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_chixi_p1,
     &                   quasiset_chixi_p2,
     &                   quasiset_chixi_p3,
     &                   quasiset_chixi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_xixi(lind)=
     &                 thirdord_extrap(
     &                   quasiset_xixi_p1,
     &                   quasiset_xixi_p2,
     &                   quasiset_xixi_p3,
     &                   quasiset_xixi_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_trace(lind)=
     &                 thirdord_extrap(
     &                   quasiset_trace_p1,
     &                   quasiset_trace_p2,
     &                   quasiset_trace_p3,
     &                   quasiset_trace_p4,
     &                   zp1,zp2,zp3,zp4,zex)
!                 quasiset_trace(lind)=
!     &           (
!     &           gamma0qssphbdy_uu_tt*quasiset_tt(lind)
!     &        +2*gamma0qssphbdy_uu_tchi*quasiset_tchi(lind)
!     &        +2*gamma0qssphbdy_uu_txi*quasiset_txi(lind)
!     &          +gamma0qssphbdy_uu_chichi*quasiset_chichi(lind)
!     &        +2*gamma0qssphbdy_uu_chixi*quasiset_chixi(lind)
!     &          +gamma0qssphbdy_uu_xixi*quasiset_xixi(lind)
!     &                   )
                 quasiset_massdensity(lind)=
     &                 thirdord_extrap(
     &                   quasiset_massdensity_p1,
     &                   quasiset_massdensity_p2,
     &                   quasiset_massdensity_p3,
     &                   quasiset_massdensity_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_angmomdensityx(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityx_p1,
     &                   quasiset_angmomdensityx_p2,
     &                   quasiset_angmomdensityx_p3,
     &                   quasiset_angmomdensityx_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_angmomdensityy(lind)=
     &                 thirdord_extrap(
     &                   quasiset_angmomdensityy_p1,
     &                   quasiset_angmomdensityy_p2,
     &                   quasiset_angmomdensityy_p3,
     &                   quasiset_angmomdensityy_p4,
     &                   zp1,zp2,zp3,zp4,zex)
                 quasiset_angmomdensityz(lind)=
     &                 thirdord_extrap(
     &           quasiset_angmomdensityz_p1,
     &           quasiset_angmomdensityz_p2,
     &           quasiset_angmomdensityz_p3,
     &           quasiset_angmomdensityz_p4,
     &                   zp1,zp2,zp3,zp4,zex)
              end if

!               quasiset_massdensity(lind)=sin(PI*chiex)*cos(2*PI*xiex)  !TEST














!            if ((abs(maxxyzp1-abs(xp1)).lt.10.0d0**(-10))) then
!             if (xp1.gt.0) then
!               write(*,*) "EXTRAP_ORDER 3
!     &       i,
!     &       x(i),
!     &       quasiset_tt_ll(i,j,k),
!     &       quasiset_tchi_ll(i,j,k),
!     &       quasiset_txi_ll(i,j,k),
!     &       quasiset_chichi_ll(i,j,k),
!     &       quasiset_chixi_ll(i,j,k),
!     &       quasiset_xixi_ll(i,j,k),
!     &       i-1*ind_distance_fixedpts,
!     &       x(i-1*ind_distance_fixedpts),
!     &       quasiset_tt_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_tchi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_txi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_chichi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_chixi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_xixi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       i-2*ind_distance_fixedpts,
!     &       x(i-2*ind_distance_fixedpts),
!     &       quasiset_tt_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_tchi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_txi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_chichi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_chixi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_xixi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       i-3*ind_distance_fixedpts,
!     &       x(i-3*ind_distance_fixedpts),
!     &       quasiset_tt_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_tchi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_txi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_chichi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_chixi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_xixi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       lind,
!     &       xex,
!     &       yex,
!     &       zex,
!     &       quasiset_tt(lind),
!     &       quasiset_tchi(lind),
!     &       quasiset_txi(lind),
!     &       quasiset_chichi(lind),
!     &       quasiset_chixi(lind),
!     &       quasiset_xixi(lind),
!     &       quasiset_trace(lind)=",
!     &       i,
!     &       x(i),
!     &       quasiset_tt_ll(i,j,k),
!     &       quasiset_tchi_ll(i,j,k),
!     &       quasiset_txi_ll(i,j,k),
!     &       quasiset_chichi_ll(i,j,k),
!     &       quasiset_chixi_ll(i,j,k),
!     &       quasiset_xixi_ll(i,j,k),
!     &       i-1*ind_distance_fixedpts,
!     &       x(i-1*ind_distance_fixedpts),
!     &       quasiset_tt_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_tchi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_txi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_chichi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_chixi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       quasiset_xixi_ll(i-1*ind_distance_fixedpts,j,k),
!     &       i-2*ind_distance_fixedpts,
!     &       x(i-2*ind_distance_fixedpts),
!     &       quasiset_tt_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_tchi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_txi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_chichi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_chixi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       quasiset_xixi_ll(i-2*ind_distance_fixedpts,j,k),
!     &       i-3*ind_distance_fixedpts,
!     &       x(i-3*ind_distance_fixedpts),
!     &       quasiset_tt_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_tchi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_txi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_chichi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_chixi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       quasiset_xixi_ll(i-3*ind_distance_fixedpts,j,k),
!     &       lind,
!     &       xex,
!     &       yex,
!     &       zex,
!     &       quasiset_tt(lind),
!     &       quasiset_tchi(lind),
!     &       quasiset_txi(lind),
!     &       quasiset_chichi(lind),
!     &       quasiset_chixi(lind),
!     &       quasiset_xixi(lind),
!     &       quasiset_trace(lind)
!
!             else
!               write(*,*) "EXTRAP_ORDER 3
!     &         i,
!     &         x(i),
!     &         quasiset_tt_ll(i,j,k),
!     &         quasiset_tchi_ll(i,j,k),
!     &         quasiset_txi_ll(i,j,k),
!     &         quasiset_chichi_ll(i,j,k),
!     &         quasiset_chixi_ll(i,j,k),
!     &         quasiset_xixi_ll(i,j,k),
!     &         i+1*ind_distance_fixedpts,
!     &         x(i+1*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_tchi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_txi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_chichi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_chixi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_xixi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         i+2*ind_distance_fixedpts,
!     &         x(i+2*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_tchi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_txi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_chichi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_chixi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_xixi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         i+3*ind_distance_fixedpts,
!     &         x(i+3*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_tchi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_txi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_chichi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_chixi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_xixi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         lind,
!     &         xex,
!     &         yex,
!     &         zex,
!     &         quasiset_tt(lind),
!     &         quasiset_tchi(lind),
!     &         quasiset_txi(lind),
!     &         quasiset_chichi(lind),
!     &         quasiset_chixi(lind),
!     &         quasiset_xixi(lind),
!     &         quasiset_trace(lind)=",
!     &         i,
!     &         x(i),
!     &         quasiset_tt_ll(i,j,k),
!     &         quasiset_tchi_ll(i,j,k),
!     &         quasiset_txi_ll(i,j,k),
!     &         quasiset_chichi_ll(i,j,k),
!     &         quasiset_chixi_ll(i,j,k),
!     &         quasiset_xixi_ll(i,j,k),
!     &         i+1*ind_distance_fixedpts,
!     &         x(i+1*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_tchi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_txi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_chichi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_chixi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         quasiset_xixi_ll(i+1*ind_distance_fixedpts,j,k),
!     &         i+2*ind_distance_fixedpts,
!     &         x(i+2*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_tchi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_txi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_chichi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_chixi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         quasiset_xixi_ll(i+2*ind_distance_fixedpts,j,k),
!     &         i+3*ind_distance_fixedpts,
!     &         x(i+3*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_tchi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_txi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_chichi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_chixi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         quasiset_xixi_ll(i+3*ind_distance_fixedpts,j,k),
!     &         lind,
!     &         xex,
!     &         yex,
!     &         zex,
!     &         quasiset_tt(lind),
!     &         quasiset_tchi(lind),
!     &         quasiset_txi(lind),
!     &         quasiset_chichi(lind),
!     &         quasiset_chixi(lind),
!     &         quasiset_xixi(lind),
!     &         quasiset_trace(lind)
!             end if
!           else if ((abs(maxxyzp1-abs(yp1)).lt.10.0d0**(-10))) then
!            if (yp1.gt.0) then
!               write(*,*) "EXTRAP_ORDER 3
!     &         j,
!     &         y(j),
!     &         quasiset_tt_ll(i,j,k),
!     &         quasiset_tchi_ll(i,j,k),
!     &         quasiset_txi_ll(i,j,k),
!     &         quasiset_chichi_ll(i,j,k),
!     &         quasiset_chixi_ll(i,j,k),
!     &         quasiset_xixi_ll(i,j,k),
!     &         j-1*ind_distance_fixedpts,
!     &         y(j-1*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_tchi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_txi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_chichi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_chixi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_xixi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         j-2*ind_distance_fixedpts,
!     &         y(j-2*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_tchi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_txi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_chichi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_chixi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_xixi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         j-3*ind_distance_fixedpts,
!     &         y(j-3*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_tchi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_txi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_chichi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_chixi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_xixi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         lind,
!     &         xex,
!     &         yex,
!     &         zex,
!     &         quasiset_tt(lind),
!     &         quasiset_tchi(lind),
!     &         quasiset_txi(lind),
!     &         quasiset_chichi(lind),
!     &         quasiset_chixi(lind),
!     &         quasiset_xixi(lind),
!     &         quasiset_trace(lind)=",
!     &         j,
!     &         y(j),
!     &         quasiset_tt_ll(i,j,k),
!     &         quasiset_tchi_ll(i,j,k),
!     &         quasiset_txi_ll(i,j,k),
!     &         quasiset_chichi_ll(i,j,k),
!     &         quasiset_chixi_ll(i,j,k),
!     &         quasiset_xixi_ll(i,j,k),
!     &         j-1*ind_distance_fixedpts,
!     &         y(j-1*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_tchi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_txi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_chichi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_chixi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         quasiset_xixi_ll(i,j-1*ind_distance_fixedpts,k),
!     &         j-2*ind_distance_fixedpts,
!     &         y(j-2*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_tchi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_txi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_chichi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_chixi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         quasiset_xixi_ll(i,j-2*ind_distance_fixedpts,k),
!     &         j-3*ind_distance_fixedpts,
!     &         y(j-3*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_tchi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_txi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_chichi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_chixi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         quasiset_xixi_ll(i,j-3*ind_distance_fixedpts,k),
!     &         lind,
!     &         xex,
!     &         yex,
!     &         zex,
!     &         quasiset_tt(lind),
!     &         quasiset_tchi(lind),
!     &         quasiset_txi(lind),
!     &         quasiset_chichi(lind),
!     &         quasiset_chixi(lind),
!     &         quasiset_xixi(lind),
!     &         quasiset_trace(lind)
!             else
!               write(*,*) "EXTRAP_ORDER 3
!     &        j,
!     &        y(j),
!     &        quasiset_tt_ll(i,j,k),
!     &        quasiset_tchi_ll(i,j,k),
!     &        quasiset_txi_ll(i,j,k),
!     &        quasiset_chichi_ll(i,j,k),
!     &        quasiset_chixi_ll(i,j,k),
!     &        quasiset_xixi_ll(i,j,k),
!     &        j+1*ind_distance_fixedpts,
!     &        y(j+1*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_tchi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_txi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_chichi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_chixi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_xixi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        j+2*ind_distance_fixedpts,
!     &        y(j+2*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_tchi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_txi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_chichi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_chixi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_xixi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        j+3*ind_distance_fixedpts,
!     &        y(j+3*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_tchi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_txi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_chichi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_chixi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_xixi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        lind,
!     &        xex,
!     &        yex,
!     &        zex,
!     &        quasiset_tt(lind),
!     &        quasiset_tchi(lind),
!     &        quasiset_txi(lind),
!     &        quasiset_chichi(lind),
!     &        quasiset_chixi(lind),
!     &        quasiset_xixi(lind),
!     &        quasiset_trace(lind)=",
!     &        j,
!     &        y(j),
!     &        quasiset_tt_ll(i,j,k),
!     &        quasiset_tchi_ll(i,j,k),
!     &        quasiset_txi_ll(i,j,k),
!     &        quasiset_chichi_ll(i,j,k),
!     &        quasiset_chixi_ll(i,j,k),
!     &        quasiset_xixi_ll(i,j,k),
!     &        j+1*ind_distance_fixedpts,
!     &        y(j+1*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_tchi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_txi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_chichi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_chixi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        quasiset_xixi_ll(i,j+1*ind_distance_fixedpts,k),
!     &        j+2*ind_distance_fixedpts,
!     &        y(j+2*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_tchi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_txi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_chichi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_chixi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        quasiset_xixi_ll(i,j+2*ind_distance_fixedpts,k),
!     &        j+3*ind_distance_fixedpts,
!     &        y(j+3*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_tchi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_txi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_chichi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_chixi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        quasiset_xixi_ll(i,j+3*ind_distance_fixedpts,k),
!     &        lind,
!     &        xex,
!     &        yex,
!     &        zex,
!     &        quasiset_tt(lind),
!     &        quasiset_tchi(lind),
!     &        quasiset_txi(lind),
!     &        quasiset_chichi(lind),
!     &        quasiset_chixi(lind),
!     &        quasiset_xixi(lind),
!     &        quasiset_trace(lind)
!             end if
!           else
!            if (zp1.gt.0) then
!               write(*,*) "EXTRAP_ORDER 3
!     &        k,
!     &        z(k),
!     &        quasiset_tt_ll(i,j,k),
!     &        quasiset_tchi_ll(i,j,k),
!     &        quasiset_txi_ll(i,j,k),
!     &        quasiset_chichi_ll(i,j,k),
!     &        quasiset_chixi_ll(i,j,k),
!     &        quasiset_xixi_ll(i,j,k),
!     &        k-1*ind_distance_fixedpts,
!     &        z(k-1*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_tchi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_txi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_chichi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_chixi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_xixi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        k-2*ind_distance_fixedpts,
!     &        z(k-2*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_tchi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_txi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_chichi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_chixi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_xixi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        k-3*ind_distance_fixedpts,
!     &        z(k-3*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_tchi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_txi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_chichi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_chixi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_xixi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        lind,
!     &        xex,
!     &        yex,
!     &        zex,
!     &        quasiset_tt(lind),
!     &        quasiset_tchi(lind),
!     &        quasiset_txi(lind),
!     &        quasiset_chichi(lind),
!     &        quasiset_chixi(lind),
!     &        quasiset_xixi(lind),
!     &        quasiset_trace(lind)=",
!     &        k,
!     &        z(k),
!     &        quasiset_tt_ll(i,j,k),
!     &        quasiset_tchi_ll(i,j,k),
!     &        quasiset_txi_ll(i,j,k),
!     &        quasiset_chichi_ll(i,j,k),
!     &        quasiset_chixi_ll(i,j,k),
!     &        quasiset_xixi_ll(i,j,k),
!     &        k-1*ind_distance_fixedpts,
!     &        z(k-1*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_tchi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_txi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_chichi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_chixi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        quasiset_xixi_ll(i,j,k-1*ind_distance_fixedpts),
!     &        k-2*ind_distance_fixedpts,
!     &        z(k-2*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_tchi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_txi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_chichi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_chixi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        quasiset_xixi_ll(i,j,k-2*ind_distance_fixedpts),
!     &        k-3*ind_distance_fixedpts,
!     &        z(k-3*ind_distance_fixedpts),
!     &        quasiset_tt_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_tchi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_txi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_chichi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_chixi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        quasiset_xixi_ll(i,j,k-3*ind_distance_fixedpts),
!     &        lind,
!     &        xex,
!     &        yex,
!     &        zex,
!     &        quasiset_tt(lind),
!     &        quasiset_tchi(lind),
!     &        quasiset_txi(lind),
!     &        quasiset_chichi(lind),
!     &        quasiset_chixi(lind),
!     &        quasiset_xixi(lind),
!     &        quasiset_trace(lind)
!             else
!               write(*,*) "EXTRAP_ORDER 3
!     &         k,
!     &         z(k),
!     &         quasiset_tt_ll(i,j,k),
!     &         quasiset_tchi_ll(i,j,k),
!     &         quasiset_txi_ll(i,j,k),
!     &         quasiset_chichi_ll(i,j,k),
!     &         quasiset_chixi_ll(i,j,k),
!     &         quasiset_xixi_ll(i,j,k),
!     &         k+1*ind_distance_fixedpts,
!     &         z(k+1*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_tchi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_txi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_chichi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_chixi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_xixi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         k+2*ind_distance_fixedpts,
!     &         z(k+2*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_tchi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_txi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_chichi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_chixi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_xixi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         k+3*ind_distance_fixedpts,
!     &         z(k+3*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_tchi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_txi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_chichi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_chixi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_xixi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         lind,
!     &         xex,
!     &         yex,
!     &         zex,
!     &         quasiset_tt(lind),
!     &         quasiset_tchi(lind),
!     &         quasiset_txi(lind),
!     &         quasiset_chichi(lind),
!     &         quasiset_chixi(lind),
!     &         quasiset_xixi(lind),
!     &         quasiset_trace(lind)=",
!     &         k,
!     &         z(k),
!     &         quasiset_tt_ll(i,j,k),
!     &         quasiset_tchi_ll(i,j,k),
!     &         quasiset_txi_ll(i,j,k),
!     &         quasiset_chichi_ll(i,j,k),
!     &         quasiset_chixi_ll(i,j,k),
!     &         quasiset_xixi_ll(i,j,k),
!     &         k+1*ind_distance_fixedpts,
!     &         z(k+1*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_tchi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_txi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_chichi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_chixi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         quasiset_xixi_ll(i,j,k+1*ind_distance_fixedpts),
!     &         k+2*ind_distance_fixedpts,
!     &         z(k+2*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_tchi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_txi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_chichi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_chixi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         quasiset_xixi_ll(i,j,k+2*ind_distance_fixedpts),
!     &         k+3*ind_distance_fixedpts,
!     &         z(k+3*ind_distance_fixedpts),
!     &         quasiset_tt_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_tchi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_txi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_chichi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_chixi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         quasiset_xixi_ll(i,j,k+3*ind_distance_fixedpts),
!     &         lind,
!     &         xex,
!     &         yex,
!     &         zex,
!     &         quasiset_tt(lind),
!     &         quasiset_tchi(lind),
!     &         quasiset_txi(lind),
!     &         quasiset_chichi(lind),
!     &         quasiset_chixi(lind),
!     &         quasiset_xixi(lind),
!     &         quasiset_trace(lind)
!
!             end if
!            end if




















            end if !closes condition on bdy_extrap_order.eq.3




           end if !closes condition on chrbdy(i,j,k).ne.ex
         end do
        end do
       end do


!       do lind=1,numbdypoints
!        write(*,*) "lind,xpbdy(lind),ypbdy(lind),zpbdy(lind)="
!     &             ,lind,xpbdy(lind),ypbdy(lind),zpbdy(lind)
!        write(*,*) "quasiset_trace(lind)=",quasiset_trace(lind)
!       end do

       return
       end
c--------------------------------------------------------------------------------------


c----------------------------------------------------------------------
c 2-point (i.e. first order) extrapolation stencil using values Tp1,Tp2 at xp1,xp2 (p1 is the closest point to the one we want to extrapolate, p2 is the furthest one)
c----------------------------------------------------------------------
        real*8 function firstord_extrap(Tp1,Tp2,xp1,xp2,xex)
        implicit none
        real*8 Tp1,Tp2,xp1,xp2,xex

        !--------------------------------------------------------------

        firstord_extrap=Tp1*(xex-xp2)/(xp1-xp2)+Tp2*(xex-xp1)/(xp2-xp1)

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c 3-point (i.e. second order) extrapolation stencil using values Tp1,Tp2,Tp3 at xp1,xp2,xp3 (p1 is the closest point to the one we want to extrapolate, p2 is the middle one, p3 is the furthest one)
c----------------------------------------------------------------------
        real*8 function secondord_extrap(Tp1,Tp2,Tp3,xp1,xp2,xp3,xex)
        implicit none
        real*8 Tp1,Tp2,Tp3,xp1,xp2,xp3,xex

        !--------------------------------------------------------------

        secondord_extrap=Tp1*(xex-xp2)*(xex-xp3)/((xp1-xp2)*(xp1-xp3))
     &                  +Tp2*(xex-xp1)*(xex-xp3)/((xp2-xp1)*(xp2-xp3))
     &                  +Tp3*(xex-xp1)*(xex-xp2)/((xp3-xp1)*(xp3-xp2))

        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c 4-point (i.e. third order) extrapolation stencil using values Tp1,Tp2,Tp3,Tp4 at xp1,xp2,xp3,xp4 (if xp1>0, the order is xp4<xp3<xp2<xp1<xex; if xp1<0, the order is xex<xp1<xp2<xp3<xp4)
c----------------------------------------------------------------------
        real*8 function thirdord_extrap(Tp1,Tp2,Tp3,Tp4
     &                                 ,xp1,xp2,xp3,xp4,xex)
        implicit none
        real*8 Tp1,Tp2,Tp3,Tp4,xp1,xp2,xp3,xp4,xex

        !--------------------------------------------------------------

        thirdord_extrap= Tp1*(xex-xp2)*(xex-xp3)*(xex-xp4)
     &                      /((xp1-xp2)*(xp1-xp3)*(xp1-xp4))
     &                  +Tp2*(xex-xp1)*(xex-xp3)*(xex-xp4)
     &                      /((xp2-xp1)*(xp2-xp3)*(xp2-xp4))
     &                  +Tp3*(xex-xp1)*(xex-xp2)*(xex-xp4)
     &                      /((xp3-xp1)*(xp3-xp2)*(xp3-xp4))
     &                  +Tp4*(xex-xp1)*(xex-xp2)*(xex-xp3)
     &                      /((xp4-xp1)*(xp4-xp2)*(xp4-xp3))
        
        return
        end
c--------------------------------------------------------------------------------------

c----------------------------------------------------------------------
c 5-point (i.e. fourth order) extrapolation stencil using values Tp1,Tp2,Tp3,Tp4,Tp5 at xp1,xp2,xp3,xp4,xp5 (if xp1>0, the order is xp5<xp4<xp3<xp2<xp1<xex; if xp1<0, the order is xex<xp1<xp2<xp3<xp4<xp5)
c----------------------------------------------------------------------
        real*8 function fourthord_extrap(Tp1,Tp2,Tp3,Tp4,Tp5
     &                                 ,xp1,xp2,xp3,xp4,xp5,xex)
        implicit none
        real*8 Tp1,Tp2,Tp3,Tp4,Tp5,xp1,xp2,xp3,xp4,xp5,xex

        !--------------------------------------------------------------

        fourthord_extrap= Tp1*(xex-xp2)*(xex-xp3)*(xex-xp4)*(xex-xp5)
     &                      /((xp1-xp2)*(xp1-xp3)*(xp1-xp4)*(xp1-xp5))
     &                   +Tp2*(xex-xp1)*(xex-xp3)*(xex-xp4)*(xex-xp5)
     &                      /((xp2-xp1)*(xp2-xp3)*(xp2-xp4)*(xp2-xp5))
     &                   +Tp3*(xex-xp1)*(xex-xp2)*(xex-xp4)*(xex-xp5)
     &                      /((xp3-xp1)*(xp3-xp2)*(xp3-xp4)*(xp3-xp5))
     &                   +Tp4*(xex-xp1)*(xex-xp2)*(xex-xp3)*(xex-xp5)
     &                      /((xp4-xp1)*(xp4-xp2)*(xp4-xp3)*(xp4-xp5))
     &                   +Tp5*(xex-xp1)*(xex-xp2)*(xex-xp3)*(xex-xp4)
     &                      /((xp5-xp1)*(xp5-xp2)*(xp5-xp3)*(xp5-xp4))

        return
        end
c--------------------------------------------------------------------------------------

c------------------------------------------------------------------------------------------------------------------------------
c The following routine performs the double integral on the S^2 boundary. We approximate the integrals via the trapezoidal rule
c------------------------------------------------------------------------------------------------------------------------------

        subroutine doubleintegralonsphere(integral,density,
     &                  xpbdy,ypbdy,zpbdy,numbdypoints,
     &                  rhobdy,chibdy,xibdy,
     &                  bdy_Nchi,bdy_Nxi)

!----------------------------------------------------------------------

        implicit none
        integer numbdypoints
        real*8 density(numbdypoints)
        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 integral

        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        integer bdy_Nchi,bdy_Nxi

        real*8 rhobdy
        real*8 chibdy(bdy_Nchi),xibdy(bdy_Nxi)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer e,i,j,lind
        integer lind_chipxip,lind_chipxipp1
        integer lind_chipp1xip,lind_chipp1xipp1
        integer additions

        real*8 rhoextrap,rho2

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        real*8 x_chipxip,y_chipxip,z_chipxip
        real*8 x_chipxipp1,y_chipxipp1,z_chipxipp1
        real*8 x_chipp1xip,y_chipp1xip,z_chipp1xip
        real*8 x_chipp1xipp1,y_chipp1xipp1,z_chipp1xipp1

        integer trig_chipxip,trig_chipxipp1
        integer trig_chipp1xip,trig_chipp1xipp1

        real*8 dist_extr_bdypp_min,dist_extr_bdyppp1_min
        real*8 dist_extr_bdypp1p_min,dist_extr_bdypp1pp1_min

        real*8 dist_extr_bdypp_a,dist_extr_bdyppp1_a
        real*8 dist_extr_bdypp1p_a,dist_extr_bdypp1pp1_a


        !-----------------------------------------------------------------------

        integral=0.0d0
        do i=1,bdy_Nchi-1
         do j=1,bdy_Nxi-1
           x_chipxip=rhobdy*cos(PI*chibdy(i))
           y_chipxip=rhobdy*sin(PI*chibdy(i))*cos(2*PI*xibdy(j))
           z_chipxip=rhobdy*sin(PI*chibdy(i))*sin(2*PI*xibdy(j))

           x_chipxipp1=rhobdy*cos(PI*chibdy(i))
           y_chipxipp1=rhobdy*sin(PI*chibdy(i))*cos(2*PI*xibdy(j+1))
           z_chipxipp1=rhobdy*sin(PI*chibdy(i))*sin(2*PI*xibdy(j+1))

           x_chipp1xip=rhobdy*cos(PI*chibdy(i+1))
           y_chipp1xip=rhobdy*sin(PI*chibdy(i+1))*cos(2*PI*xibdy(j))
           z_chipp1xip=rhobdy*sin(PI*chibdy(i+1))*sin(2*PI*xibdy(j))

           x_chipp1xipp1=rhobdy*cos(PI*chibdy(i+1))
           y_chipp1xipp1=rhobdy*sin(PI*chibdy(i+1))*cos(2*PI*xibdy(j+1))
           z_chipp1xipp1=rhobdy*sin(PI*chibdy(i+1))*sin(2*PI*xibdy(j+1))

             dist_extr_bdypp_min=sqrt((xpbdy(1)-x_chipxip)**2
     &                              +(ypbdy(1)-y_chipxip)**2
     &                              +(zpbdy(1)-z_chipxip)**2)
             lind_chipxip=1

             dist_extr_bdyppp1_min=sqrt((xpbdy(1)-x_chipxipp1)**2
     &                              +(ypbdy(1)-y_chipxipp1)**2
     &                              +(zpbdy(1)-z_chipxipp1)**2)
             lind_chipxipp1=1

             dist_extr_bdypp1p_min=sqrt((xpbdy(1)-x_chipp1xip)**2
     &                              +(ypbdy(1)-y_chipp1xip)**2
     &                              +(zpbdy(1)-z_chipp1xip)**2)
             lind_chipp1xip=1

             dist_extr_bdypp1pp1_min=
     &                             sqrt((xpbdy(1)-x_chipp1xipp1)**2
     &                              +(ypbdy(1)-y_chipp1xipp1)**2
     &                              +(zpbdy(1)-z_chipp1xipp1)**2)
             lind_chipp1xipp1=1

           do lind=2,numbdypoints

              dist_extr_bdypp_a=sqrt((xpbdy(lind)-x_chipxip)**2
     &                              +(ypbdy(lind)-y_chipxip)**2
     &                              +(zpbdy(lind)-z_chipxip)**2)
              if (dist_extr_bdypp_a.lt.dist_extr_bdypp_min) then
                  dist_extr_bdypp_min=dist_extr_bdypp_a
                  lind_chipxip=lind
              end if

              dist_extr_bdyppp1_a=sqrt((xpbdy(lind)-x_chipxipp1)**2
     &                              +(ypbdy(lind)-y_chipxipp1)**2
     &                              +(zpbdy(lind)-z_chipxipp1)**2)
              if (dist_extr_bdyppp1_a.lt.dist_extr_bdyppp1_min) then
                  dist_extr_bdyppp1_min=dist_extr_bdyppp1_a
                  lind_chipxipp1=lind
              end if

              dist_extr_bdypp1p_a=sqrt((xpbdy(lind)-x_chipp1xip)**2
     &                              +(ypbdy(lind)-y_chipp1xip)**2
     &                              +(zpbdy(lind)-z_chipp1xip)**2)
              if (dist_extr_bdypp1p_a.lt.dist_extr_bdypp1p_min) then
                  dist_extr_bdypp1p_min=dist_extr_bdypp1p_a
                   lind_chipp1xip=lind
              end if

             dist_extr_bdypp1pp1_a=
     &        sqrt((xpbdy(lind)-x_chipp1xipp1)**2
     &         +(ypbdy(lind)-y_chipp1xipp1)**2
     &         +(zpbdy(lind)-z_chipp1xipp1)**2)
              if (dist_extr_bdypp1pp1_a.lt.dist_extr_bdypp1pp1_min) then
                  dist_extr_bdypp1pp1_min=dist_extr_bdypp1pp1_a
                   lind_chipp1xipp1=lind
              end if
            end do

              integral=integral+
     &              (chibdy(i+1)-chibdy(i))/2 * (xibdy(j+1)-xibdy(j))/2
     &     *(2*PI**2*sin(PI*chibdy(i))*density(lind_chipxip)
     &      +2*PI**2*sin(PI*chibdy(i))*density(lind_chipxipp1)
     &      +2*PI**2*sin(PI*chibdy(i+1))*density(lind_chipp1xip)
     &      +2*PI**2*sin(PI*chibdy(i+1))*density(lind_chipp1xipp1))

         end do
        end do
 
        return
        end
!----------------------------------------------------------------------

c--------------------------------------------------------------------------------------
c Routine to calculate the value of angular coordinates chi and xi at the boundary points where we extrapolate the quasi-local stress energy tensor
c--------------------------------------------------------------------------------------

        subroutine chixiextrap(
     &                         rhoextrap,chiextrap,xiextrap,
     &                         xpbdy,ypbdy,zpbdy,
     &                         numbdypoints)

!----------------------------------------------------------------------

        integer numbdypoints

        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 rhoextrap
        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer e

                !--------------------------------------------------------------

        rhoextrap=1.0d0
        do e=1,numbdypoints
         chiextrap(e)=(1/PI)*acos(xpbdy(e)/rhoextrap)
         if (zpbdy(e).lt.0) then
             xiextrap(e)=(1/(2*PI))
     &        *(atan2(zpbdy(e),ypbdy(e))+2*PI)
         else
             xiextrap(e)=(1/(2*PI))*atan2(zpbdy(e),ypbdy(e))
         end if
        end do

        return
        end
!----------------------------------------------------------------------


c---------------------------------------------------------------------------------------------------------------------------------------------------------
c The following routine calculate the number of different values taken by chi and xi at the boundary points where the stress-energy tensor is extrapolated
c---------------------------------------------------------------------------------------------------------------------------------------------------------

        subroutine bdyn(
     &                  bdy_Nchi,bdy_Nxi,
     &                  numbdypoints,
     &                  chiextrap,xiextrap)

!----------------------------------------------------------------------

        implicit none
        integer numbdypoints
        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        integer bdy_Nchi,bdy_Nxi

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer e

        real*8 rhoextrap

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        !-----------------------------------------------------------------------

        chiextrap_min=minval(chiextrap)-1
        chiextrap_max=maxval(chiextrap)
        xiextrap_min=minval(xiextrap)-1
        xiextrap_max=maxval(xiextrap)

!calculate the number of different values of the chi and xi coordinates taken on extrapolated bdy points
        bdy_Nchi=0
        do while (chiextrap_min.lt.chiextrap_max)
            bdy_Nchi=bdy_Nchi+1
            chiextrap_min=
     &                minval(chiextrap,mask=chiextrap.gt.chiextrap_min)
        end do

        bdy_Nxi=0
        do while (xiextrap_min.lt.xiextrap_max)
            bdy_Nxi=bdy_Nxi+1
            xiextrap_min=minval(xiextrap,mask=xiextrap.gt.xiextrap_min)
        end do

        return
        end
!----------------------------------------------------------------------

c---------------------------------------------------------------------------------------------------------------------------------------------------------
c The following routine returns the arrays chibdy(bdy_Nchi) and xibdy(bdy_Nxi) containing the values of chi and xi at the extrapolated boundary points, in increasing order
c---------------------------------------------------------------------------------------------------------------------------------------------------------

        subroutine chibdy_xibdy(
     &                  chibdy,xibdy,
     &                  xpbdy,ypbdy,zpbdy,numbdypoints,
     &                  chiextrap,xiextrap,
     &                  bdy_Nchi,bdy_Nxi)

        implicit none
        integer numbdypoints
        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)

        real*8 chiextrap(numbdypoints)
        real*8 xiextrap(numbdypoints)

        integer bdy_Nchi,bdy_Nxi

        real*8 chibdy(bdy_Nchi)
        real*8 xibdy(bdy_Nxi)

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer i,j

        real*8 chiextrap_min,chiextrap_max
        real*8 xiextrap_min,xiextrap_max

        !-----------------------------------------------------------------------

        chiextrap_min=minval(chiextrap)-1
        chiextrap_max=maxval(chiextrap)
        xiextrap_min=minval(xiextrap)-1
        xiextrap_max=maxval(xiextrap)

        do i=1,bdy_Nchi
            chiextrap_min=
     &                minval(chiextrap,mask=chiextrap.gt.chiextrap_min)
            chibdy(i)=chiextrap_min
        end do

        do j=1,bdy_Nxi 
            xiextrap_min=minval(xiextrap,mask=xiextrap.gt.xiextrap_min)
            xibdy(j)=xiextrap_min
        end do

        return
        end
!----------------------------------------------------------------------
