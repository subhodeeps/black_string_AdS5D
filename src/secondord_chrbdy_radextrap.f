c----------------------------------------------------------------------
c routine to select chrbdy mask, which is not ex only at the outermost points
c that we use for second order radial extrapolation
c
c If we use derivatives to define near boundary quantities, we will only define 
c them at points between is and ie (js and je, ks and ke).  
c Therefore, for extrapolation, we can only select near boundary points whose  
c neighbors used for extrapolation are within that range
c We also need to make sure that those neighbours are not excised.
c The condition (chrbdy2(i+1,j,k).ne.ex) makes sure that (i,j,k) is the 
c outmost point satisfying the conditions of the previous for-loop, which 
c sets chrbdy2 as well as chrbdy. 
c In other words, if there's an outer point w.r.t. (i,j,k) that satisfies 
c those conditions, then we don't want to use (i,j,k) for extrapolation, 
c but we will use that other point. 
c----------------------------------------------------------------------

        subroutine secondord_chrbdy_radextrap(
     &                  chrbdy,
     &                  chrbdy2,
     &                  is,ie,js,je,ks,ke,
     &                  i,j,k,
     &                  xp1,yp1,zp1,rhop1,chip1,xip1,
     &                  chr,ex,Nx,Ny,Nz)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz),chrbdy2(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer ip2a,jp2a,kp2a
        integer ip2b,jp2b,kp2b
        integer ip2c,jp2c,kp2c
        integer ip2d,jp2d,kp2d

        integer ip3a,jp3a,kp3a
        integer ip3b,jp3b,kp3b
        integer ip3c,jp3c,kp3c
        integer ip3d,jp3d,kp3d

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1,rhop1,chip1,xip1

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

              if ((abs(xp1).ge.abs(yp1)).and.
     &            (abs(xp1).ge.abs(zp1))) then !(i.e., |xp1|>=|yp1|,|zp1|, so xp1 cannot be 0)
            
                if (xp1.gt.0) then !(i.e., we are in the upper part of 
                                    !the spatial grid, called "a" sector)

                  ip2a=i-1
                  ip2b=i-1
                  ip2c=i-1
                  ip2d=i-1

                  ip3a=i-2
                  ip3b=i-2
                  ip3c=i-2
                  ip3d=i-2

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((i+1).le.Nx) then 
                    if (chrbdy2(i+1,j,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      jp3a=j
                      jp3b=j
                      jp3c=j-1
                      jp3d=j-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((i-2).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                             ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((i-2).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jp2a=j
                      jp2b=j
                      jp2c=j+1
                      jp2d=j+1

                      jp3a=j
                      jp3b=j
                      jp3c=j+1
                      jp3d=j+1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((i-2).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((i-2).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kp2a=k
                      kp2b=k
                      kp2c=k-1
                      kp2d=k-1

                      kp3a=k
                      kp3b=k
                      kp3c=k-1
                      kp3d=k-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((i-2).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((i-2).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      kp3a=k
                      kp3b=k
                      kp3c=k+1
                      kp3d=k+1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((i-2).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((i-2).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                else !(i.e., xp1<0, i.e., we are in the lower part of 
                     !the spatial grid, called "b" sector)

                  ip2a=i+1
                  ip2b=i+1
                  ip2c=i+1
                  ip2d=i+1

                  ip3a=i+2
                  ip3b=i+2
                  ip3c=i+2
                  ip3d=i+2


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((i-1).ge.1) then 
                    if (chrbdy2(i-1,j,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      jp3a=j
                      jp3b=j
                      jp3c=j-1
                      jp3d=j-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((i+2).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((i+2).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jp2a=j
                      jp2b=j
                      jp2c=j+1
                      jp2d=j+1

                      jp3a=j
                      jp3b=j
                      jp3c=j+1
                      jp3d=j+1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((i+2).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((i+2).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kp2a=k
                      kp2b=k
                      kp2c=k-1
                      kp2d=k-1

                      kp3a=k
                      kp3b=k
                      kp3c=k-1
                      kp3d=k-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((i+2).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((i+2).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      kp3a=k
                      kp3b=k
                      kp3c=k+1
                      kp3d=k+1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((i+2).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((i+2).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                end if !closes condition on (xp1.gt.0)











              else if ((abs(yp1).ge.abs(zp1)).and.
     &                 (abs(yp1).ge.abs(xp1))) then !(i.e., |yp1|>=|zp1|,|xp1|, so yp1 cannot be 0)


                if (yp1.gt.0) then !(i.e., we are in either Q.Ia or Q.Ib or Q.IVa or Q.IVb)

                  jp2a=j-1
                  jp2b=j-1
                  jp2c=j-1
                  jp2d=j-1

                  jp3a=j-2
                  jp3b=j-2
                  jp3c=j-2
                  jp3d=j-2

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((j+1).le.Ny) then 
                    if (chrbdy2(i,j+1,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      kp2a=k
                      kp2b=k
                      kp2c=k-1
                      kp2d=k-1

                      kp3a=k
                      kp3b=k
                      kp3c=k-1
                      kp3d=k-1

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((j-2).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant Ib)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((j-2).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IVa or Q.IVb)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      kp3a=k
                      kp3b=k
                      kp3c=k+1
                      kp3d=k+1

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((j-2).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IVb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((j-2).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1| - xp1=zp1=0 is not an issue)

                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      ip3a=i
                      ip3b=i
                      ip3c=i-1
                      ip3d=i-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((j-2).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((j-2).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IVb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      ip3a=i
                      ip3b=i
                      ip3c=i+1
                      ip3d=i+1

                      if (zp1.gt.0) then !(i.e., quadrant Ib)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((j-2).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVb)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((j-2).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                else !(i.e., yp1<0, so we are in either Q.IIa or Q.IIb or Q.IIIa or Q.IIIb)

                  jp2a=j+1
                  jp2b=j+1
                  jp2c=j+1
                  jp2d=j+1

                  jp3a=j+2
                  jp3b=j+2
                  jp3c=j+2
                  jp3d=j+2


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((j-1).ge.1) then 
                    if (chrbdy2(i,j-1,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant IIa or IIb)

                      kp2a=k
                      kp2b=k
                      kp2c=k-1
                      kp2d=k-1

                      kp3a=k
                      kp3b=k
                      kp3c=k-1
                      kp3d=k-1

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((j+2).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((j+2).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IIIb)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      kp3a=k
                      kp3b=k
                      kp3c=k+1
                      kp3d=k+1

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((j+2).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((j+2).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1|, it's possible that xp1=zp1=0)

                    if (xp1.gt.0) then !(i.e., either quadrant IIa or IIIa)

                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      ip3a=i
                      ip3b=i
                      ip3c=i-1
                      ip3d=i-1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((j+2).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((j+2).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.IIb or Q.IIIb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      ip3a=i
                      ip3b=i
                      ip3c=i+1
                      ip3d=i+1

                      if (zp1.gt.0) then !(i.e., quadrant IIb)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((j+2).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIb)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((j+2).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                end if !closes condition on (yp1.gt.0)











              else !(i.e., |zp1|>=|xp1|,|yp1|, so zp1 cannot be 0)

                if (zp1.gt.0) then !(i.e., we are in either Q.Ia or Q.Ib or Q.IIa or Q.IIb)

                  kp2a=k-1
                  kp2b=k-1
                  kp2c=k-1
                  kp2d=k-1

                  kp3a=k-2
                  kp3b=k-2
                  kp3c=k-2
                  kp3d=k-2

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((k+1).le.ke) then 
                    if (chrbdy2(i,j,k+1).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      ip3a=i
                      ip3b=i
                      ip3c=i-1
                      ip3d=i-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks).or.
     &                      ((k-2).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks).or.
     &                      ((k-2).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IIb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      ip3a=i
                      ip3b=i
                      ip3c=i+1
                      ip3d=i+1

                      if (yp1.gt.0) then !(i.e., quadrant Ib)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks).or.
     &                      ((k-2).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIb)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks).or.
     &                      ((k-2).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1| - yp1=xp1=0 is not an issue)

                    if (yp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      jp3a=j
                      jp3b=j
                      jp3c=j-1
                      jp3d=j-1

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks).or.
     &                      ((k-2).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant Ib)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks).or.
     &                      ((k-2).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIb)

                      jp2a=j
                      jp2b=j
                      jp2c=j+1
                      jp2d=j+1

                      jp3a=j
                      jp3b=j
                      jp3c=j+1
                      jp3d=j+1

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks).or.
     &                      ((k-2).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks).or.
     &                      ((k-2).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))


                else !(i.e., zp1<0, so we are in either Q.IIIa or Q.IIIb or Q.IVa or Q.IVb)

                  kp2a=k+1
                  kp2b=k+1
                  kp2c=k+1
                  kp2d=k+1

                  kp3a=k+2
                  kp3b=k+2
                  kp3c=k+2
                  kp3d=k+2


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((k-1).ge.ks) then 
                    if (chrbdy2(i,j,k-1).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant IIIa or IVa)

                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      ip3a=i
                      ip3b=i
                      ip3c=i-1
                      ip3d=i-1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke).or.
     &                      ((k+2).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke).or.
     &                      ((k+2).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.IIIb or IVb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      ip3a=i
                      ip3b=i
                      ip3c=i+1
                      ip3d=i+1

                      if (yp1.gt.0) then !(i.e., quadrant IVb)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke).or.
     &                      ((k+2).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIb)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke).or.
     &                      ((k+2).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1|, it's possible that yp1=xp1=0)

                    if (yp1.gt.0) then !(i.e., either quadrant IVa or IVb)

                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      jp3a=j
                      jp3b=j
                      jp3c=j-1
                      jp3d=j-1

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke).or.
     &                      ((k+2).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IVb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke).or.
     &                      ((k+2).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIIa or IIIb)

                      jp2a=j
                      jp2b=j
                      jp2c=j+1
                      jp2d=j+1

                      jp3a=j
                      jp3b=j
                      jp3c=j+1
                      jp3d=j+1

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke).or.
     &                      ((k+2).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke).or.
     &                      ((k+2).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ip2a,jp2a,kp2a).eq.ex).or.
     &                           (chr(ip2b,jp2b,kp2b).eq.ex).or. 
     &                           (chr(ip2c,jp2c,kp2c).eq.ex).or. 
     &                           (chr(ip2d,jp2d,kp2d).eq.ex).or.
     &                           (chr(ip3a,jp3a,kp3a).eq.ex).or.
     &                           (chr(ip3b,jp3b,kp3b).eq.ex).or. 
     &                           (chr(ip3c,jp3c,kp3c).eq.ex).or. 
     &                           (chr(ip3d,jp3d,kp3d).eq.ex)
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))

                end if !closes condition on (zp1.gt.0)






              end if !closes condition on ((abs(xp1).gt.abs(yp1)).and.(abs(xp1).gt.abs(zp1)))

        return
        end


c------------------------------------------------------------------------------------------------------