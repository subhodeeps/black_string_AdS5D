c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the second order radially extrapolated value of the
c leadordcoeff_f grid function at the AdS boundary
c----------------------------------------------------------------------

        real*8 function secondord_func_radextrap(
     &                  leadordcoeff_f,
     &                  leadordcoeff_f_p1,
     &                  lind,
     &                  i,j,k,
     &                  xp1,yp1,zp1,
     &                  rhop1,chip1,xip1,
     &                  xex,yex,zex,
     &                  rhoex,chiex,xiex,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        integer i,j,k,is,ie,js,je,ks,ke,lind
        integer ip2a,jp2a,kp2a
        integer ip2b,jp2b,kp2b
        integer ip2c,jp2c,kp2c
        integer ip2d,jp2d,kp2d

        integer ip3a,jp3a,kp3a
        integer ip3b,jp3b,kp3b
        integer ip3c,jp3c,kp3c
        integer ip3d,jp3d,kp3d

        integer a,b,c,d

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz


        real*8 x0,y0,z0,rho0,q,chi0,xi0
        real*8 xp2a,xp2b,xp2c,xp2d
        real*8 yp2a,yp2b,yp2c,yp2d
        real*8 zp2a,zp2b,zp2c,zp2d

        real*8 xp3a,xp3b,xp3c,xp3d
        real*8 yp3a,yp3b,yp3c,yp3d
        real*8 zp3a,zp3b,zp3c,zp3d

        real*8 dx,dy,dz

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 leadordcoeff_f(Nx,Ny,Nz)
        real*8 leadordcoeff_f_p1

        real*8 leadordcoeff_f_p2a
        real*8 leadordcoeff_f_p2b
        real*8 leadordcoeff_f_p2c
        real*8 leadordcoeff_f_p2d

        real*8 leadordcoeff_f_p3a
        real*8 leadordcoeff_f_p3b
        real*8 leadordcoeff_f_p3c
        real*8 leadordcoeff_f_p3d

        real*8 leadordcoeff_f_p2
        real*8 leadordcoeff_f_p3
        real*8 leadordcoeff_f_p4

        real*8 xp1,yp1,zp1
        real*8 xp2,yp2,zp2
        real*8 xp3,yp3,zp3
        real*8 xp4,yp4,zp4
        real*8 xex,yex,zex,rhoex,chiex,xiex
        real*8 rhop1,chip1,xip1
        real*8 rhop2,chip2,xip2
        real*8 rhop3,chip3,xip3
        real*8 maxxyzp1

        real*8 bilinear_interp
        real*8 secondord_extrap
!----------------------------------------------------------------------

              chip2=chip1
              xip2=xip1

              chip3=chip1
              xip3=xip1


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

                  xp2a=x(ip2a)
                  xp2b=x(ip2b)
                  xp2c=x(ip2c)
                  xp2d=x(ip2d)
                  xp2=x(ip2a)

                  xp3a=x(ip3a)
                  xp3b=x(ip3b)
                  xp3c=x(ip3c)
                  xp3d=x(ip3d)
                  xp3=x(ip3a)

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

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      yp3a=y(jp3a)
                      yp3b=y(jp3b)
                      yp3c=y(jp3c)
                      yp3d=y(jp3d)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = -abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      yp3a=y(jp3a)
                      yp3b=y(jp3b)
                      yp3c=y(jp3c)
                      yp3d=y(jp3d)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = -abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = -abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = -abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      zp3a=z(kp3a)
                      zp3b=z(kp3b)
                      zp3c=z(kp3c)
                      zp3d=z(kp3d)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)


                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = -abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      zp3a=z(kp3a)
                      zp3b=z(kp3b)
                      zp3c=z(kp3c)
                      zp3d=z(kp3d)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = -abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., yp1<=0: quadrant IIIa)

                        if ((abs(yp1).lt.10.0d0**(-10)).and.
     &                      (abs(zp1).lt.10.0d0**(-10))) then
                          yp2a=0.0d0
                          yp2b=0.0d0
                          yp2c=0.0d0
                          yp2d=0.0d0

                          yp3a=0.0d0
                          yp3b=0.0d0
                          yp3c=0.0d0
                          yp3d=0.0d0

                          yp2=0.0d0
                          zp2=0.0d0
                          rhop2=abs(xp2)

                          yp3=0.0d0
                          zp3=0.0d0
                          rhop3=abs(xp3)

                          leadordcoeff_f_p2=
     &                        leadordcoeff_f(i-1,j,k)

                          leadordcoeff_f_p3=
     &                        leadordcoeff_f(i-2,j,k)


                        else

                          jp2a=j
                          jp2b=j+1
                          jp2c=j+1
                          jp2d=j

                          jp3a=j
                          jp3b=j+1
                          jp3c=j+1
                          jp3d=j

                          yp2a=y(jp2a)
                          yp2b=y(jp2b)
                          yp2c=y(jp2c)
                          yp2d=y(jp2d)

                          yp3a=y(jp3a)
                          yp3b=y(jp3b)
                          yp3c=y(jp3c)
                          yp3d=y(jp3d)


                          yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                          zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                          rhop2=abs(xp2/cos(PI*chip2))

                          yp3  = -abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                          zp3  = -abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                          rhop3=abs(xp3/cos(PI*chip3))

                          leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                          leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                          leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                          leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                          leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                          leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                          leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                          leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                          leadordcoeff_f_p2=
     &                        bilinear_interp(
     &                             leadordcoeff_f_p2a,
     &                             leadordcoeff_f_p2b,
     &                             leadordcoeff_f_p2c,
     &                             leadordcoeff_f_p2d,
     &                             yp2a,yp2b,yp2c,yp2d,
     &                             zp2a,zp2b,zp2c,zp2d,
     &                             yp2,
     &                             zp2)

                          leadordcoeff_f_p3=
     &                        bilinear_interp(
     &                             leadordcoeff_f_p3a,
     &                             leadordcoeff_f_p3b,
     &                             leadordcoeff_f_p3c,
     &                             leadordcoeff_f_p3d,
     &                             yp3a,yp3b,yp3c,yp3d,
     &                             zp3a,zp3b,zp3c,zp3d,
     &                             yp3,
     &                             zp3)

                        end if

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


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

                  xp2a=x(ip2a)
                  xp2b=x(ip2b)
                  xp2c=x(ip2c)
                  xp2d=x(ip2d)
                  xp2=x(ip2a)

                  xp3a=x(ip3a)
                  xp3b=x(ip3b)
                  xp3c=x(ip3c)
                  xp3d=x(ip3d)
                  xp3=x(ip3a)

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

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      yp3a=y(jp3a)
                      yp3b=y(jp3b)
                      yp3c=y(jp3c)
                      yp3d=y(jp3d)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)
                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = -abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      yp3a=y(jp3a)
                      yp3b=y(jp3b)
                      yp3c=y(jp3c)
                      yp3d=y(jp3d)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)
                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = -abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = -abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = -abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      zp3a=z(kp3a)
                      zp3b=z(kp3b)
                      zp3c=z(kp3c)
                      zp3d=z(kp3d)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)


                        yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = -abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      zp3a=z(kp3a)
                      zp3b=z(kp3b)
                      zp3c=z(kp3c)
                      zp3d=z(kp3d)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)


                        yp2  = abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                        zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                        rhop2=abs(xp2/cos(PI*chip2))

                        yp3  = abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                        zp3  = -abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                        rhop3=abs(xp3/cos(PI*chip3))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           yp2,
     &                           zp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           yp3,
     &                           zp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., yp1<=0: quadrant IIIa)

                        if ((abs(yp1).lt.10.0d0**(-10)).and.
     &                      (abs(zp1).lt.10.0d0**(-10))) then
                          yp2a=0.0d0
                          yp2b=0.0d0
                          yp2c=0.0d0
                          yp2d=0.0d0

                          yp3a=0.0d0
                          yp3b=0.0d0
                          yp3c=0.0d0
                          yp3d=0.0d0

                          yp2=0.0d0
                          zp2=0.0d0
                          rhop2=abs(xp2)

                          yp3=0.0d0
                          zp3=0.0d0
                          rhop3=abs(xp3)

                          leadordcoeff_f_p2=
     &                        leadordcoeff_f(i+1,j,k)

                          leadordcoeff_f_p3=
     &                        leadordcoeff_f(i+2,j,k)


                        else

                          jp2a=j
                          jp2b=j+1
                          jp2c=j+1
                          jp2d=j

                          jp3a=j
                          jp3b=j+1
                          jp3c=j+1
                          jp3d=j

                          yp2a=y(jp2a)
                          yp2b=y(jp2b)
                          yp2c=y(jp2c)
                          yp2d=y(jp2d)

                          yp3a=y(jp3a)
                          yp3b=y(jp3b)
                          yp3c=y(jp3c)
                          yp3d=y(jp3d)


                          yp2  = -abs(xp2*tan(PI*chip2)*cos(2*PI*xip2))
                          zp2  = -abs(xp2*tan(PI*chip2)*sin(2*PI*xip2))
                          rhop2=abs(xp2/cos(PI*chip2))

                          yp3  = -abs(xp3*tan(PI*chip3)*cos(2*PI*xip3))
                          zp3  = -abs(xp3*tan(PI*chip3)*sin(2*PI*xip3))
                          rhop3=abs(xp3/cos(PI*chip3))

                          leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                          leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                          leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                          leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                          leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                          leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                          leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                          leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                          leadordcoeff_f_p2=
     &                        bilinear_interp(
     &                             leadordcoeff_f_p2a,
     &                             leadordcoeff_f_p2b,
     &                             leadordcoeff_f_p2c,
     &                             leadordcoeff_f_p2d,
     &                             yp2a,yp2b,yp2c,yp2d,
     &                             zp2a,zp2b,zp2c,zp2d,
     &                             yp2,
     &                             zp2)

                          leadordcoeff_f_p3=
     &                        bilinear_interp(
     &                             leadordcoeff_f_p3a,
     &                             leadordcoeff_f_p3b,
     &                             leadordcoeff_f_p3c,
     &                             leadordcoeff_f_p3d,
     &                             yp3a,yp3b,yp3c,yp3d,
     &                             zp3a,zp3b,zp3c,zp3d,
     &                             yp3,
     &                             zp3)

                        end if

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


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

                  yp2a=y(jp2a)
                  yp2b=y(jp2b)
                  yp2c=y(jp2c)
                  yp2d=y(jp2d)
                  yp2=y(jp2a)

                  yp3a=y(jp3a)
                  yp3b=y(jp3b)
                  yp3c=y(jp3c)
                  yp3d=y(jp3d)
                  yp3=y(jp3a)

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

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      zp3a=z(kp3a)
                      zp3b=z(kp3b)
                      zp3c=z(kp3c)
                      zp3d=z(kp3d)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = abs(yp3*tan(2*PI*xip3))
                        xp3  = abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., xp1<=0: quadrant Ib)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = abs(yp3*tan(2*PI*xip3))
                        xp3  = -abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      end if

                    else !(i.e., zp1<=0, so Q.IVa or IVb)

                      kp2a=k
                      kp2b=k
                      kp2c=k+1
                      kp2d=k+1

                      kp3a=k
                      kp3b=k
                      kp3c=k+1
                      kp3d=k+1

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      zp3a=z(kp3a)
                      zp3b=z(kp3b)
                      zp3c=z(kp3c)
                      zp3d=z(kp3d)

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = -abs(yp3*tan(2*PI*xip3))
                        xp3  = abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., xp1<=0: quadrant IVb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = -abs(yp3*tan(2*PI*xip3))
                        xp3  = -abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      xp3a=x(ip3a)
                      xp3b=x(ip3b)
                      xp3c=x(ip3c)
                      xp3d=x(ip3d)

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = abs(yp3*tan(2*PI*xip3))
                        xp3  = abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., zp1<=0: quadrant IVa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = -abs(yp3*tan(2*PI*xip3))
                        xp3  = abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      xp3a=x(ip3a)
                      xp3b=x(ip3b)
                      xp3c=x(ip3c)
                      xp3d=x(ip3d)

                      if (zp1.gt.0) then !(i.e., quadrant Ib)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = abs(yp3*tan(2*PI*xip3))
                        xp3  = -abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., zp1<=0: quadrant IVb)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = -abs(yp3*tan(2*PI*xip3))
                        xp3  = -abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


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

                  yp2a=y(jp2a)
                  yp2b=y(jp2b)
                  yp2c=y(jp2c)
                  yp2d=y(jp2d)
                  yp2=y(jp2a)

                  yp3a=y(jp3a)
                  yp3b=y(jp3b)
                  yp3c=y(jp3c)
                  yp3d=y(jp3d)
                  yp3=y(jp3a)

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

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      zp3a=z(kp3a)
                      zp3b=z(kp3b)
                      zp3c=z(kp3c)
                      zp3d=z(kp3d)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)
                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = abs(yp3*tan(2*PI*xip3))
                        xp3  = abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., xp1<=0: quadrant IIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = abs(yp3*tan(2*PI*xip3))
                        xp3  = -abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      zp2a=z(kp2a)
                      zp2b=z(kp2b)
                      zp2c=z(kp2c)
                      zp2d=z(kp2d)

                      zp3a=z(kp3a)
                      zp3b=z(kp3b)
                      zp3c=z(kp3c)
                      zp3d=z(kp3d)

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)
                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = -abs(yp3*tan(2*PI*xip3))
                        xp3  = abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., xp1<=0: quadrant IIIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = -abs(yp3*tan(2*PI*xip3))
                        xp3  = -abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1| - the case xp1=zp1=0 is not an issue)

                    if (xp1.gt.0) then !(i.e., either quadrant IIa or IIIa)

                      ip2a=i
                      ip2b=i
                      ip2c=i-1
                      ip2d=i-1

                      ip3a=i
                      ip3b=i
                      ip3c=i-1
                      ip3d=i-1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      xp3a=x(ip3a)
                      xp3b=x(ip3b)
                      xp3c=x(ip3c)
                      xp3d=x(ip3d)

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = abs(yp3*tan(2*PI*xip3))
                        xp3  = abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = -abs(yp3*tan(2*PI*xip3))
                        xp3  = abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      end if

                    else !(i.e., xp1<=0, so Q.IIb or IIIb)

                      ip2a=i
                      ip2b=i
                      ip2c=i+1
                      ip2d=i+1

                      ip3a=i
                      ip3b=i
                      ip3c=i+1
                      ip3d=i+1

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      xp3a=x(ip3a)
                      xp3b=x(ip3b)
                      xp3c=x(ip3c)
                      xp3d=x(ip3d)

                      if (zp1.gt.0) then !(i.e., quadrant IIb)

                        kp2a=k
                        kp2b=k-1
                        kp2c=k-1
                        kp2d=k

                        kp3a=k
                        kp3b=k-1
                        kp3c=k-1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)


                        zp2  = abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = abs(yp3*tan(2*PI*xip3))
                        xp3  = -abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., zp1<=0: quadrant IIIb)

                        kp2a=k
                        kp2b=k+1
                        kp2c=k+1
                        kp2d=k

                        kp3a=k
                        kp3b=k+1
                        kp3c=k+1
                        kp3d=k

                        zp2a=z(kp2a)
                        zp2b=z(kp2b)
                        zp2c=z(kp2c)
                        zp2d=z(kp2d)

                        zp3a=z(kp3a)
                        zp3b=z(kp3b)
                        zp3c=z(kp3c)
                        zp3d=z(kp3d)


                        zp2  = -abs(yp2*tan(2*PI*xip2))
                        xp2  = -abs(yp2*(1/tan(PI*chip2))
     &                     /cos(2*PI*xip2))
                        rhop2=abs(yp2/(sin(PI*chip2)*cos(2*PI*xip2)))

                        zp3  = -abs(yp3*tan(2*PI*xip3))
                        xp3  = -abs(yp3*(1/tan(PI*chip3))
     &                     /cos(2*PI*xip3))
                        rhop3=abs(yp3/(sin(PI*chip3)*cos(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           zp2a,zp2b,zp2c,zp2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           zp2,
     &                           xp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           zp3a,zp3b,zp3c,zp3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           zp3,
     &                           xp3)


                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


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

                  zp2a=z(kp2a)
                  zp2b=z(kp2b)
                  zp2c=z(kp2c)
                  zp2d=z(kp2d)
                  zp2=z(kp2a)

                  zp3a=z(kp3a)
                  zp3b=z(kp3b)
                  zp3c=z(kp3c)
                  zp3d=z(kp3d)
                  zp3=z(kp3a)

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

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      xp3a=x(ip3a)
                      xp3b=x(ip3b)
                      xp3c=x(ip3c)
                      xp3d=x(ip3d)

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., yp1<=0: quadrant IIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = -abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      xp3a=x(ip3a)
                      xp3b=x(ip3b)
                      xp3c=x(ip3c)
                      xp3d=x(ip3d)

                      if (yp1.gt.0) then !(i.e., quadrant Ib)

                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = -abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., yp1<=0: quadrant IIb)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = -abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = -abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      yp3a=y(jp3a)
                      yp3b=y(jp3b)
                      yp3c=y(jp3c)
                      yp3d=y(jp3d)

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., xp1<=0: quadrant Ib)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = -abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      yp3a=y(jp3a)
                      yp3b=y(jp3b)
                      yp3c=y(jp3c)
                      yp3d=y(jp3d)

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = -abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., xp1<=0: quadrant IIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = -abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = -abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


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

                  zp2a=z(kp2a)
                  zp2b=z(kp2b)
                  zp2c=z(kp2c)
                  zp2d=z(kp2d)
                  zp2=z(kp2a)

                  zp3a=z(kp3a)
                  zp3b=z(kp3b)
                  zp3c=z(kp3c)
                  zp3d=z(kp3d)
                  zp3=z(kp3a)

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

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      xp3a=x(ip3a)
                      xp3b=x(ip3b)
                      xp3c=x(ip3c)
                      xp3d=x(ip3d)

                      if (yp1.gt.0) then !(i.e., quadrant IVa)
                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., yp1<=0: quadrant IIIa)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = -abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      xp2a=x(ip2a)
                      xp2b=x(ip2b)
                      xp2c=x(ip2c)
                      xp2d=x(ip2d)

                      xp3a=x(ip3a)
                      xp3b=x(ip3b)
                      xp3c=x(ip3c)
                      xp3d=x(ip3d)

                      if (yp1.gt.0) then !(i.e., quadrant IVb)
                        jp2a=j
                        jp2b=j-1
                        jp2c=j-1
                        jp2d=j

                        jp3a=j
                        jp3b=j-1
                        jp3c=j-1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = -abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      else !(i.e., yp1<=0: quadrant IIIb)

                        jp2a=j
                        jp2b=j+1
                        jp2c=j+1
                        jp2d=j

                        jp3a=j
                        jp3b=j+1
                        jp3c=j+1
                        jp3d=j

                        yp2a=y(jp2a)
                        yp2b=y(jp2b)
                        yp2c=y(jp2c)
                        yp2d=y(jp2d)

                        yp3a=y(jp3a)
                        yp3b=y(jp3b)
                        yp3c=y(jp3c)
                        yp3d=y(jp3d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = -abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = -abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1| - yp1=xp1=0 is not an issue)

                    if (yp1.gt.0) then !(i.e., either quadrant IVa or IVb)

                      jp2a=j
                      jp2b=j
                      jp2c=j-1
                      jp2d=j-1

                      jp3a=j
                      jp3b=j
                      jp3c=j-1
                      jp3d=j-1

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      yp3a=y(jp3a)
                      yp3b=y(jp3b)
                      yp3c=y(jp3c)
                      yp3d=y(jp3d)

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., xp1<=0: quadrant IVb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = -abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)

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

                      yp2a=y(jp2a)
                      yp2b=y(jp2b)
                      yp2c=y(jp2c)
                      yp2d=y(jp2d)

                      yp3a=y(jp3a)
                      yp3b=y(jp3b)
                      yp3c=y(jp3c)
                      yp3d=y(jp3d)

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ip2a=i
                        ip2b=i-1
                        ip2c=i-1
                        ip2d=i

                        ip3a=i
                        ip3b=i-1
                        ip3c=i-1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        xp2  = abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = -abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ip2a=i
                        ip2b=i+1
                        ip2c=i+1
                        ip2d=i

                        ip3a=i
                        ip3b=i+1
                        ip3c=i+1
                        ip3d=i

                        xp2a=x(ip2a)
                        xp2b=x(ip2b)
                        xp2c=x(ip2c)
                        xp2d=x(ip2d)

                        xp3a=x(ip3a)
                        xp3b=x(ip3b)
                        xp3c=x(ip3c)
                        xp3d=x(ip3d)

                        xp2  = -abs(zp2*(1/tan(PI*chip2))
     &                     /sin(2*PI*xip2))
                        yp2  = -abs(zp2/tan(2*PI*xip2))
                        rhop2=abs(zp2/(sin(PI*chip2)*sin(2*PI*xip2)))

                        xp3  = -abs(zp3*(1/tan(PI*chip3))
     &                     /sin(2*PI*xip3))
                        yp3  = -abs(zp3/tan(2*PI*xip3))
                        rhop3=abs(zp3/(sin(PI*chip3)*sin(2*PI*xip3)))

                        leadordcoeff_f_p2a=
     &                        leadordcoeff_f(ip2a,jp2a,kp2a)
                        leadordcoeff_f_p2b=
     &                        leadordcoeff_f(ip2b,jp2b,kp2b)
                        leadordcoeff_f_p2c=
     &                        leadordcoeff_f(ip2c,jp2c,kp2c)
                        leadordcoeff_f_p2d=
     &                        leadordcoeff_f(ip2d,jp2d,kp2d)

                        leadordcoeff_f_p3a=
     &                        leadordcoeff_f(ip3a,jp3a,kp3a)
                        leadordcoeff_f_p3b=
     &                        leadordcoeff_f(ip3b,jp3b,kp3b)
                        leadordcoeff_f_p3c=
     &                        leadordcoeff_f(ip3c,jp3c,kp3c)
                        leadordcoeff_f_p3d=
     &                        leadordcoeff_f(ip3d,jp3d,kp3d)

                        leadordcoeff_f_p2=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p2a,
     &                           leadordcoeff_f_p2b,
     &                           leadordcoeff_f_p2c,
     &                           leadordcoeff_f_p2d,
     &                           xp2a,xp2b,xp2c,xp2d,
     &                           yp2a,yp2b,yp2c,yp2d,
     &                           xp2,
     &                           yp2)

                        leadordcoeff_f_p3=
     &                      bilinear_interp(
     &                           leadordcoeff_f_p3a,
     &                           leadordcoeff_f_p3b,
     &                           leadordcoeff_f_p3c,
     &                           leadordcoeff_f_p3d,
     &                           xp3a,xp3b,xp3c,xp3d,
     &                           yp3a,yp3b,yp3c,yp3d,
     &                           xp3,
     &                           yp3)

                        secondord_func_radextrap=
     &                      secondord_extrap(
     &                           leadordcoeff_f_p1,
     &                           leadordcoeff_f_p2,
     &                           leadordcoeff_f_p3,
     &                           rhop1,rhop2,rhop3,rhoex)


                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))

                end if !closes condition on (zp1.gt.0)



              end if !closes condition on ((abs(xp1).gt.abs(yp1)).and.(abs(xp1).gt.abs(zp1)))

        return
        end
c--------------------------------------------------------------------------------------