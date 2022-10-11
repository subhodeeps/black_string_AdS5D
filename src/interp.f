c-----------------------------------------------------------------------
c Sets f at i,j via xy interpolation between interior point and AdS boundary
c-----------------------------------------------------------------------
!        subroutine interp_from_ads_bdy(f,x,y,z,L,i,j,k,chr,ex,Nx,Ny,Nz)
!        implicit none
!        integer Nx,Ny,Nz,i,j,k
!        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,x(Nx),y(Ny),z(Nz),L
!
!        real*8 dx,dy,dz,rho_bdy,xi
!        real*8 d_to_bdy
!        real*8 f_bdy,f_int
!
!        real*8 PI
!        parameter (PI=3.141592653589793d0)
!
!        ! initialize fixed-size variables
!        data dx,dy,rho_bdy,xi/0.0,0.0,0.0,0.0/
!        data d_to_bdy/0.0/
!        data f_bdy,f_int/0.0,0.0/
!
!        !--------------------------------------------------------------
!
!        dx=x(2)-x(1)
!        dy=y(2)-y(1)
!        dz=z(2)-z(1)
!
!        rho_bdy=1.0d0
!        f_bdy=0.0d0        
!
!        ! see if linear interpolation does the trick ... either in x or y
!
!        if (abs(x(i)).gt.abs(y(j))) then
!           xi=asin(y(j)/rho_bdy)
!           if (x(i).lt.0) then
!              xi=PI-xi
!              d_to_bdy=x(i)-rho_bdy*cos(xi)
!              if (i.eq.Nx) then
!                 write(*,*) 'interp_from_ads_bdy, error: out of bounds'
!                 write(*,*) 'i=Nx'
!                 stop
!              end if
!              f_int=f(i+1,j,k)
!           else
!              d_to_bdy=rho_bdy*cos(xi)-x(i)
!              if (i.eq.1) then
!                 write(*,*) 'interp_from_ads_bdy, error: out of bounds'
!                 write(*,*) 'i=1'
!                 stop
!              end if
!              f_int=f(i-1,j,k)
!           end if
!           f(i,j,k)=(dx*f_bdy+d_to_bdy*f_int)/(dx+d_to_bdy)
!        else
!           xi=acos(x(i)/rho_bdy)
!           d_to_bdy=rho_bdy*sin(xi)-y(j)
!           if (j.eq.1) then
!              write(*,*) 'interp_from_ads_bdy, error: out of bounds'
!                 write(*,*) 'j=1'
!              stop
!           end if
!           f_int=f(i,j-1,k)
!           f(i,j,k)=(dy*f_bdy+d_to_bdy*f_int)/(dy+d_to_bdy)
!        end if
!
!        return
!        end

!3+1 version!!!
c-----------------------------------------------------------------------
c Sets f at i,j via xy interpolation between interior point and AdS
c boundary
c-----------------------------------------------------------------------
        subroutine interp_from_ads_bdy(f,x,y,z,L,i,j,k,chr,ex,Nx,Ny,Nz)
        implicit none
        integer Nx,Ny,Nz,i,j,k
        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,x(Nx),y(Ny),z(Nz),L

        real*8 rhoijk
        real*8 dx,dy,dz,rho_bdy,xi,phi
        real*8 d_to_bdy
        real*8 f_bdy,f_int

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 maxxyz

        ! initialize fixed-size variables
        data dx,dy,dz,rho_bdy,xi,phi/0.0,0.0,0.0,0.0,0.0,0.0/
        data d_to_bdy/0.0/
        data f_bdy,f_int/0.0,0.0/

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        rhoijk=sqrt(x(i)**2+y(j)**2+z(k)**2)
        rho_bdy=1.0d0
        f_bdy=0.0d0

        maxxyz=max(abs(x(i)),abs(y(j)),abs(z(k)))

        ! see if linear interpolation does the trick ... either in x or
        ! y

        xi=acos(x(i)/rhoijk)

        if (abs(x(i)).eq.maxxyz) then
           if (x(i).lt.0) then
              d_to_bdy=x(i)-rho_bdy*cos(xi)
              if (i.eq.Nx) then
                 write(*,*) 'interp_from_ads_bdy, error: out of bounds'
                 write(*,*) 'i=Nx'
                 stop
              end if
              f_int=f(i+1,j,k)
           else
              d_to_bdy=rho_bdy*cos(xi)-x(i)
              if (i.eq.1) then
                 write(*,*) 'interp_from_ads_bdy, error: out of bounds'
                 write(*,*) 'i=1'
                 stop
              end if
              f_int=f(i-1,j,k)
           end if
           f(i,j,k)=(dx*f_bdy+d_to_bdy*f_int)/(dx+d_to_bdy)
!           write(*,*) "x(i),y(j),z(k),rho="
!     &                ,x(i),y(j),z(k),sqrt(x(i)**2+y(j)**2+z(k)**2)
!           write(*,*) "d_to_bdy,f_int,dx+d_to_bdy="
!     &                ,d_to_bdy,f_int,dx+d_to_bdy
        else if (abs(y(j)).eq.maxxyz) then

           if (y(j).ne.0) then
              if (z(k).lt.0) then
                 phi=atan2(z(k),y(j))+2*PI
              else
                 phi=atan2(z(k),y(j))
              end if

           else
              write(*,*) 'interp_from_ads_bdy, error: 
     &                     center of the grid'
              write(*,*) 'i,j,k,x(i),y(j),z(k)=',i,j,k,x(i),y(j),z(k)
              stop
           end if

           if (y(j).lt.0) then
              d_to_bdy=y(j)-rho_bdy*sin(xi)*cos(phi)
              if (j.eq.Ny) then
                 write(*,*) 'interp_from_ads_bdy, error: out of bounds'
                 write(*,*) 'j=Ny'
                 stop
              end if
              f_int=f(i,j+1,k)
           else
              d_to_bdy=rho_bdy*sin(xi)*cos(phi)-y(j)
              if (j.eq.1) then
                 write(*,*) 'interp_from_ads_bdy, error: out of bounds'
                 write(*,*) 'j=1'
                 stop
              end if
              f_int=f(i,j-1,k)
           end if
           f(i,j,k)=(dy*f_bdy+d_to_bdy*f_int)/(dy+d_to_bdy)

        else !i.e. abs(z(k)).eq.maxxyz

           if (z(k).ne.(0.0d0)) then

              if (z(k).lt.0) then
                 phi=atan2(z(k),y(j))+2*PI
              else
                 phi=atan2(z(k),y(j))
              end if

              if (z(k).lt.0) then
                 d_to_bdy=z(k)-rho_bdy*sin(xi)*sin(phi)
                 if (k.eq.Nz) then
                  write(*,*) 'interp_from_ads_bdy, error: out of bounds'
                  write(*,*) 'k=Nz'
                  stop
                 end if
                 f_int=f(i,j,k+1)
!           write(*,*) "x(i),y(j),z(k),rho="
!     &                ,x(i),y(j),z(k),sqrt(x(i)**2+y(j)**2+z(k)**2)
!           write(*,*) "d_to_bdy,f_int,dz+d_to_bdy="
!     &                ,d_to_bdy,f_int,dz+d_to_bdy
              else
                 d_to_bdy=rho_bdy*sin(xi)*sin(phi)-z(k)
                 if (k.eq.1) then
                  write(*,*) 'interp_from_ads_bdy, error: out of bounds'
                  write(*,*) 'k=1'
                  stop
                 end if
                 f_int=f(i,j,k-1)
              end if

           else
              write(*,*) 'interp_from_ads_bdy, error: 
     &                    center of the grid'
              write(*,*) 'i,j,k,x(i),y(j),z(k)=',i,j,k,x(i),y(j),z(k)
              stop
           end if

           f(i,j,k)=(dz*f_bdy+d_to_bdy*f_int)/(dz+d_to_bdy)
        end if

        return
        end


!THE FOLLOWING ROUTINES ARE NOT USED IN THE CODE FOR NOW
!c----------------------------------------------------------------------
!c Sets f at i,j via rho interpolation between interior point and AdS boundary 
!c----------------------------------------------------------------------
!        subroutine interp_from_ads_bdy_rho(f,x,y,z,L,i,j,k,chr,ex,
!     &         Nx,Ny,Nz,
!     &         shift,ghost_width)
!        implicit none
!        integer Nx,Ny,Nz,i,j,k,i_shift,j_shift
!        integer shift
!        real*8 rho_shift
!        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,f_x,x(Nx),y(Ny),z(Nz),L
!
!        real*8 x0,y0,z0,rho0,f_int,f_bdy,rhodist
!        real*8 d_x,d_y
!        real*8 dx,dy,dz,da,xi,q,PI,rhotobdy,rhotoint
!        parameter (PI=3.141592653589793d0)
!
!        integer ghost_width(6)
!
!        dx=x(2)-x(1)
!        dy=y(2)-y(1)
!        dz=z(2)-z(1)
!        x0=x(i)
!        y0=y(j)
!        z0=z(k)
!        rho0=sqrt(x0**2+y0**2)
!
!        rhotobdy=1.0d0-rho0
!        f_bdy=0.0d0
!
!        d_x=abs(x0)/rho0-abs(x0)
!        d_y=y0/rho0-y0
!
!        rho_shift=shift*dx
!
!        ! calculate i_shift,j_shift according to the  
!        ! distance rhodist of (x0,y0) from the interior
!        if (rho_shift.eq.0) then
!          i_shift=0
!          j_shift=0
!        else
!          rhodist=rho0-(1-rho_shift)
!        end if
!        i_shift=int((rhodist*abs(x0)/rho0)/dx)
!        j_shift=int((rhodist*y0/rho0)/dy)
!
!        ! correct i_shift for exceptional |x|>|y| points
!        if ( sqrt(x(i-1-i_shift)**2+y(j-j_shift)**2).gt.1-rho_shift
!     &       .and. x0.ge.0 .and. abs(x0).gt.abs(y0) ) then
!          i_shift=i_shift+1
!        end if
!        if ( sqrt(x(i+1+i_shift)**2+y(j-j_shift)**2).gt.1-rho_shift
!     &       .and. x0.lt.0 .and. abs(x0).gt.abs(y0) ) then
!          i_shift=i_shift+1
!        end if
!
!        ! correct j_shift for exceptional |x|<|y| points
!        if ( sqrt(x(i-i_shift)**2+y(j-1-j_shift)**2).gt.1-rho_shift
!     &       .and. x0.ge.0 .and. abs(x0).le.abs(y0) ) then
!          j_shift=j_shift+1
!        end if
!        if ( sqrt(x(i+i_shift)**2+y(j-1-j_shift)**2).gt.1-rho_shift
!     &       .and. x0.lt.0 .and. abs(x0).le.abs(y0) ) then
!          j_shift=j_shift+1
!        end if
!
!        ! implement interpolation
!        if (abs(x0).gt.abs(y0)) then 
!          da=d_y*(dx+i_shift*dx)/d_x-j_shift*dy
!          rhotoint=sqrt((dx+i_shift*dx)**2+(da+j_shift*dy)**2)
!          if (x0.ge.0) then
!            f_int=(da*f(i-1-i_shift,j-1-j_shift,k)
!     &            +(dy-da)*f(i-1-i_shift,j-j_shift,k))
!     &            /dy
!          else
!            f_int=(da*f(i+1+i_shift,j-1-j_shift,k)
!     &            +(dy-da)*f(i+1+i_shift,j-j_shift,k))
!     &            /dy
!          endif
!        else 
!          da=d_x*(dy+j_shift*dy)/d_y-i_shift*dx
!          rhotoint=sqrt((dy+j_shift*dy)**2+(da+i_shift*dx)**2)
!          if (x0.ge.0) then
!            f_int=(da*f(i-1-i_shift,j-1-j_shift,k)
!     &            +(dx-da)*f(i-i_shift,j-1-j_shift,k))
!     &            /dx
!          else
!            f_int=(da*f(i+1+i_shift,j-1-j_shift,k)
!     &            +(dx-da)*f(i+i_shift,j-1-j_shift,k))
!     &            /dx
!          endif
!        end if
!        f(i,j,k)=(rhotobdy*f_int+rhotoint*f_bdy)/(rhotobdy+rhotoint)
!
!        return
!        end
!
!
!c----------------------------------------------------------------------
!c Sets f at i,j via rho interpolation between two interior pts and AdS boundary 
!c with Neville's algorithm for Lagrange interpolating polynomials defined
!c using three pts 
!c (x1,y1=f_int1)
!c (x2,y2=f_int2)
!c (x3,y3=f_bdy)
!c and polynomials going through these points
!c p12(x0) = ( (x2-x0)*y1 + (x0-x1)*y2 ) / (x2 - x1)
!c p23(x0) = ( (x3-x0)*y2 + (x0-x2)*y3 ) / (x3 - x2)
!c p13(x0) = ( (x3-x0)*p12(x0) + (x0-x1)*p23(x0) ) / (x3 - x1)
!c where the desired interpolated value is f=p13(x0)
!c----------------------------------------------------------------------
!        subroutine interp_from_ads_bdy_rho_3pt(f,x,y,z,L,i,j,k,chr,
!     &          ex,Nx,Ny,Nz,
!     &         rho_shift,ghost_width)
!        implicit none
!        integer Nx,Ny,Nz,i,j,k
!        integer i_shift
!
!        real*8 rho_shift
!
!        integer i_shift1,j_shift1
!        real*8 rhodist1,da1
!
!        integer i_shift2,j_shift2
!        real*8 rhodist2,da2
!
!        real*8 f(Nx,Ny,Nz),chr(Nx,Ny,Nz),ex,f_x,x(Nx),y(Ny),z(Nz),L
!
!        real*8 x0,y0,z0,rho0,rho0_bdy,f_int1,f_int2,f_bdy
!        real*8 x20,x01,x30,x02,x21,x32,x31,y1,y2,y3,p12,p23
!        real*8 dx,dy,dz,d_x,d_y,xi,q,PI,rho1,rho2,rho_bdy
!        parameter (PI=3.141592653589793d0)
!
!        integer ghost_width(6)
!
!        dx=x(2)-x(1)
!        dy=y(2)-y(1)
!        dz=z(2)-z(1)
!        rho0_bdy=1
!        x0=x(i)
!        y0=y(j)
!        z0=z(k)
!        rho0=sqrt(x0**2+y0**2)
!        d_x=rho0_bdy*abs(x0)/rho0-abs(x0)
!        d_y=rho0_bdy*y0/rho0-y0
!        rho_bdy=sqrt(d_x**2+d_y**2)
!        f_bdy=0
!
!        ! calculate i_shift1,j_shift1 according to the distance rhodist1 of (x0,y0)  
!        ! from the dx into the interior
!        if (rho_shift.eq.0) then
!          i_shift1=0
!          j_shift1=0
!        else
!          rhodist1=rho0-(1-rho_shift-dx)
!        end if
!        i_shift1=int((rhodist1*abs(x0)/rho0)/dx)
!        j_shift1=int((rhodist1*y0/rho0)/dy)
!
!        ! correct i_shift1 for exceptional |x|>|y| points
!        if ( sqrt(x(i-1-i_shift1)**2+y(j-j_shift1)**2).gt.1-rho_shift-dx
!     &       .and. x0.ge.0 .and. abs(x0).gt.abs(y0) ) then
!          i_shift1=i_shift1+1
!        end if
!        if ( sqrt(x(i+1+i_shift1)**2+y(j-j_shift1)**2).gt.1-rho_shift-dx
!     &       .and. x0.lt.0 .and. abs(x0).gt.abs(y0) ) then
!          i_shift1=i_shift1+1
!        end if
!
!        ! correct j_shift1 for exceptional |x|<|y| points
!        if ( sqrt(x(i-i_shift1)**2+y(j-1-j_shift1)**2).gt.1-rho_shift-dx
!     &       .and. x0.ge.0 .and. abs(x0).le.abs(y0) ) then
!          j_shift1=j_shift1+1
!        end if
!        if ( sqrt(x(i+i_shift1)**2+y(j-1-j_shift1)**2).gt.1-rho_shift-dx
!     &       .and. x0.lt.0 .and. abs(x0).le.abs(y0) ) then
!          j_shift1=j_shift1+1
!        end if
!
!        ! calculate i_shift2,j_shift2 according to the distance rhodist2 
!        ! of (x0,y0) from the interior
!        if (rho_shift.eq.0) then
!          i_shift2=0
!          j_shift2=0
!        else
!          rhodist2=rho0-(1-rho_shift)
!        end if
!        i_shift2=int((rhodist2*abs(x0)/rho0)/dx)
!        j_shift2=int((rhodist2*y0/rho0)/dy)
!
!        ! correct i_shift2 for exceptional |x|>|y| points
!        if ( sqrt(x(i-1-i_shift2)**2+y(j-j_shift2)**2).gt.1-rho_shift
!     &       .and. x0.ge.0 .and. abs(x0).gt.abs(y0) ) then
!          i_shift2=i_shift2+1
!        end if
!        if ( sqrt(x(i+1+i_shift2)**2+y(j-j_shift2)**2).gt.1-rho_shift
!     &       .and. x0.lt.0 .and. abs(x0).gt.abs(y0) ) then
!          i_shift2=i_shift2+1
!        end if
!
!        ! correct j_shift2 for exceptional |x|<|y| points
!        if ( sqrt(x(i-i_shift2)**2+y(j-1-j_shift2)**2).gt.1-rho_shift
!     &       .and. x0.ge.0 .and. abs(x0).le.abs(y0) ) then
!          j_shift2=j_shift2+1
!        end if
!        if ( sqrt(x(i+i_shift2)**2+y(j-1-j_shift2)**2).gt.1-rho_shift
!     &       .and. x0.lt.0 .and. abs(x0).le.abs(y0) ) then
!          j_shift2=j_shift2+1
!        end if
!
!        if (abs(x0).gt.abs(y0) .and. i_shift1.eq.i_shift2) 
!     &    i_shift1=i_shift1+1
!        if (abs(x0).le.abs(y0) .and. j_shift1.eq.j_shift2) 
!     &    j_shift1=j_shift1+1
!
!        ! implement interpolation
!        if (abs(x0).gt.abs(y0)) then 
!          da1=d_y*(dx+i_shift1*dx)/d_x-j_shift1*dy
!          da2=d_y*(dx+i_shift2*dx)/d_x-j_shift2*dy
!          rho1=sqrt((dx+i_shift1*dx)**2+(da1+j_shift1*dy)**2)
!          rho2=sqrt((dx+i_shift2*dx)**2+(da2+j_shift2*dy)**2)
!          if (x0.ge.0) then
!            f_int1=(da1*f(i-1-i_shift1,j-1-j_shift1,k)
!     &            +(dy-da1)*f(i-1-i_shift1,j-j_shift1,k))
!     &            /dy
!            f_int2=(da2*f(i-1-i_shift2,j-1-j_shift2,k)
!     &            +(dy-da2)*f(i-1-i_shift2,j-j_shift2,k))
!     &            /dy
!          else
!            f_int1=(da1*f(i+1+i_shift1,j-1-j_shift1,k)
!     &            +(dy-da1)*f(i+1+i_shift1,j-j_shift1,k))
!     &            /dy
!            f_int2=(da2*f(i+1+i_shift2,j-1-j_shift2,k)
!     &            +(dy-da2)*f(i+1+i_shift2,j-j_shift2,k))
!     &            /dy
!          endif
!        else 
!          da1=d_x*(dy+j_shift1*dy)/d_y-i_shift1*dx
!          da2=d_x*(dy+j_shift2*dy)/d_y-i_shift2*dx
!          rho1=sqrt((dy+j_shift1*dy)**2+(da1+i_shift1*dx)**2)
!          rho2=sqrt((dy+j_shift2*dy)**2+(da2+i_shift2*dx)**2)
!          if (x0.ge.0) then
!            f_int1=(da1*f(i-1-i_shift1,j-1-j_shift1,k)
!     &            +(dx-da1)*f(i-i_shift1,j-1-j_shift1,k))
!     &            /dx
!            f_int2=(da2*f(i-1-i_shift2,j-1-j_shift2,k)
!     &            +(dx-da2)*f(i-i_shift2,j-1-j_shift2,k))
!     &            /dx
!          else
!            f_int1=(da1*f(i+1+i_shift1,j-1-j_shift1,k)
!     &            +(dx-da1)*f(i+i_shift1,j-1-j_shift1,k))
!     &            /dx
!            f_int2=(da2*f(i+1+i_shift2,j-1-j_shift2,k)
!     &            +(dx-da2)*f(i+i_shift2,j-1-j_shift2,k))
!     &            /dx
!          endif
!        end if
!        x02=rho2
!        x01=rho1
!        x30=rho_bdy
!        x20=-x02
!        x21=x01-x02
!        x32=x30+x02
!        x31=x30+x01
!        y1=f_int1
!        y2=f_int2
!        y3=f_bdy
!        p12=(x20*y1+x01*y2)/x21
!        p23=(x30*y2+x02*y3)/x32
!        f(i,j,k)=(x30*p12+x01*p23)/x31
!
!        if (rho1.eq.rho2) then
!          write(*,*) '--------------------------------'
!          write(*,*) 'x,y=',x0,y0
!          write(*,*) 'i_shift1,j_shift1=',i_shift1,j_shift1
!          write(*,*) 'i_shift2,j_shift2=',i_shift2,j_shift2
!          write(*,*) '--------------------------------'
!        end if
!
!        return
!        end
