C-----------------------------------------------------------------------
C  nek5000 user-file template
C
C  user specified routines:
C     - userbc : boundary conditions
C     - useric : initial conditions
C     - uservp : variable properties
C     - userf  : local acceleration term for fluid
C     - userq  : local source term for scalars
C     - userchk: general purpose routine for checking errors etc. 
C
C-----------------------------------------------------------------------

      subroutine uservp(ix,iy,iz,eg) ! set variable properties
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      udiff  = 0.0
      utrans = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      ffx = 0.0
      ffy = 0.0
      ffz = 0.0
      return
      end

c-----------------------------------------------------------------------

      subroutine userq(ix,iy,iz,eg) ! set source term
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      qvol   = -2.0*uz
      source = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userbc(ix,iy,iz,iside,ieg) 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      ux = 0.0
      uy = 0.0
      uz = 0.0
      temp=0.0
      
      return
      end

c-----------------------------------------------------------------------

      subroutine useric(ix,iy,iz,ieg) ! set up initial conditions

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'NEKUSE_DEF'
      include 'NEKUSE'
      include 'PARALLEL_DEF'
      include 'PARALLEL'

      integer ix,iy,iz,ieg
      real AR_ic, H_ic, w_ic, a_ic, L_ic, ys_ic, xs_ic
      real alpha1_ic,alpha2_ic, wso_ic, ws1_ic, ws2_ic, wstotal_ic
      real pi, coeff_ic

!     Define pi
      pi= 4.*atan(1.0)

!     Duct aspect ratio
      AR_ic=1.

!     Coefficient for Ub=1
      coeff_ic=7.15150691475340

!     Geometry
      H_ic=2.
      w_ic=AR_ic*H_ic
      a_ic=1./AR_ic
      L_ic=w_ic/2.
      ys_ic=y/L_ic
      xs_ic=x/L_ic

!     Laminar duct profile from Panton
!     Leading order term
      wso_ic=1./2.*(a_ic**2.-ys_ic**2.)

!     Parameter alpha for first and second order terms
      alpha1_ic=(2.*1.-1.)*pi/(2.*a_ic)
      alpha2_ic=(2.*2.-1.)*pi/(2.*a_ic)

!     First and second order terms
      ws1_ic=2./a_ic*(-1.)/(alpha1_ic)**3.*cos(alpha1_ic*ys_ic)
     &     *(exp(alpha1_ic*(xs_ic-1.))+exp(-alpha1_ic*(xs_ic+1.)))
     &     /(1+exp(-2.*alpha1_ic))
      ws2_ic=2./a_ic*(1.)/(alpha2_ic)**3.*cos(alpha2_ic*ys_ic)
     &     *(exp(alpha2_ic*(xs_ic-1.))+exp(-alpha2_ic*(xs_ic+1.)))
     &     /(1+exp(-2.*alpha2_ic))

!     Sum of leading, first and second order terms
      wstotal_ic=wso_ic+ws1_ic+ws2_ic

!     Initial conditions: laminar profile
      ux=0.
      uy=0.
      uz=coeff_ic*wstotal_ic

      temp= -0.71/0.066*uz ! -Pr *  u_z / u_t
      
      return
      end

c-----------------------------------------------------------------------

      subroutine userchk()

      implicit none

      include 'SIZE_DEF'      
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'                    ! xm1, ym1, zm1
      include 'SOLN_DEF'
      include 'SOLN'                    ! T
      include 'MASS_DEF'
      include 'MASS'                    !BM1 for lambda2
      include 'TSTEP_DEF'
      include 'TSTEP'                   ! ISTEP
      include 'INPUT_DEF'
      include 'INPUT'                   ! PARAM(12) (DT)
      include 'USERPAR'         ! l2freq, FIXGEOM, NEW_DT
      include 'TRIPF'
      
c      real temp_scaled(lx1*ly1*lz1*lelt)
c      integer ntot,j

c      ntot=lx1*ly1*lz1*lelt
      
      call checkpoint                   ! Restart check


         call stat_avg_2D
      
      return
      end
      
c-----------------------------------------------------------------------
      
      real function step(x)
c
c     Smooth step function:
c     x<=0 : step(x) = 0
c     x>=1 : step(x) = 1
c     Non-continuous derivatives at x=0.02 and x=0.98
c     
      implicit none

      real x

      if (x.le.0.02) then
         step = 0.0
      else
         if (x.le.0.98) then
            step = 1./( 1. + exp(1./(x - 1.) + 1./x) )
         else
            step = 1.
         end if
      end if

      end function step
      
c-----------------------------------------------------------------------

      subroutine readwallfile



      include 'SIZE_DEF'
      include 'SIZE'
      include 'TRIPF'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
            
      integer len,ierr
      character*132 wallfname
      character*1 wallfnam1(132)
      equivalence (wallfnam1,wallfname)

      ierr = 0
      if (nid .eq. 0) then
         call blank(wallfname,132)
         len = ltrunc(SESSION,132)
         call chcopy(wallfnam1(1),SESSION,len)
         call chcopy(wallfnam1(len+1),'.wall',5)
         open(unit=75,file=wallfname,err=30,status='old')
         read(75,*,err=30)
         read(75,*,err=30) nwalls
         read(75,*,err=30) nwallpar
         read(75,*,err=30) npwallpar
         goto 31
 30      ierr=1
 31      continue
      endif
      call err_chk(ierr,'Error reading .wall file.$')
      call bcast(nwalls, ISIZE)
      call bcast(nwallpar, ISIZE)
      call bcast(npwallpar, ISIZE)

      if(nwalls .gt. maxwalls .or. nwallpar .gt. maxwallpar
     $           .or. npwallpar .gt. maxpwallpar ) then
         if(nid .eq. 0) then
           write(6,*) 'Too many walls/parameters in ',wallfname
         endif
         call exitt
      endif
      if (nid .eq. 0) then
c   read global parameters
        read(75,*,err=32)
        do i=1,nwallpar
          read(75,*,err=32) wallpar(i)
        end do
c   read wall definitions and parameters
        read(75,*,err=32)
        read(75,*,err=32)
        do i=1,nwalls
          read(75,*,err=32) direction(i)
          read(75,*,err=32) tripx(i),tripy(i),tripz(i)
          do j=1,npwallpar
            read(75,*,err=32) (pwallpar(k,j,i), k=1,3)
          end do
        end do
        goto 33
 32     ierr=1
 33     continue
      endif
      call err_chk(ierr,'Not enough walls.$')
      call bcast(wallpar,nwallpar*WDSIZE)
      call bcast(direction,nwalls*CSIZE)
      call bcast(tripx,nwalls*WDSIZE)
      call bcast(tripy,nwalls*WDSIZE)
      call bcast(tripz,nwalls*WDSIZE)
      call bcast(pwallpar(1,1,1),3*npwallpar*nwalls*WDSIZE)

c      if (nid.eq.0) write(*,*) 'Directions',direction
     
      return
      end
c----------------------------------------------------------------------
      subroutine znekgen(wall)

      implicit none
      
c      include 'SIZE_DEF'
c      include 'SIZE'
c      include 'TOTAL'
c      include 'TRIPF'

      include 'SIZE_DEF'
      include 'SIZE'
      include 'GEOM_DEF'
      include 'GEOM'            ! xm1, ym1, zm1
      include 'SOLN_DEF'
      include 'SOLN'            ! T
      include 'MASS_DEF'
      include 'MASS'            !BM1 for lambda2
      include 'TSTEP_DEF'
      include 'TSTEP'           ! ISTEP
      include 'INPUT_DEF'
      include 'INPUT'           ! PARAM(12) (DT)
      include 'USERPAR'         ! l2freq, FIXGEOM, NEW_DT
      include 'TRIPF'
      include 'PARALLEL_DEF'
      include 'PARALLEL'
            
      real dx1(lx1,ly1,lz1,lelv), bouxm1(lx1,ly1,lz1,lelv)
      real dx2(lx1,ly1,lz1,lelv), bouxm2(lx1,ly1,lz1,lelv)
      real dr2(lx1,ly1,lz1,lelv), bouxm3(lx1,ly1,lz1,lelv)
      real tripx1, tripx2

      real vals(maxlxyz,nelv), valsSort(maxlxyz,nelv)
      real valsf(nelv), valsfSort(nelv)
      real gvalsf(np*lelv),gvalsfw(np*lelv),gvalsfSort(np*lelv)
      integer valsfw(nelv), gvalsfi(np*lelv), wall
      real gCloseGLL(2), lCloseGLL, realTemp
      real znekw(lelv*maxlxyz)
      integer lCloseGLLid, cGLLnid, intTemp
      integer cvals, cvals1, cvals2, myit, itx, ity, itz
      integer i, j, ix, iy, iz, k

      if (nid.eq.0) write(*,*) 'Wall:',wall,direction(wall)
      
      itx = lx1
      ity = ly1
      itz = lz1
c   compute the differences
      if (direction(wall) .eq. 'x') then
         bouxm1 = ym1
         bouxm2 = zm1
         bouxm3 = xm1
         myit = itx
         tripx1 = tripy(wall)
         tripx2 = tripz(wall)
         itx = 1
c         write(*,*) 'wallx',tripx1,tripx2,itx,nid
      elseif (direction(wall) .eq. 'y') then
         bouxm1 = xm1
         bouxm2 = zm1
         bouxm3 = ym1
         myit = ity
         tripx1 = tripx(wall)
         tripx2 = tripz(wall)
         ity = 1
c         write(*,*) 'wally',tripx1,tripx2,ity,nid
      else
         bouxm1 = xm1
         bouxm2 = ym1
         bouxm3 = zm1
         tripx1 = tripx(wall)
         tripx2 = tripy(wall)
         myit = itz
         itz = 1
      endif

      dx1 = tripx1 - bouxm1
      dx2 = tripx2 - bouxm2
      dr2 = dx1*dx1 + dx2*dx2
      lCloseGLL = dr2(1,1,1,1)
      lCloseGLLid = 1

c      if (nid.eq.0) write(*,*) 'LCLOSE',lCloseGLL
      
c   calculate the local minimum distance
      do j = 1, nelv
      do iz = 1,itz
        do iy = 1,ity
          do ix = 1,itx
            if (dr2(ix,iy,iz,j) .lt. lCloseGLL) then
              lCloseGLL = dr2(ix,iy,iz,j)
              lCloseGLLid = ix + lx1*(iy-1) + lx1*ly1*(iz-1)
     $                         + (j-1)*lx1*ly1*lz1
            end if
          end do
        end do
      end do
      end do
      gCloseGLL(1) = lCloseGLL

c      if (nid.eq.0) write(*,*) 'GCLOSE',gCloseGLL

c      print *, 'ISM', nid, lCloseGLL
c   pick the global minimum distance
      call gop(gCloseGLL(1),realTemp,'m  ',1)

c      write(*,*) 'Global', gCloseGLL
      
c   chose a proc who has this distance
      if (lCloseGLL .eq. gCloseGLL(1)) then
         cGLLnid = nid
      else
         cGLLnid = 0
      end if
      call igop(cGLLnid,intTemp,'M  ',1)

c   share its x,y value to everyone
      if (cGLLnid .eq. nid) then
         gCloseGLL(1) = bouxm1(lCloseGLLid,1,1,1)
         gCloseGLL(2) = bouxm2(lCloseGLLid,1,1,1)
      else
         gCloseGLL(1) = 0.
         gCloseGLL(2) = 0.
      end if
c      print *, 'ISM', nid, gCloseGLL(1), gCloseGLL(2),tripx1,tripx2
      call bcastn0(gCloseGLL(1),2*WDSIZE,cGLLnid)

      write(*,*) '(gCloseGLL(1)',gCloseGLL(1),nid
      
c      print *, 'ISM', nid, gCloseGLL(1), gCloseGLL(2),tripx1,tripx2
c   sort the first z-value of each element containing the tripping points
      cvals = 0
      do j = 1,nelv
       do iz = 1, itz
         do iy = 1, ity
           do ix = 1, itx
             if (bouxm1(ix,iy,iz,j) .eq. gCloseGLL(1)
     $         .and. bouxm2(ix,iy,iz,j) .eq. gCloseGLL(2)) then
               cvals = cvals + 1
               if (direction(wall) .eq. 'x') then
c                  vals(1:myit,cvals) = bouxm3(:,iy,iz,j)
                  do k=1,myit
                     vals(k,cvals) = bouxm3(k,iy,iz,j)
                  enddo
               elseif (direction(wall) .eq. 'y') then
c                  vals(1:myit,cvals) = bouxm3(ix,:,iz,j)
                  do k=1,myit
                     vals(k,cvals) = bouxm3(ix,k,iz,j)
                  enddo
               else
c                  vals(1:myit,cvals) = bouxm3(ix,iy,:,j)
                  do k=1,myit
                     vals(k,cvals) = bouxm3(ix,iy,k,j)
                  enddo
               end if
               valsf(cvals) = bouxm3(ix,iy,iz,j)
               goto 100
             end if
           end do
         end do
       end do
 100   continue
      end do
      call sorts(valsfSort,valsf,valsfw,cvals)

c      print *, 'ISM', nid, cvals, valsfSort(1:cvals)
c   remove duplicate and share with everyone
      gvalsf = huge(1.0)
      if (cvals .gt. 0) then
         cvals1 = 1
c         valsSort(:,1) = vals(:,valsfw(1))
         do k=1,maxlxyz
            valsSort(k,1) = vals(k,valsfw(1))
         enddo
         gvalsf(1 + nid*lelv) = valsfSort(1)
      else
         cvals1 = 0
      end if
      do i = 2, cvals
         if(valsfSort(i) .ne. valsfSort(i-1)) then
           cvals1 = cvals1 + 1
c           valsSort(:,cvals1) = vals(:,valsfw(i))
           do k=1,maxlxyz
              valsSort(k,cvals1) = vals(k,valsfw(i))
           enddo
           gvalsf(cvals1 + nid*lelv) = valsfSort(i)
         end if
      end do
      call gop(gvalsf, gvalsfw,'m  ', np*lelv)

c      kpts(wall)=96
c      nnelx1x2(wall)=64
      
c   define kpts (lx * number of direction elements), nnelx1x2 (nx1*nx2), znek
      call sorts(gvalsfSort,gvalsf,gvalsfi,np*lelv)
c      print *, 'ISM1', nid, cvals, gvalsfSort(1:10)
      cvals2 = 1
      do i = 1,np*lelv
        if (gvalsfSort(i) .ne. gvalsfSort(cvals2)) then
           cvals2 = cvals2 + 1
           gvalsfSort(cvals2) = gvalsfSort(i)
        endif
        if (i .ne. cvals2) gvalsfSort(i) = huge(1.0)
c        gvalsfSort(i) = huge(1.0)
      end do
      
c      print *, 'ISM2', nid, cvals, gvalsfSort(1:10)

c      znek(:,wall) = huge(1.0)

      do k=1,lelv*maxlxyz
         znek(k,wall) = huge(1.0)
      enddo
      
      cvals2 = 1
      do i = 1, lelv
        if (gvalsfSort(i) .eq. huge(1.0)) then
c          print *, 'ISM fini', i, lz1
          kpts(wall) = (i-1)*myit
          nnelx1x2(wall) = nelgv/(i-1)
c          kpts(wall)=lelx*lx1 ! lelx*lx1
c          write(*,*) 'kpts',kpts(wall),nid,wall,nnelx1x2(wall),
c     &         gvalsfSort(i)
c          
          exit
        end if
        if (gvalsf(cvals2 + nid*lelv) .eq. gvalsfSort(i)) then
          do j = 1,myit !i*lz1,(i+1)*lz1
             znek((i-1)*myit+j,wall) = valsSort(j,cvals2)
c             write(*,*) 'znek_all',znek(1,wall),nid,wall
          end do
          cvals2 = cvals2 + 1
        end if
      end do
      call gop(znek(1,wall), znekw,'m  ', kpts(wall))

c      write(*,*) 'znek_after',znek(1,wall),nid,wall

      write(*,*) 'kpts',kpts(wall),nid,wall,nnelx1x2(wall),
     &     gvalsfSort(i)
      
      
      if (nid.eq.0) write(*,*) 'Test',znek(1,wall)
      if (nid.eq.0) write(*,*) 'Test2',kpts(wall)

      if (nid .eq. 0) then
      do i=1,kpts(wall)
         print *,'ISM', znek(i,wall)
      end do
      end if
c   Done!

ccccccccccccccc
c   print the values for script processing (no more needed)

c      if (nid .eq. 0) write (6,"('ISM2',x,i7)") nelgv
c      if (nid .eq. 0) write (6,"('ISM3',x,i7)") lz1
c      write(clz1,"(i2)") lz1
c      do k = 1,nelv
c        do i = 1, lx1*ly1
c          if (xm1(i,1,1,k) .eq. gCloseGLL(1)
c     $       .and. ym1(i,1,1,k) .eq. gCloseGLL(2)) then
c           write (6,"('ISM1'," // adjustr(clz1) // "(x,g25.16E4))")
c     $            (zm1(i,1,j,k),j=1,lz1)
c          end if
c       end do
c      end do
c      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine bcastn0(buf,len,proc)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)

      call mpi_bcast (buf,len,mpi_byte,proc,nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine tripf

      include 'SIZE_DEF'
      include 'SIZE'
      include 'TSTEP_DEF'
      include 'TSTEP'
      include 'TRIPF'
      
      integer k

      integer z,i
      real p,b
      real tamps, tampt, tdt
      integer num_modes

      tamps = wallpar(1)
      tampt = wallpar(2)
      tdt   = wallpar(3)
      num_modes = int(wallpar(4))

c
c     Generate the time independent part fzt(z,2)
c     at first iteration

      if (istep.eq.1) then
c
c     Get random distribution and rescale
c
         do k=1,nwalls
c            if (k.le.4) then
c               call rand_func(fzt2(1,k),znek(1,k),kpts(k),seed1,
c     &              num_modes)
c            elseif (k.gt.4.and.k.le.8) then
c               call rand_func(fzt2(1,k),znek(1,k),kpts(k),seed2,
c     &              num_modes)
c            else
c               call rand_func(fzt2(1,k),znek(1,k),kpts(k),seed3,
c     &              num_modes)
c            endif
            call rand_func(fzt2(1,k),znek(1,k),kpts(k),seed,num_modes)
            do z=1,kpts(k)
               fzt2(z,k)=tamps*fzt2(z,k)
            end do
         enddo
         ntdt=-2
      end if
c
c     Generate new time dependent part if necessary
c     to be able to recreate the trip of restarted simulations,
c     loop from ntdt=-1 up to present trip count.
c
      do i=ntdt+1,int(time/tdt)
         do k=1,nwalls
            do z=1,kpts(k)
               fzt3(z,k)=fzt4(z,k)
            end do
         enddo
c
c     Get random distribution and rescale
c
         do k=1,nwalls
c            if (k.le.4) then
c               call rand_func(fzt4(1,k),znek(1,k),kpts(k),seed1,
c     &              num_modes)
c            elseif (k.gt.4.and.k.le.8) then
c               call rand_func(fzt4(1,k),znek(1,k),kpts(k),seed2,
c     &              num_modes)
c            else
c               call rand_func(fzt4(1,k),znek(1,k),kpts(k),seed3,
c     &              num_modes)
c            endif
            call rand_func(fzt4(1,k),znek(1,k),kpts(k),seed,num_modes)
            do z=1,kpts(k)
               fzt4(z,k)=tampt*fzt4(z,k)
            enddo
         enddo
      enddo
c
c     Update trip count as actual time divided by time scale
c
      ntdt=int(time/tdt)
c
c     Generate the z-dependence of the trip
c     as a smooth transition between old and new trip vectors
c     p is varying from 0 to 1 for a given trip count.
c
      p=(time-real(ntdt)*tdt)/tdt
      b=p*p*(3.-2.*p)
      do k=1,nwalls
         do z=1,kpts(k)
            fzt1(z,k)=fzt2(z,k)+(1.-b)*fzt3(z,k)+b*fzt4(z,k)
         enddo
      enddo

      end subroutine tripf
c-----------------------------------------------------------------------
      subroutine rand_func(rand_vec,zvec,zpts,seed,num_modes)

      implicit none

      integer seed,k
      integer zpts
      real zvec(1:zpts),bb
      real rand_vec(zpts)

      real pi
      parameter (pi = 3.1415926535897932385)
      integer num_modes
c      parameter (num_modes = 10)
c
c     Local variables
c
      integer z,m
      real zlength
      real phase
      real theta
c
c     External function
c
      real ran2
c
c     Compute length of z-interval
c
      zlength = zvec(zpts) - zvec(1)
      if (zlength .eq. 0.) zlength = 1.
      do z=1,zpts
         rand_vec(z) = 0.0
      enddo
c
c     Compute m sinus modes
c
      do m=1,num_modes
         bb = ran2(seed)
         phase = 2.*pi*bb
         do z=1,zpts
            theta = 2.*pi*m*zvec(z)/zlength
            rand_vec(z) = rand_vec(z) + sin(theta + phase)
         enddo
      enddo

      end subroutine rand_func



      real function ran2(idum)
c
c     A simple portable random number generator
c
c     Requires 32-bit integer arithmetic
c     Taken from Numerical Recipes, William Press et al.
c     gives correlation free random numbers but does not have a very large
c     dynamic range, i.e only generates 714025 different numbers
c     for other use consult the above
c     Set idum negative for initialization
c
      implicit none

      integer idum,ir(97),m,ia,ic,iff,iy,j
      real rm
      parameter (m=714025,ia=1366,ic=150889,rm=1./m)
      save iff,ir,iy
      data iff /0/

      if (idum.lt.0.or.iff.eq.0) then
c
c     Initialize
c
         iff=1
         idum=mod(ic-idum,m)
         do j=1,97
            idum=mod(ia*idum+ic,m)
            ir(j)=idum
         end do
         idum=mod(ia*idum+ic,m)
         iy=idum
      end if
c
c     Generate random number
c
      j=1+(97*iy)/m
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum

      end function ran2
c----------------------------------------------------------------------
      
      subroutine usrdat()
c     Read user-defined boundary conditions into a common array bocoarray
c     Read input data for fringe region

      include 'SIZE'

c     variables used for user-defined boundary conditions

      real dum(3)

!      call user_param
      call uprm_read            ! New user parameter read function
      
      return
      end
      
c-----------------------------------------------------------------------

      subroutine usrdat2()  
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
     
      return
      end

c-----------------------------------------------------------------------
      
      subroutine usrdat3()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      return
      end
c---------------------------------------------------------------------- 
c
c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end

C=======================================================================

c
c automatically added by makenek
      subroutine usrflt(rmult) ! user defined filter
      include 'SIZE'
      real rmult(lx1)
      call rone(rmult,lx1)
      return
      end
c
c automatically added by makenek
      subroutine userflux ! user defined flux
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real fluxout(lx1*lz1)

      return
      end
c
c automatically added by makenek
      subroutine userEOS ! user defined EOS 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      return
      end
