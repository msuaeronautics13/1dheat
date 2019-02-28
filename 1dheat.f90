      program 1dheat      

          implicit none

! ******* Misc          
          integer :: j,jl
          integer :: ier,tmp
          integer :: itime
          integer,dimension(1) :: start,length
          
! ******* Grid
          real*8, dimension(:), pointer :: x
          integer, dimension(:), pointer :: el,narr
          integer :: nPts,nEl
          character(len=1024) :: gridfile

! ******* 1D heat conduction          
          real*8, dimension(:), pointer :: T1,Tb,Tg
          real*8, dimension(:), pointer :: a,b,c,d,T,flux,cfdflux
          real*8 :: T2,k,L,alpha,rho,cp,Fo,dt,q,m
          integer :: n,kk

! ******* VTK          
          character(len=1024) :: vtkfile

          real*8 :: timestep
      
          call getarg(1,gridfile)
            open(100,file=gridfile)
            read(100,*) nPts
            allocate(x(3*nPts),narr(nPts))
            do jl=1,nPts
              read(100,*) narr(jl), x(3*(jl-1)+1), x(3*(jl-1)+2),x(3*(jl-1)+3)
              narr(jl)=jl
            enddo
           ! do jl=1,nEl
           !  read(100,*) tmp,el(jl),el(jl+nEl),el(jl+2*nEl),el(jl+3*nEl)
!              write(*,*) tmp,el(jl),el(jl+nEl),el(jl+2*nEl),el(jl+3*nEl)
           ! enddo
          close(100)         
          allocate(Tb(nPts),Tg(npts),T1(nPts),cfdflux(nPts))

!********* Solve 1D heat relation (CN) ************
!          FTCS ~ O(dt**1) 
!          AB ~ O(dt**2)
!          Through-thickness model
!        * Tg - ghost node
!        |
!--------* Tb - boundary node-------
!        | 
!        * T1 - interior node
!        | 
!--------* backside node (300K)-----

          ! Properties
          k=215
          L=0.01
          T2=300 ! BC
          rho=2700
          cp=900
          alpha=k/(rho*cp)
          dt=0.1d0
          Fo=alpha*dt/L**2

          ! Initial conditions
          if (itime .eq. 1) then 
            T1(:) = 300.0
            Tb(:) = 300.0
          endif

          do j=1,nPts
           ! 300K backside
            Tg(j) = 2*Tb(j)-T1(j) ! Ghost node
          enddo
           
          allocate(T(4*nPts),a(4*nPts),b(4*nPts),c(4*nPts),d(4*nPts),flux(4*nPts))

!         Construct T           
          n=2
          flux(:) = 0.0d0
          do j=1,nPts
            T(n)=Tb(j)
            T(n-1)=Tg(j)
            T(n+1)=T1(j)
            !T(n+2)=300.0d0 ! T2 (300K)            
            T(n+2)=T(n+1) ! T2 dT/dX=0
            flux(n)=cfdflux(j)! Q_load
            flux(n+1)=(T(n)-T(n+1))/L ! Q_cond  
            flux(n+2)=(T(n+1)-T(n+2))/L  
            n=n+4
          enddo

          nn=4*nPts ! Total number of nodes, real and through thickness

          ! Assemble coefficients
          do j=1,nn
            a(j)=-alpha/(2*L**2)
            b(j)=(1/dt) + (alpha/L**2)
            c(j)=a(j)
            d(j)=a(j)*T(j-1)+(1/dt+a(j)+c(j))*T(j)+c(j)*T(j+1)+flux(j) 
          enddo

          do j=2,nn
            m = a(j)/b(j-1)
            b(j) = b(j) - m*c(j-1)
            d(j) = d(j) - m*d(j-1)
          enddo

          ! Back substitution
          T(nn) = d(nn)/b(nn)
          kk=nn
          !do j=2,nn
          do kk=nn-1,1,-1
           ! kk=nn+1-j
            T(kk) = (d(kk) - c(kk)*T(kk+1))/b(kk)
          enddo
            
          n=0
          do j=1,nPts            
            T1(j) = T(4*n+3)
            Tg(j) = T(4*n+1)
            Tb(j) = T(4*n+2)
            n=n+1
          enddo
        enddo ! Time loop

      end program










          
