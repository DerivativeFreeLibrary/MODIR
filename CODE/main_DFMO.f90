subroutine main_DFMO(n,m_vin,q,bl,bu,num_funct,first_time_call,maxnf_DFMO,istop)

	use eps_mod
	use filtro
	use popolazione
	use cache
	use vincoli
	use alfa_mod
	implicit none

	integer :: n,p,q,m_vin,maxnf_DFMO
	integer			 :: n_int
	real*8	:: bl(n), bu(n)
	logical :: first_time_call

	integer	:: icontr,numvar,j,istop,nobjn,icheck
	real		:: tbegin, tend
	character*30 :: nomefun
	integer ::i,k

	real*8, allocatable :: x(:),step(:)    
	real*8, allocatable :: punti(:,:)
	real*8, allocatable :: ceq(:), ciq(:)
	real*8, allocatable :: fob(:), finiz(:), alfaciniz(:),f(:)
	real*8, allocatable :: fob_media(:), ciq_media(:)

	integer ::            num_funct
	real*8             :: delta 
	real*8			   :: fapp, alfainiz
	real*8			   :: violiniz, violint, fex
	real*8             :: alfa_stop
	integer            :: nf_max,iprint
	logical			:: flag
    integer t

	  p = 0
	  write(*,*) 'n=',n,' q=',q,' m=',m,' mm=',mm,' p=',p
	  !pause
	  nelem = 0

      t=npop
      
	  allocate(alfaciniz(n))
	  allocate(x(n),step(n))
	  allocate(punti(n,t))
        
	  allocate(fob(q),finiz(q), f(q))
	  allocate(fob_media(q))

	if(m.ge.1) then
		allocate(ciq(m),ciq_media(m))
	endif

	  write(*,*) 'after allocate''s',n,q,m

	  !call check_feasibility(n,m,q)

	  fob = 1.d0
	  finiz = 0.d0

 	  call cpu_time(tbegin)

	  alfainiz=1.0d1

	  alfainiz=min(10.d0,maxval(bu-bl)/10.d0)
      !alfainiz=min(10.d0,maxval(bu-bl)/2.d0)

	  do i = 1,n
        	alfaciniz(i) = min(10.d0,(bu(i)-bl(i))/10.d0)
	  	!alfaciniz(i) = min(10.d0,(bu(i)-bl(i))/2.d0)
	  enddo

	  !==========================================
      !    SOSTITUIRE CON INIZIALIZZAZIONE LISTA
	  !==========================================

!----------- calcolo esp con medie ---------------------------

		if(.false.) then
		
				write(*,*) 'call diagonale...'
				call diagonale(n,t,bl,bu,punti)

				fob_media = 0.d0
				if (m.ge.1) then
				    ciq_media = 0.d0
				endif

				if (n>1) then 

				    do i = 1,n
				        x = punti(:,i)
				        call functs(n,x,q,fob)
				        if(m.ge.1) then
				            call fconstriq(n,m,x,ciq)
				            if(maxval(ciq) >= 0.d0) then
					            fob_media = fob_media + fob
					            ciq_media = ciq_media + ciq
				            endif
				        endif
				    enddo
				    fob_media = fob_media / dble(n)
				    if(m.ge.1) then
				        ciq_media = ciq_media / dble(n)
				    endif
				else

				    x=punti(1,1)

				    call functs(n,x,q,fob)
				    if(m.ge.1) then
				        call fconstriq(n,m,x,ciq)
				        if(maxval(ciq) <= 0.d0) then
				            fob_media =  fob
				            ciq_media =  ciq
				        endif
				    endif

				endif


			!	write(*,*) '--------------------------------------'
			!	write(*,*) '--------- EPSILON INIZIALI -----------'
			!	write(*,*) '--------------------------------------'
				do k=1,q
				    write(*,112) k		
				    do j = 1,m
				        epsiq(k,j)=1.0d-2*dmax1( 1.d-3, dmax1( 0.0d0, ciq_media(j) ) ) / dmax1( 1.0d-3, dabs( fob_media(k) ) )
				        write(*,111) epsiq(k,j), ciq_media(j), fob_media(k)
				    end do
				    write(*,*)
				    !!!pause
				end do
			112 format(' f.ob. ',i3,': ',$)
			111 format(es16.9,1x,'(ciq_medio=',es13.6,')','( f_medio=',es13.6,')',1x,$) 
		endif
! ------   fina cocolo eps con medie  ------------------------------------------------------------                
      
           
    if(.false.) then      !false per eliminare la popolazione iniziale
        
       write(*,*) 'call latin hypercube for Lpop init...'
	   call latin(t,n,bl,bu,punti)
	   do i = 1,t
			x = punti(:,i)
			write(*,113) x

			if(.not.find_in_cache(n,q,m,x,fob,ciq)) then   
				call functs(n,x,q,fob)
				if(m.ge.1) call fconstriq(n,m,x,ciq)
				call insert_in_cache(n,q,m,x,fob,ciq)
			endif
			call functs_pen(n, q, m, x, finiz, fob, ciq)	
			call insert_point(n,q,x,finiz,alfainiz, alfaciniz)
			write(*,*) 'after ins: finiz', finiz
       enddo
       
    else  
        
	    nelem = 0
        
    endif
    
    
113	format(10(1x,es16.9))

	if(first_time_call) then
    
		call startp(n,x)
		write(*,113) x

		write(*,*) 'call diagonale...'
		call diagonale(n,t,bl,bu,punti)

		if(.not.find_in_cache(n,q,m,x,fob,ciq)) then   
		call functs(n,x,q,fob)
			if(m.ge.1) call fconstriq(n,m,x,ciq)
			call insert_in_cache(n,q,m,x,fob,ciq)
		endif
		call functs_pen(n, q, m, x, finiz, fob, ciq)	
		call update_point(n,q,x,finiz,alfainiz, alfaciniz, flag)
	 	!write(*,*) flag,fob,finiz !,max(0.d0,ciq)

		if(.true.) then
			t = n
			do i = 1,t
				x = punti(:,i)
				write(*,113) x

				if(.not.find_in_cache(n,q,m,x,fob,ciq)) then   
					call functs(n,x,q,fob)
					if(m.ge.1) call fconstriq(n,m,x,ciq)
					call insert_in_cache(n,q,m,x,fob,ciq)
				endif
				call functs_pen(n, q, m, x, finiz, fob, ciq)	
				if(intorno(n,q,x,finiz,1.d-2)) then
					call update_point(n,q,x,finiz,alfainiz, alfaciniz, flag)
					!write(*,*) flag,fob,finiz !,max(0.d0,ciq)
				endif
			enddo
		endif
	endif

	!call print_filter(n,q)
	!pause

        write(*,*) ' ------------------------------------------------- '

        write(*,*) '      f(xo) = ',fob

        viol=0.d0

        do i = 1,m
	   	   viol=viol+max(0.d0,ciq(i))
	    enddo
        do i = 1,p
	   	   viol=viol+abs(ceq(i))
	    enddo

        write(*,*) '  cviol(xo) = ',viol

        write(*,*) ' ------------------------------------------------- '
	    do i=1,n
		   write(*,*) ' xo(',i,') =',x(i)
        enddo

        write(*,*) ' ------------------------------------------------- '

!-----------------------------------------------------------------------
!	choice of starting penalty parameter values
!-----------------------------------------------------------------------


	   num_funct   = 0 
	  

	   if(m >= 1) epsiq_in     = epsiq
!	   if(p >= 1) epseq_in     = epseq

	   finiz       = fob
	   violiniz    = viol

1	   alfa_stop=1.d-9
  	   nf_max=maxnf_DFMO !-500*n !20000
	   iprint=0


		!====================================================================================================

	   !write(*,*)
	   !write(*,*) ' call the optimizer ...',n,q,m
	   !pause

		istop = 1
        call sd_box(n,q,bl,bu,alfa_stop,nf_max,num_funct,iprint,istop,2)
	!write(3,*) n

	   !write(*,*) ' ... done'
	   !write(*,*)


	   call cpu_time(tend)

      
	   write(*,*) ' ------------------------------------------------------------------------------'     
	   if(istop.eq.1) then
         write(*,*)  ' STOP - stopping criterion satisfied. alfa_max <= ',alfa_stop
	   endif
       if(istop.eq.2) then
         write(*,*)  ' STOP - max number function evaluations exceeded. nf > ',nf_max
	   endif

       if(istop.eq.3) then
         write(*,*)  ' STOP - max number iterations exceeded. ni > ',nf_max
	   endif

	   write(*,*) ' total time:',tend-tbegin
       write(*,*) ' number of function evaluations = ',num_funct 
	   write(*,*) ' ------------------------------------------------------------------------------'  

		write(*,*)

        write(*,*) ' ------------------------------------------------- '
   
	   write(*,*) 

133 format(a30,' & ',i6,' & ',es15.6,'\\\hline')

	  deallocate(alfaciniz)
	  deallocate(x,step)
	  deallocate(punti)

	  deallocate(fob,finiz, f, fob_media)

	if(m.ge.1) then
		deallocate (ciq,ciq_media)
	endif

end subroutine main_DFMO

subroutine copy_filter(n,q,noutdir,nindfmo)
	use mod_type
	use alfa_mod
	use cache
	use filtro
	use paretoset, only : Ldir,maxdim

	implicit none
	integer			:: n, q, i, noutdir, nindfmo
	real*8			:: deltaf(ndim)
	real*8			:: ysud(n), ytemp(n), xsud(n)
	real*8			:: xbar(n), lbs(n), ubs(n)
	logical 		:: flag
	
	write(*,*) 'copy_filter: INGRESSO'

	xbar = 0.5d0
	lbs  = 0.d0
	ubs  = 1.d0

	noutdir = 0

	do i=1,maxdim

		if(Ldir(i)%flag_busy) then
			noutdir = noutdir + 1
			!write(*,*) 'transfer_filter: PRIMA di temp%cent'
			ysud = Ldir(i)%x
			xsud = ysud
			!call unscalevars(n,ysud,ytemp,xbar,lbs,ubs)
			!call unscalevars_direct(n,ytemp,xsud)
			!write(*,*) 'transfer_filter:  xsud = ',xsud
			!write(*,*) 'transfer_filter:  n,q,mm = ',n,q,mm

			if(.not.find_in_cache(n,q,mm,xsud,fob_m,ciq_m)) then   
				!write(*,*) 'transfer_filter:  NOT IN CACHE'
				call functs(n,xsud,q,fob_m)
				!write(*,*) 'transfer_filter:  AFTER FUNCTS'
				if(mm.ge.1) call fconstriq(n,mm,xsud,ciq_m)
				!write(*,*) 'transfer_filter:  AFTER FCONSTRIQ'
				call insert_in_cache(n,q,mm,xsud,fob_m,ciq_m)
			endif
			!write(*,*) 'transfer_filter:  CALL FUNCTS_PEN'
			call functs_pen(n, q, mm, xsud, finiz_m, fob_m, ciq_m)	
			!if(intorno(n,q,xsud,finiz_m,1.d-9)) then
			!	write(*,*) 'transfer_filter:  INTORNO'
				call update_point(n,q,xsud,finiz_m,alfainiz_m, alfaciniz_m, flag)
			!	write(*,*) flag, xsud,finiz_m
			!endif
			!---------------------------------------------------
		endif
	enddo

!	call spread_ord(n,q,deltaf,flag)
!	
!	nindfmo = 0
!
!	do i = 1,ndim
!		if(i <= 10*n) then
!			if(Lkappa(i,n+q+1) > 50.d0) then
!				nindfmo = nindfmo + 1
!			endif
!		else
!			Lkappa(i,n+q+1) = 0.d0
!		endif
!	enddo

	write(*,*) 'copy_filter: USCITA'

	return
end subroutine copy_filter

subroutine funct1(n,x,f)
	implicit none
	integer n,indfun
	real*8  x(n),f

	return
end subroutine funct1

subroutine functs_pen(n, q, m, x, fpen, fob, ciq)
	use eps_mod
	implicit none
	integer n,i,j,m,q
	real*8  fpen(q)
	real*8  x(n),f(q),fob(q),fmax(q),ciq(m), viol

	fmax = 0.d0
	viol = 0.d0

	do i = 1,m
		do j=1,q
                  fmax(j) = fmax(j) + max(0.d0,ciq(i))/epsiq(j,i)
		end do
		viol=viol+max(0.d0,ciq(i))
	enddo
        
	fpen = fob + fmax

	return
end subroutine functs_pen

subroutine fconstreq(n,p,x,ceq)
	implicit none
	integer n,p
	real*8 x(n), ceq(p)
	return
end subroutine fconstreq

subroutine diagonale(n,t,l,u,punti)
	implicit none
	integer		:: n,t
	real*8			:: l(n), u(n)
	real*8			:: punti(n,t)
	integer		:: i

	if(n > 1) then
		do i = 1,n
			punti(:,i) = l + dble(i-1)/dble(n-1)*(u-l)
		enddo
	else
		punti(1,1) = l(1)
	endif

	return	
end subroutine diagonale

subroutine latin(t,n,l,u,punti)
	implicit none
        integer         :: j
	integer		:: t,n
	real*8			:: l(n), u(n)
	real*8			:: punti(n,t)
	integer		:: i, seed
        !call get_seed ( seed )
	seed = 1
        call random_initialize ( seed )
	call latin_random ( n, t, seed, punti )
        
        do i=1,n
            do j=1,t
	       punti(i,j)= l(i)+0.01+(u(i)-0.01-l(i)-0.01) *punti(i,j)
            end do
        end do

	return	
end subroutine latin
