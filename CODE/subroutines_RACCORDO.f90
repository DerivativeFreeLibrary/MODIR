!============================================================================================
!    MODIR - DIRECT Algorithm for derivative-free multiobjective optimization.  
!    A Derivative-Free algorithm for bound  constrained multiobjective global 
!	 optimization problems proposed by 
!	 E.F.Campana, M.Diez, G.Liuzzi, S.Lucidi, R.Pellegrini, V.Piccialli, F.Rinaldi, A.Serani 
!	 (see Ref. below)
!
!    Copyright (C) 2017  E.F.Campana, M.Diez, G.Liuzzi, S.Lucidi, 
!						 R.Pellegrini, V.Piccialli, F.Rinaldi, A.Serani
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    E.F.Campana, M.Diez, G.Liuzzi, S.Lucidi, R.Pellegrini, V.Piccialli, F.Rinaldi, A.Serani. 
!	 A Multi-objective DIRECT algorithm for ship hull optimization, 
!	 submitted for pubblication on Computational Optimization and Applications (2017).
!
!==============================================================================================
subroutine wrapper(maxnf_MODIR,maxnf_DFMO)
	use mod_type
	use alfa_mod
	use eps_mod
	use vincoli
	use filtro
	use paretoset
	use paretodim, only: init_paretodim
	use popolazione
	use cache
	use mod_box
	implicit none

	interface
		subroutine main_MODIRECT(n,nobj, lb0, ub0, nftot, root,first_time_call,maxnf_MODIRECT)
			use mod_type
			implicit none

			type(colonna), pointer	:: root
			integer			:: n, nobj, nftot, maxnf_MODIRECT
			real*8			:: lb0(n), ub0(n)
			logical			:: first_time_call
		end subroutine main_MODIRECT
		subroutine dealloca_struct(root)
			use mod_type
			implicit none

			type(colonna), pointer	:: root
		end subroutine dealloca_struct
		subroutine crea_nondom_DIRECT(n,q,root)
			use mod_type
			implicit none

			integer						:: n, q
			type(colonna), pointer		:: root
		end subroutine crea_nondom_DIRECT
		subroutine reset_nondom_DIRECT(root)
			use mod_type
			implicit none

			type(colonna), pointer		:: root
		end subroutine reset_nondom_DIRECT
	end interface

	type(colonna), pointer		:: root
	integer						:: which_first, which_tocall
	!							   which_first = 1 --> parte prima il locale (DFMO)
	!							   which_first = 2 --> parte prima il globale (MODIRECT)
	logical						:: first_time_call, stop_DFMO
	character*40				:: display_line
	integer						:: n, mv, qq
	integer						:: i, k, nfparz, nftot, j, pn, istop
	integer						:: noutdir, nindfmo
	integer						:: maxnf, maxnf_DFMO, maxnf_MODIR
	integer						:: ncall_MODIR, ncall_DFMO
	real*8, allocatable			:: bl(:), bu(:), x(:), ciq(:), fob(:)
	integer						:: icicli
	integer, parameter			:: ncicli = 2
	integer						:: maxfevals(ncicli)

	!-----------------------------------------------
	! CHECK DI CONSISTENZA TRA FILTRO(DFMO) E PARETOSET(MODIRECT)
	!   ndim: dimensione del filtro
	! maxdim: dimensione di paretoset
	!
	! DEVE RISULTARE  ndim = maxdim 
	!-----------------------------------------------
	if( ndim/=maxdim ) then
		write(*,*) 
		write(*,*) 'ATTENZIONE: ndim /= maxdim'
		write(*,*) '   ndim = ',ndim
		write(*,*) ' maxdim = ',maxdim
		write(*,*) 
		stop
	endif

	!-----------------------------------------------
	! DIMENSIONAMENTO DEL PROBLEMA
	!-----------------------------------------------
	write(*,*) 'set dimensions ...'
	call setdim(n,mv,qq)
	m  = mv
	mm = mv
	q  = qq

	write(*,*) 'main: n=',n,' m=',m,' mm=',mm

	!-----------------------------------------------
	! ALLOCAZIONE
	!-----------------------------------------------
	allocate(lb(n),ub(n),xtemp(n),ytemp(n),xbar(n),lbs(n),ubs(n))
	allocate(bl(n),bu(n),x(n))
	allocate(fob(q))
	allocate(fob_m(q),finiz_m(q))
	allocate(alfaciniz_m(n))

	if(m.ge.1) then
		allocate (epsiq(q,m),epsiq_in(q,m))
		allocate(eps(q,m),constr(m),epsiniz(q,m))
		allocate(ciq(m),ciq_m(mm))
	endif

	!-----------------------------------------------
	! SET DEI BOUNDS
	!-----------------------------------------------
	write(*,*) 'set bounds ...'
	call setbounds(n,bl,bu)

	lb       = bl
	ub       = bu
	lbs      = 0.d0
	ubs      = 1.d0
	xbar     = (ubs+lbs)/2.d0
	fattore  = 1.0d0
	ampiezza = 1.0d0

	write(*,*) 'bl = ',bl
	write(*,*) 'bu = ',bu

	!-----------------------------------------------
	! INIZIALIZZAZIONI
	!-----------------------------------------------
	call setup_cache(n,q,m)
	call init_paretoset(n,q)
	call init_paretodim(n,q)

	nullify(root)
	npop = n 
	allocate(Lkappa(ndim,n+q+2+n),Lnew(ndim,n+q+2+n),Ltilde(ndim,n+q+2+n))
	allocate(Lpop(npop,n+q+2+n))
	Lkappa(:,n+q+1) = 0.d0
	Lnew  (:,n+q+1) = 0.d0
	Ltilde(:,n+q+1) = 0.d0
	Lpop  (:,n+q+1) = 0.d0
	nftot           = 0

	which_first = 2

	maxfevals(1) = 500*n
	maxfevals(2) = max(0,maxnf - maxfevals(1))

	alfainiz_m=min(10.d0,maxval(bu-bl)/10.d0)

	do i = 1,n
		alfaciniz_m(i) = min(10.d0,(bu(i)-bl(i))/10.d0)
	enddo

	!-----------------------------------------------
	! CALCOLO EPSILON INIZIALI (per vincolato)
	!-----------------------------------------------
	if(.true..and.(m >= 1)) then   ! .false.   eps fisso

		write(*,*) 'call startp...'
		call startp(n,x)

		do i=1,n
			if((x(i).lt.bl(i)).or.(x(i).gt.bu(i))) then
				write(*,*) 'The starting point violates bound constraints, STOP'
				write(*,*) bl
				write(*,*) x
				write(*,*) bu
				stop
			endif
		enddo

		call functs(n,x,q,fob)
		if(m.ge.1) call fconstriq(n,m,x,ciq)
		nftot = nftot+1

        !-----------------------------------------------------------------------
        !     print starting point info
        !-----------------------------------------------------------------------
		write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
		write(*,*) 'fob = ',fob
		write(*,*) 'ciq = ',ciq
		write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'


        2002    format(2d20.10)

		write(*,*) '--------------------------------------'
		write(*,*) '--------- EPSILON INIZIALI -----------'
		write(*,*) '--------------------------------------'
		do k = 1,q
			do i = 1,m
				if(max(0.d0,ciq(i)) < 1.d-0) then
					epsiq(k,i) = 1.d-3
				else
					epsiq(k,i) = 1.d-1
					!epsiq(k,i) = 1.d-3 !eps fissato
				endif
			enddo
		enddo
        
	else
	  
		if(m>=1) epsiq=1.d-0
		
	endif !fine if eps fisso
       

	select case (which_first)
	case (1)
		display_line = "Parte per primo il locale  (DFMO)"
	case (2)
		display_line = "Parte per primo il globale (MODIRECT)"
	case default
		display_line = "Parte per primo il globale (MODIRECT)"
	end select

	!write(*,*) 
	!write(*,*) display_line
	!write(*,*)
	!pause

	which_tocall    = which_first
	first_time_call = .true.
	stop_DFMO       = .false.
	ncall_MODIR     = 0
	ncall_DFMO      = 0

	if(maxnf_MODIR > 0) then
		call main_MODIRECT(n,q,bl,bu,nfparz,root,first_time_call,maxnf_MODIR)
		ncall_MODIR = ncall_MODIR + 1
		write(*,*) 'MODIRECT completed. Passing nondom to DFMO...',n,q
		!pause
		!call init_nondom_DFMO(n,q,bl,bu)
		call copy_filter(n,q,noutdir,nindfmo)
		nftot = nftot + nfparz
	endif

	if(maxnf_DFMO > 0) then
		call main_DFMO(n,mv,q,bl,bu,nfparz,first_time_call,maxnf_DFMO,istop)
		ncall_DFMO = ncall_DFMO + 1
		if(nfparz < maxnf_DFMO) then
			stop_DFMO = .true.
		endif
		write(*,*) 'DFMO completed.'
		call reset_nondom_DIRECT(root)
		call crea_nondom_DIRECT(n,q,root)
		first_time_call = .false.
		nftot = nftot + nfparz
	endif

	!-----------------------------------------------
	! STAMPA DEI RISULTATI
	!-----------------------------------------------
	!write(2,*) '     LISTA  '
	!write(*,*) '     LISTA  '

	pn = 0
	j = 0

	write(*,*) 'Nondom solutions printed to file: out_nondom.txt'
	open(2,file='out_nondom.txt',status='replace')

	do i=1,ndim
		if(Lkappa(i,n+q+1)>50.0d0) then
			j = j+1
			if(m.ge.1) then
				call fconstriq(n,m,Lkappa(i,1:n),ciq)
				viol = max(0.d0,maxval(ciq))
			else
				viol = 0.d0
			endif

			!write(*,*) i, Lkappa(i,1:n), viol
			if(viol <= 1.0d-3) then
				pn = pn+1
				!call functs(n,Lkappa(i,1:n),q,Lkappa(i,n+1:n+q))
				write(2,*) Lkappa(i,n+1:n+q)
				!write(*,*) i, Lkappa(i,n+1:n+q)
			endif
		end if
	end do
	close(2)
	!write(3,700) n, q, noutdir,nindfmo,j, pn, nftot
	write(*,*) 'Nondom. = ',j,' di cui amm.',pn

700   format(' & ',I5,' & ',I5,' & ',I7,' & ',I7,' & ',I7,' & ',I7,' & ',I7,' \\\hline')

	!-----------------------------------------------
	! DEALLOCAZIONE
	!-----------------------------------------------
	call destroy_cache()
	call canc_paretoset()
	call dealloca_struct(root)

	deallocate(lb,ub,xtemp,ytemp,xbar,lbs,ubs)
	deallocate(bl,bu,x)
	deallocate(fob)
	deallocate(fob_m,finiz_m)
	deallocate(alfaciniz_m)
	if(m.ge.1) then
		deallocate (epsiq,epsiq_in)
		deallocate(eps,constr,epsiniz)
		deallocate(ciq,ciq_m)
	endif
	deallocate(Lkappa,Lnew,Ltilde)
	deallocate(Lpop)

end subroutine wrapper

!-----------------------------------------------------------------------
! La subroutine reset_nondom_DIRECT dealloca la lista dei non-dominati
! associati a ciascun intervallo della partizione corrente
!
! Idealmente, la sequenza di chiamata delle subroutine
! reset_nondom_DIRECT e crea_nondom_DIRECT dovrebbe essere:
!
!   call reset_nondom_DIRECT(root)
!   call crea_nondom_DIRECT(n,q,root)
!
! cioè, PRIMA si azzerano le liste di nondominati e
! POI si creano quelle nuove.
!-----------------------------------------------------------------------
subroutine reset_nondom_DIRECT(root)
	use mod_type
	implicit none

	type(colonna), pointer		:: root, temp
	type(intervallo), pointer	:: curr
	type(typ_lista), pointer	:: pnondom
	
	temp => root

	do while (associated(temp))
		curr => temp%int
		do while (associated(curr))

			curr%nnondom = 0
			pnondom => curr%list_nondom
			do while (associated(pnondom))
				pnondom => pnondom%next
				nullify(curr%list_nondom%next)
				deallocate(curr%list_nondom%fobs,curr%list_nondom%x)
				deallocate(curr%list_nondom)
				curr%list_nondom => pnondom
			enddo

			curr => curr%next	
		enddo
		temp => temp%next
	enddo

end subroutine reset_nondom_DIRECT

!-----------------------------------------------------------------------
! La subroutine crea_nondom_DIRECT considera in ingresso una lista
! di non-dominati Lkappa (generata da DFMO) e con questa aggiorna la
! struttura dati di DIRECT aggiungendo a ciascun intervallo i punti
! non-dominati di Lkappa che sono in lui contenuti
!-----------------------------------------------------------------------
subroutine crea_nondom_DIRECT(n,q,root)
	use filtro, only : Lkappa, ndim
	use paretoset, only: Ldir, update_point, print_filter
	use paretodim, only: Ldim, update_point_dim
	use mod_type
	use mod_box
	implicit none

	type(colonna), pointer		:: root, temp
	type(intervallo), pointer	:: curr
	type(typ_lista), pointer	:: templist
	integer				:: n, q
	integer				:: ip, j
	!real*8				:: xtemp(n), lb(n), ub(n)
	real*8				:: lb1(n), ub1(n), fdim(q+1)
	real*8				:: norma
	logical				:: soddisfa, flag, trovato

	if(.not.associated(root)) then
	! la struttura dati di MODIRECT non è allocata, allocala
		write(*,*) 'root pointer non associato'
		allocate(root)
		nullify(root%next)
		nullify(root%pred)

		allocate(root%int)
		call alloca_intervallo(root%int,n,q)
		root%int%cent    = 0.5d0
		root%int%dimen   = 1.d0
		root%int%maxdim  = 1.d0
		root%int%der     = 0.d0
		root%int%xbars   = 0.5d0
		root%int%lbs     = 0.d0
		root%int%ubs     = 1.d0

		root%int%diam    = norma(n,root%int%dimen)/2.d0
		root%diam        = norma(n,root%int%dimen)/2.d0
		root%int%flagloc = .false.
		root%int%flagdiv = .true.
		root%int%flagcon = .false.
		root%int%flagopt = .false.
		root%int%id      = 1
		root%int%nnondom = 0
		nullify(root%int%list_nondom)
		nullify(root%int%next)
		nullify(root%int%pred)

		call unscalevars(n,root%int%cent,ytemp,xbar,lbs,ubs)
		call unscalevars_direct(n,ytemp,xtemp)
	    	ytemp = xtemp

		write(*,*) 'calling funct...'
		call funct(n,xtemp,root%int%fint,root%int%fobs)

		write(*,*) 'calling update_point...'
		call update_point(n,q,xtemp,root%int%fobs,root%int,flag)
		!nf          = 1
		!nint        = 1
	endif

	do ip = 1,ndim

		if(Lkappa(ip,n+q+1) > 50.d0) then
			!write(*,*) 'insert del punto pn=',ip
			!cerca nella struttura di DIRECT l'intervallo di riferimento

			temp => root
			trovato = .false.

			do while (associated(temp))
				curr => temp%int
				do while (associated(curr))

					call unscalevars(n,curr%cent - curr%dimen/2.d0,xtemp, curr%xbars,curr%lbs,curr%ubs)
					call unscalevars_direct(n,xtemp,lb1)

					call unscalevars(n,curr%cent + curr%dimen/2.d0,xtemp,curr%xbars,curr%lbs,curr%ubs)
					call unscalevars_direct(n,xtemp,ub1)
					
					!write(*,*) 'id = ',curr%id
					!write(*,*) 'ls = ',curr%lbs
					!write(*,*) 'us = ',curr%ubs
					!write(*,*) 'xs = ',curr%xbars
					!write(*,*) 'lb = ',lb,curr%cent - curr%dimen/2.d0
					!write(*,*) 'ub = ',ub,curr%cent + curr%dimen/2.d0
					!write(*,*) ' x = ',Lkappa(ip,1:n)

					!se il punto Lkappa(ip,1:n) soddisfa lb e ub aggiungi
					soddisfa = .true.
					do j = 1,n
						if( (Lkappa(ip,j) > ub1(j)).or.(Lkappa(ip,j) < lb1(j)) ) then
							soddisfa = .false.
							exit
						endif
					enddo
					if( soddisfa ) then
						!l'intervallo corrente contiene il punto Lkappa(ip,1:n)
						!quindi il punto va aggiunto alla lista dei non-dominati dell'intervallo
						allocate(templist)
						allocate(templist%fobs(q))
						allocate(templist%x(n))
						templist%x    = Lkappa(ip,1:n)
						templist%fobs = Lkappa(ip,n+1:n+q)
						nullify(templist%prev)
						nullify(templist%next)
						if(associated(curr%list_nondom)) then
							templist%next => curr%list_nondom
							curr%list_nondom%prev => templist
							curr%list_nondom => templist
						else
							curr%list_nondom => templist
						endif
						nullify(templist)
						curr%nnondom          = curr%nnondom+1
						trovato               = .true.
						call update_point(n,q,Lkappa(ip,1:n),Lkappa(ip,n+1:n+q),curr,flag)
						fdim(1:q) = Lkappa(ip,n+1:n+q)
						fdim(q+1) = -curr%diam
						call update_point_dim(n,q,Lkappa(ip,1:n),fdim,curr,flag)
					endif
					if(trovato) exit
					curr => curr%next	
				enddo
				if(trovato) exit
				temp => temp%next
			enddo

		endif

	enddo
	!write(*,*) 'crea_nondom_DIRECT: call print_filter(n,q)'
	!call print_filter(n,q)
	!write(*,*) 'crea_nondom_DIRECT: done call print_filter(n,q)'

end subroutine crea_nondom_DIRECT

subroutine init_nondom_DFMO(n,q,bl,bu)
	use filtro
	use paretoset, only: maxdim, Ldir
	implicit none

	integer					:: n, q
	integer					:: i, j, k
	real*8					:: alfainiz
	real*8					:: bl(n), bu(n), alfaciniz(n) 

	write(*,*) 'init_nondom_DFMO: n,q ',n,q
	write(*,*) 'init_nondom_DFMO:  bl ',bl
	write(*,*) 'init_nondom_DFMO:  bu ',bu

	alfainiz=min(10.d0,maxval(bu-bl)/10.d0)

	do i = 1,n
		alfaciniz(i) = min(10.d0,(bu(i)-bl(i))/10.d0)
	enddo
	write(*,*) 'init_nondom_DFMO:  alfa ',alfainiz,alfaciniz

	write(*,*) 'init_nondom_DFMO: printing filter ...'
	!call print_filter(n,q)
	!pause

	write(*,*) 'Building Ltilde ...'
	j = 0
	do i = 1,maxdim
		if( Ldir(i)%flag_busy ) then
			j = j+1
			do k = 1,n
				Ltilde(j,  k  ) = Ldir(i)%x(k)
				Ltilde(j,n+q+2+k) = alfaciniz(k)
			enddo
			do k = 1,q
				Ltilde(j,n+k) = Ldir(i)%f(k)
			enddo
			Ltilde(j,n+q+1  ) = 100.d0
			Ltilde(j,n+q+2  ) = alfainiz
			write(*,*) 'init_nondom_DFMO: Ltilde ',Ltilde(j,:)
		endif		
	enddo

	write(*,*) 'Merging Ltilde and Lkappa ...'
	call merge(n,q,ndim,Ltilde,Lkappa)
	
	write(*,*) 'init_nondom_DFMO: printing filter ...'
	!call print_filter(n,q)
	!pause
	
	return
end subroutine init_nondom_DFMO
