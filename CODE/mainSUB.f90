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
subroutine nuovoglobale(root,n,lbb,ubb,xbest,fbest,nf,nminloc,fglob,maxint,maxnf,	      &
			iprint,trovato,nint,mindiam,maxdiam,maxL,xmin,fmin,imin,jmin,	&
			tolglob,trigLS,first_time_call)
	use mod_type
	use mod_box
	use mod_globale
	use mod_suddividi 
	use mod_mem
	use mod_multiobj
	use paretoset
	use paretodim
	implicit none

	interface
		subroutine transfer_filter(root,n,q)
			use mod_type
			implicit none
			type(colonna),pointer	:: root
			integer			:: n, q
		end subroutine transfer_filter
		subroutine aggiorna_struttura(root,Ltilde,fdir,nint)
			use mod_type
			use mod_multiobj
			implicit none

			type(colonna),pointer		:: root
			real*8				:: Ltilde,fdir(q)
			integer			:: nint
		end subroutine aggiorna_struttura
		subroutine genera_partizione(root,n,tol,iprint)
			use mod_type
			implicit none

			type(colonna),pointer		:: root
			integer			:: n,iprint
			real*8				:: tol
		end subroutine genera_partizione
		subroutine elimina_colonna(currcol,num)
			use mod_type
			implicit none
			type(colonna),pointer		:: currcol
			integer			:: num
		end subroutine elimina_colonna
		subroutine ricintervallo_dx(root,convexhull,nconv,iprint,Ltilde,eps,fmin)
			use mod_type
			use mod_mem
			implicit none
 
			type(colonna),pointer		:: root
			type(vertice),   pointer	:: convexhull
			real*8				:: Ltilde, eps
			real*8				:: fmin
			integer			:: iprint, nconv
		end subroutine ricintervallo_dx
		subroutine ricintervallo(root,convexhull,nconv,iprint,maxL)
			use mod_type
			implicit none

			type(colonna),pointer		:: root
			type(vertice),pointer		:: convexhull
			integer			:: nconv,iprint
			real*8				:: maxL
		end subroutine ricintervallo
		subroutine ricintervallo_altre(root,convexhull,nconv,iprint,maxL)
			use mod_type
			use mod_multiobj
			implicit none

			type(colonna),pointer		:: root
			type(vertice),pointer		:: convexhull
			integer			:: nconv,iprint
			real*8				:: maxL
		end subroutine ricintervallo_altre
		subroutine ricintervallo_nnondom(root,convexhull,nconv,iprint,maxL)
			use mod_type
			use mod_multiobj
			implicit none

			type(colonna),pointer		:: root
			type(vertice),pointer		:: convexhull
			integer			:: nconv,iprint
			real*8				:: maxL
		end subroutine ricintervallo_nnondom
		subroutine ricintervallo_maxnnondom(root,convexhull,nconv,iprint,maxL)
			use mod_type
			use mod_multiobj
			implicit none

			type(colonna),pointer		:: root
			type(vertice),pointer		:: convexhull
			integer			:: nconv,iprint
			real*8				:: maxL
		end subroutine ricintervallo_maxnnondom
		subroutine riduciconvexhull(convexhull,nelim,eps,toldiam,fmin)
			use mod_type
			implicit none

			type(vertice),pointer		:: convexhull
			integer			:: nelim
			real*8				:: eps, fmin, toldiam
		end subroutine riduciconvexhull
		subroutine riduciconvexhull_altre(convexhull,nelim,eps,toldiam,fmin)
			use mod_type
			use mod_multiobj
			implicit none

			type(vertice),pointer		:: convexhull
			integer			:: nelim
			real*8				:: eps, fmin, toldiam
		end subroutine riduciconvexhull_altre
		subroutine suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
			use mod_type
			use mod_suddividi
			use mod_multiobj
			implicit none

			type(vertice), pointer	:: currch
			type(colonna),    pointer	:: root
			integer			:: n, nf, nint, cont
			real*8				:: xdir(n), fdir(q), xdir_unscaled(n)
			real*8				:: maxL,Ltilde,tol
		end subroutine suddividi
		subroutine stampa_intervalli(n,root)
			use mod_type
			implicit none

			type(colonna),   pointer	:: root
			integer			:: n
		end subroutine stampa_intervalli
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8				:: tol
		end subroutine find_colonna
	end interface

	type(colonna),    pointer		:: root, currcol, tempcol
	type(intervallo), target		:: primo
	type(intervallo), pointer		:: curr, temp
	type(vertice),    pointer		:: convexhull, currch, currch1, currchmax
	type(vertice),	  pointer		:: filter
	type(fpunt)				:: ottimo

	external				:: funct

	integer				:: n, iprint
	real*8				:: tolglob, trigLS
 
	integer				:: nf, ng, nint, num, i, nelim, nconv, k ,ktot, maxiter, numnf, numng
	integer				:: maxnf, maxint, iexit, nminloc, cont, maxcont
	logical				:: halt, direct_puro, trovato, lista_vuota, minric!minric=true fatte min loc per ric
	logical                         :: vicino, flag, find_in_cache, intorno, first_time_call
	real*8				:: norma, maxdiam, mindiam, toldiam, basdiam, fglob
	real*8				:: eps
	real*8				:: lbb(n),ubb(n)
	real*8				:: xbest(n), fbest, ftemp, tmpder, minder, maxL, Ltilde, gg(n), bl(n), bu(n), blsd(n), busd(n)
	real*8				:: xdir(n),  xx(n), ff, fdir(q), tol, xdir_unscaled(n), fminloc,  xtemp_sc(n)
	real*8				:: xmin(n,100), fmin(100), fdim(q+1)
	real*8				:: fobs(q)
	real*8				:: Lr, epsglob, tau
	integer				:: imin, jmin
	integer				:: seed(2), istop
	real				:: rr		

	! parametri ed inizializzazione

	cont        = 0
	maxcont     = 0
	tol         = 1.d-12
	toldiam     = 1.d+1*sqrt(dble(n))/2.d0 
	basdiam     = 0.0d0
	eps         = 1.d-4
	minder      = 0.0d0
	maxL        = 0.0d0
	Ltilde      = 1.d+60
	halt        = .false.
	trovato     = .false.
	lista_vuota = .false.
	memerror    = .false.
	minric      = .false.
	nf          = 0
	i_obj       = 0
	!questo è il max numero di iterazioni del locale, attenzione è nel modulo!

	nf = 0

	if(.not.associated(root)) then
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

		call funct(n,xtemp,root%int%fint,root%int%fobs)

		call update_point(n,q,xtemp,root%int%fobs,root%int,flag)
		fdim(1:q) = root%int%fobs
		fdim(q+1) = -root%int%diam
		call update_point_dim(n,q,xtemp,fdim,root%int,flag)
		nf          = 1
		nint        = 1
	else
		write(*,*) 'root pointer gia associato'
	endif

	fdir        = root%int%fobs
	xdir        = root%int%cent

	Lr          = 1.d0
	tau         = 0.05d0
	epsglob     = 0.1d0*tau*dble(n)

	ng          = 0
	nconv       = 1
	nelim       = 0
	nminloc     = 0

	xdir_unscaled = xbest

	direct_puro = .true.

	do while (.not.halt) 

		!if( (nf >= maxnf*0.1).or.first_time_call ) then 
		if(.true.) then
			i_obj = 2
			!i_obj = i_obj + 1
			!if(i_obj > q+1) i_obj = 1
		else
			!deve usare il criterio 1/nnondom
			i_obj = 0
		endif

		cont    = cont + 1
		maxcont = maxcont + 1

		if(iprint > 0) then
			write(*,*) 'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo'
			write(*,*) '  id          diam           fint'
		endif
		
		currcol => root

		do while (associated(currcol))
			if(currcol%diam < basdiam) then
				!------------------------------------------------
				! Elimina tutta la colonna corrispondente
				!------------------------------------------------
				call elimina_colonna(currcol,num)
				nint = nint - num
				if(.not.associated(currcol%next)) then
					lista_vuota = .true.
					deallocate(currcol)
					nullify(currcol)
					nullify(root)
					exit
				else
					currcol => currcol%next
					deallocate(root)
					root => currcol
					nullify(root%pred)
				endif
			else
				if(.not.associated(currcol%next)) then
					maxdiam = currcol%diam
					! non ci vuole un exit?					
					! exit
				endif
				curr => currcol%int
				do while (associated(curr))
					call unscalevars(n,curr%cent,ytemp,xbar,lbs,ubs)
					call unscalevars_direct(n,ytemp,xtemp)
					fdim(1:q) = curr%fobs
					fdim(q+1) = -curr%diam
					if(intorno_dim(n,q,xtemp,fdim,currcol%diam)) then
						call update_point_dim(n,q,xtemp,fdim,curr,flag)
					endif
					curr => curr%next
				enddo
				currcol => currcol%next
			endif
		enddo

!------------------------------------------------------------------------------

		if(tau*maxdiam**2.d0 < epsglob) epsglob = tau*maxdiam**2.d0

		if(associated(root)) then
			nullify(root%pred)
			if((.not.associated(root%int)).and.(.not.associated(root%next))) lista_vuota = .true.
		else
			lista_vuota = .true.
		endif
		if(lista_vuota) then
			halt = .true.
			write(*,*) 'lista vuota'
			exit
		endif

		mindiam = root%diam

		!--------------------------
		! stopping criterion
		!--------------------------

		trovato = (abs(fbest-fglob)/dmax1(1.0d0,abs(fglob)) < tolglob)
		if((nf >= maxnf).or.(nint >= maxint).or.trovato) then
			halt = .true.
			cycle
		endif

		!call stampa_intervalli(n,root)
		!call print_filter(n,q)
		!pause

		!if( (nf < maxnf*0.1).and..not.first_time_call ) then
		if( .false. ) then
			!----------------------------------------------------
			! divide gli intervalli usando il crit. 1/nnondom
			!----------------------------------------------------
			if(iprint > 0) write(*,*) 'estrae gli intervalli usando 1/nnondom'
			nelim = 0
			!call ricintervallo_nnondom(root,convexhull,nconv,iprint,Ltilde)
			call ricintervallo_maxnnondom(root,convexhull,nconv,iprint,Ltilde)
			currch => convexhull

		else
			if(i_obj == q+1) then
			!if(.true.) then
				!---------------------------------------
				! divide gli intervalli nel paretoset
				!---------------------------------------
				if(iprint > 0) write(*,*) 'estrae gli intervalli dei non-dominati'
				nelim = 0
				!call ricintervallo_pareto(root,convexhull,nconv,iprint,Ltilde)
				call ricintervallo_paretodim(n,q,root,convexhull,nconv,iprint,maxL,maxdiam)
				currch => convexhull
			else
				!---------------------------------------
				! ric. intervallo potenzialmente ottimo
				!---------------------------------------

				if(iprint > 0) write(*,*) 'inizio ric. potenzialmente ottimi'

				if(i_obj == 1) then
					!----------------------------------------------------
					! chiama la vecchia subroutine perche' la struttura
					! e' ordinata rispetto alla prima f.ob.
					!----------------------------------------------------
					call ricintervallo       (root,convexhull,nconv,iprint,Ltilde)
				else
					call ricintervallo_altre (root,convexhull,nconv,iprint,Ltilde)
				endif

				if (memerror) then
					write(*,*) '  fine memoria disponibile'
					halt = .true.
				endif
				if(iprint > 0) write(*,*) '  fine ric. potenzialmente ottimi'

				!write(*,*) 'nconv=',nconv

				!----------------------------------------------
				! riduci il convex hull con il criterio su
				! fmin
				!----------------------------------------------

				if (direct_puro) toldiam = 0.0d0
		
				if(iprint > 0) write(*,*) 'inizio riduzione convex hull'

				if(i_obj == 1) then
					call riduciconvexhull(convexhull,nelim,eps,toldiam,fdir(1))
				else
					call riduciconvexhull_altre(convexhull,nelim,eps,toldiam,fdir(i_obj))
				endif

				if(iprint > 0) write(*,*) '  fine riduzione convex hull'

				currch => convexhull
		
				if(iprint > 0) write(*,*) 'nconv = ',nconv,' nelim = ',nelim

				do i = 1,nelim
					currch => currch%next
				enddo

			endif
		endif
		!-----------------------------------------------------
		! Se non hai esaurito il budget di calcoli di funzione
		! prosegui DIRECT con le suddivisioni
		!-----------------------------------------------------

		!-----------------------------------------------
		! rimuovo i rimanenti intervalli sul convexhull
		! dalla struttura dati per inserirli in seguito
		! nella posizione corretta
		!-----------------------------------------------
		if(iprint > 0) then
			write(*,*)  'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			write(*,*)  '  id          diam           fint'
			write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			write(24,*) '     diam           fint',fdir
		endif

		currch1 => currch
		do while (associated(currch1))
			if(iprint > 0) then
				write(*,*) currch1%int%id, currch1%int%diam, currch1%int%fobs
				write(24,*) currch1%int%diam, currch1%int%fobs
			endif

			call find_colonna(root,currch1%int,tol,tempcol)
			if(tempcol%int%id == currch1%int%id) then
				!------------------------------------------------------------------
				! l'intervallo corrente sul convex hull e' il primo della colonna
				!------------------------------------------------------------------
				if(associated(currch1%int%next)) then
					tempcol%int => currch1%int%next
					nullify(tempcol%int%pred)
				else
					!write(*,*) 'la colonna si e'' svuotata quindi la elimino'
					!--------------------------------
					! la colonna si e' svuotata
					! quindi la elimino
					!--------------------------------
					if((.not.associated(tempcol%pred)) .and. (.not.associated(tempcol%next))) then
						!--------------------------------
						! E' l'unica colonna
						!--------------------------------
						nullify(tempcol)
						deallocate(root)
						nullify(root)
						!write(*,*) 'la col e'' unica e la elimino',associated(root)
						exit
					elseif(.not.associated(tempcol%pred)) then
						!write(*,*) 'la col e'' la prima ma non unica'
						!--------------------------------
						! E' la prima ma non l'unica
						!--------------------------------
						tempcol => root
						root => root%next
						deallocate(tempcol)
						nullify(root%pred)
					elseif(.not.associated(tempcol%next)) then
						!--------------------------------
						! E' l'ultima ma non l'unica
						!--------------------------------
						!write(*,*) 'la col e'' l''ultima ma non unica',tempcol%int%id,root%int%id
						nullify(tempcol%pred%next)
						deallocate(tempcol)
					else
						!--------------------------------
						! E' in mezzo
						!--------------------------------
						!write(*,*) 'la col e'' in mezzo',tempcol%int%id,root%int%id
						tempcol%pred%next => tempcol%next
						tempcol%next%pred => tempcol%pred
						deallocate(tempcol)
					endif 
				endif
			else
				!-------------------------------------------------------------------------
				! l'intervallo corrente sul convex hull NON e' il primo della colonna
				! quindi i_obj > 1 (NON PUO' ESSERE LA PRIMA F.OB.)
				!-------------------------------------------------------------------------
				!write(*,*) 'shit shit !! il primo non e'' il primo',tempcol%int%id,currch1%int%id
				!stop
				if(associated(currch1%int%next).and.associated(currch1%int%pred)) then
					!write(*,*) '1'
					!write(*,*) 'non e'' l''ultimo'
					!------------------------------------------
					! l'intervallo NON e' l'ultimo
					!------------------------------------------
					currch1%int%next%pred => currch1%int%pred
					currch1%int%pred%next => currch1%int%next
					nullify(currch1%int%pred)
					nullify(currch1%int%next)
				elseif(associated(currch1%int%pred)) then
					!write(*,*) '3'
					!write(*,*) 'e'' l''ultimo'
					nullify(currch1%int%pred%next)
					nullify(currch1%int%pred)
				elseif(associated(currch1%int%next)) then
					!write(*,*) '2'
					!write(*,*) 'e'' l''ultimo'
					nullify(currch1%int%next%pred)
					nullify(currch1%int%next)
				endif
			endif
			currch1 => currch1%next
		enddo

		if(iprint > 0) then
			write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		endif


		if(iprint > 0) write(*,*) 'ciclo principale: dopo esegui'

		if(associated(root)) then
			nullify(root%pred)
			if((.not.associated(root%int)).and.(.not.associated(root%next))) lista_vuota = .true.
		else
			lista_vuota = .true.
		endif
		if(lista_vuota) then
			write(*,*) 'lista vuota'
			lista_vuota = .false.
		endif

		if(iprint > 0) write(*,*) 'inizio suddivisioni'

		do while (associated(currch).and. (.not. memerror))
			!write(*,*) 'mainsub: ',currch%int%id, currch%int%flagdiv
			!----------------------------------------------------
			! suddivisione dell'intervallo potenzialmente ottimo
			!----------------------------------------------------
			if(currch%int%flagdiv) call suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
!			if(currch%int%flagdiv) call suddividi(currch,root,n,nf,nint,xdir,fbest,xdir_unscaled,maxL,Ltilde,cont,tol)
			if(cont == 0) maxcont = 0
			if (memerror) then
				write(*,*)'fine memoria, esco al momento della suddivisione '
				halt = .true.
				exit
			endif

			currch => currch%next
		enddo

		if(iprint > 0) write(*,*) '  fine suddivisioni'

		if(iprint > 0) write(*,*) 'ciclo principale: dopo suddividi'

		!-------------------------------------------
		! dealloca la lista che memorizza gli
		! intervalli sul convex hull
		!-------------------------------------------

		currch => convexhull
		do while (associated(currch))
			!if(associated(currch%int)) currch%int%flagcon = .false.
			currch => currch%next
			nullify(convexhull%int)
			nullify(convexhull%next)
			deallocate(convexhull)
			convexhull => currch
		enddo

		if(iprint > 0) write(*,*) '----- ',num_el_L

		if(associated(root)) then
			nullify(root%pred)
			if((.not.associated(root%int)).and.(.not.associated(root%next))) lista_vuota = .true.
		else
			lista_vuota = .true.
		endif
		if(lista_vuota) then
			halt = .true.
			write(*,*) 'lista vuota'
			exit
		endif

		write(*,800) nf,nelim,nconv-nelim,mindiam,maxdiam
		write(*,810) maxL*maxval(ubb-lbb),Ltilde,cont,nf-nint,nint

		!call stampa_intervalli(n,root)
	enddo

	return

798		FORMAT('fdir = ',$)
799		FORMAT(es11.4,1x,$)
800		FORMAT(' nf=',i11,   '  nelim=',i11,' ncnv-nelim=',i8,' diam=',es11.4,' DIAM=',es11.4)
810		FORMAT('  maxL=',es11.4,'  Ltilde=',es11.4,'    cont=',I11,   ' num_el=',I11,    '   nint=',i11)
830		FORMAT(' fbest=',es11.4)
end subroutine nuovoglobale
