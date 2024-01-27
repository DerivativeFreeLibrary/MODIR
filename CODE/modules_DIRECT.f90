module mod_type
	type typ_lista
		real*8, allocatable		:: fobs(:), x(:)
		type(typ_lista), pointer	:: prev, next
	end type typ_lista
	type intervallo
		real*8, allocatable		:: cent(:), dimen(:), lbs(:), ubs(:), xbars(:)
		real*8, allocatable		:: fobs(:)
		real*8				:: fint, diam, maxdim, der
		logical				:: flagloc, flagdiv, flagcon, flagopt
		integer				:: id,nnondom
		type(typ_lista), pointer	:: list_nondom
		type(intervallo), pointer	:: next, pred
	end type intervallo
	type vertice
		type(intervallo), pointer	:: int
		type(vertice),    pointer	:: next
	end type vertice
	type colonna
		real*8				:: diam
		type(colonna),    pointer	:: next, pred
		type(intervallo), pointer	:: int
	end type colonna
	type fpunt
		real*8				:: f
		type(intervallo), pointer	:: punt
	end type fpunt
end module mod_type

module mod_globale
	real*8, allocatable			:: globxbest(:)
	real*8					:: globfbest
	integer					:: globnftot, globnf, globmaxnf
end module mod_globale

module mod_box
	real*8, allocatable			:: lb(:), ub(:), xbar(:)
	real*8, allocatable			:: lbs(:), ubs(:)
	real*8, allocatable			:: xtemp(:), ytemp(:)
	real*8					:: ampiezza, fattore  !metà dell'ampiezza
end module mod_box

module mod_suddividi
	real*8,  allocatable			:: vetf1(:), vetf2(:)
	real*8,  allocatable			:: matf1(:,:), matf2(:,:)
	real*8,  allocatable			:: xsud(:), ysud(:)
	logical, allocatable			:: mask(:)
end module mod_suddividi

module mod_mem
	logical					:: memerror
	integer					:: num_el_L
end module mod_mem

module mod_multiobj
	!---------------------------------------------------------------------------
	! q     e' il numero di funzioni obiettivo
	!
	! i_obj dice l'indice della f.ob. rispetto a cui determinare
	!       gli iperintervalli potenzialmente ottimi
	!---------------------------------------------------------------------------
	integer					:: q
	integer					:: i_obj
end module mod_multiobj

module paretoset
	use mod_type
	use mod_mem
	use mod_multiobj
	integer, parameter :: maxdim = 20000
	!--------------------------------------------------------------------------------------------------
	! il filtro e' un array la cui generica riga di dimensione n+q+1+1+1 
	! e' intesa come segue:
	!
	!   ------------------------------------------------------------------------------------
	!   | x(1)  ....   x(n) | f(1) .... f(q) f(q+1) | flag_busy | puntatore_all_intervallo |
	!   ------------------------------------------------------------------------------------
	! maxdim indica il numero di righe nella matrice
	! flag_busy = 100 se la riga corrispondente e' occupata da un non-dominato
	!           = 0   se la riga corrispondente e' vuota/libera
	! puntatore_all_intervallo è il puntatore all'intervallo che nella
	!	struttura dati di MODIRECT *contiene* il punto x
	! N.B. il punto x potrebbe ed anzi spesso NON è il centroide dell'intervallo
	!	puntato da puntatore_all_intervallo. Questo perché i non-dominati
	!	possono essere stati ereditati dalla lista su cui opera DFMO
	! 	cioè l'algoritmo LOCALE multiobiettivo
	!
	! x è il punto non-dominato
	! f è il vettore di funzioni obiettivo e risulta: f(1:q) = F(x)
	! NO    f(q+1) = ampiezza dell'intervallo puntato da puntatore_all_intervallo quindi
	! NO	f(q+1) = punt%diam
	! ATTENZIONE: quando si fanno le divisioni e quindi quando cambiano i diametri degli intervalli
	!	è necessario ricontrollare la dominanza dei punti dentro Ldir !!!!
	!
	!--------------------------------------------------------------------------------------------------
	type paretoset_elem
		real*8, allocatable		:: x(:), f(:)
		logical				:: flag_busy
		type(intervallo), pointer	:: punt
	end type paretoset_elem

	type(paretoset_elem) 			:: Ldir(maxdim)

	contains 

	subroutine init_paretoset(n,q)
		implicit none
		integer				:: n, q
		integer				:: i

		do i = 1,maxdim
			allocate(Ldir(i)%x(n))
			allocate(Ldir(i)%f(q))
			Ldir(i)%flag_busy = .false.
			nullify(Ldir(i)%punt)
		enddo

		return
	end subroutine init_paretoset

	subroutine canc_paretoset()
		implicit none
		integer				:: i

		do i = 1,maxdim
			deallocate(Ldir(i)%x)
			deallocate(Ldir(i)%f)
		enddo

		return
	end subroutine canc_paretoset

	logical function domina(q,f1,f2)
		implicit none
		integer q,i
		real*8 f1(q),f2(q)
		!---------------------------------------------------
		! return TRUE   if f1 domina (sec. Pareto, <=_P) f2
		!        FALSE  otherwise
		!---------------------------------------------------
		domina = .true.
		do i = 1,q
			if(f1(i) > f2(i)) then
				domina = .false.
				exit
			endif
		enddo
		if(.not.domina) then
			!write(*,*) 'domina: f1 = ',f1
			!write(*,*) 'domina: f2 = ',f2
			!write(*,*) 'domina = ',domina
			return
		endif
		domina = .false.
		do i = 1,q
			if(f1(i) < f2(i)) then
				domina = .true.
				exit
			endif
		enddo

		!write(*,*) 'domina: f1 = ',f1
		!write(*,*) 'domina: f2 = ',f2
		!write(*,*) 'domina = ',domina

		return
	end function domina

	logical function minore_uguale(q,f1,f2)
		implicit none
		integer q,i
		real*8 f1(q),f2(q)
		!-----------------------------------
		! return TRUE   if f1 <= f2
		!        FALSE  otherwise
		!-----------------------------------
		minore_uguale = .true.
		do i = 1,q
			if(f1(i) > f2(i)) then
				minore_uguale = .false.
				exit
			endif
		enddo

		return
	end function minore_uguale

	subroutine update_point(n,q,x,f,punt,flag_change)
		implicit none
		integer 			:: n,q
		real*8 				:: x(n), f(q)
		type(intervallo), pointer 	:: punt
		integer 			:: i, j, indfree
		integer 			:: flag
		logical 			:: flag_change

		!--------------------------------------------------------------
		! flag assume tre valori: 1,2,3
		! flag = 1 : Ldir(i,:) < (x,f) 
		! flag = 2 : (x,f) < Ldir(i,:)
		! flag = 3 : (x,f) non domina e non e' dominato Ldir(i,:)
		!--------------------------------------------------------------
		flag_change = .false.
		indfree = -1
		flag = 0
		!write(*,*) 'update_point: f = ',f
		do i = 1,maxdim
			flag = 0
			if((.not.Ldir(i)%flag_busy).and.(indfree == -1)) then
				indfree = i
			elseif(Ldir(i)%flag_busy) then
				if( (minore_uguale(q,Ldir(i)%f,f)) .and. &
				    (minore_uguale(q,f,Ldir(i)%f)) ) then
					flag = 1
					exit
				endif
				!write(*,*) '    i = ',i,' flag_busy = ',Ldir(i)%flag_busy,' f = ',Ldir(i)%f
				if(domina(q,Ldir(i)%f,f)) then
					flag = 1
				elseif(domina(q,f,Ldir(i)%f)) then
					flag = 2
				else
					flag = 3
				endif
				!write(*,*) '    flag = ',flag
				select case (flag)
				case(1)
					!------------------------------------------------
					!il punto passato e' dominato da uno nel filtro
					!------------------------------------------------
					exit
				case(2)
					!------------------------------------------------
					!il punto passato domina uno nel filtro
					!------------------------------------------------
					Ldir(i)%flag_busy = .false.
					if(indfree == -1) indfree = i
				case(3)
				end select
			endif
		enddo
		!write(*,*) 'update_point:    flag = ',flag,' indfree = ',indfree
		if(flag /= 1) then
			flag_change = .true.
			if(indfree /= -1) then
				Ldir(indfree)%x = x
				Ldir(indfree)%f = f
				Ldir(indfree)%flag_busy = .true.
				Ldir(indfree)%punt => punt
				!write(*,*) 'update_point: Ldir%f = ',Ldir(indfree)%f
			else
				write(*,*) 'WARNING: Filtro pieno'
				!pause
			endif
		endif
		!write(*,*) 'update_point: -----------------------------------'
		return
	end subroutine update_point

	subroutine ricintervallo_pareto(root,convexhull,nconv,iprint,maxL)
		implicit none

		type(colonna),pointer		:: root, tempcol
		type(intervallo), pointer	:: curr
		type(vertice),pointer		:: convexhull, currch
		integer				:: nconv, iprint, i, i1, sv
		integer				:: counter
		real*8				:: maxL
		
		allocate(convexhull,stat=sv)
		if(sv.ne.0) then
			memerror = .true.
			return
		endif
		! cerca il primo non-dominato
		i1 = 0
		nconv = 0
		do i = 1,maxdim
			if((Ldir(i)%flag_busy).and.(Ldir(i)%punt%diam > 1.d-2).and.(.not.Ldir(i)%punt%flagcon)) then
				i1 = i
				exit
			endif
		enddo

		if(i1 > 0) then
			convexhull%int => Ldir(i1)%punt
			Ldir(i1)%punt%flagcon = .true.
		    	currch => convexhull
			nconv = 1

			if(iprint > 0) then
				write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
				write(*,*) '         i_obj = ', i_obj
				write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
				write(*,*) '  id          diam           fint'
				write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fobs
			endif

			do i = i1+1,maxdim
				if((Ldir(i)%flag_busy).and.(Ldir(i)%punt%diam > 1.d-2).and.(.not.Ldir(i)%punt%flagcon)) then
					allocate(currch%next, stat = sv)
					if(sv.ne.0) then
						memerror = .true.
						return
					endif
					currch%next%int => Ldir(i)%punt
					Ldir(i)%punt%flagcon = .true.
					currch => currch%next
					nconv  = nconv + 1
					if (iprint > 0)	then
						write(*,*) currch%int%id, currch%int%diam, currch%int%fobs
					endif
				endif
			enddo
		endif

		tempcol => root

		do while(associated(tempcol%next))
			tempcol => tempcol%next
		enddo
		
		if(nconv <= 0) then
			counter = 0
			curr => tempcol%int
			!write(*,*) '--------- QUI ----------'
			!do while ((counter <= 10*q).and.(associated(curr)))
			do while ((counter <= q).and.(associated(curr)))
				!write(*,*) counter, q
				if(.not.curr%flagcon) then
					if((counter == 0).and.(i1 == 0)) then
						allocate(convexhull,stat=sv)
						if(sv.ne.0) then
							memerror = .true.
							return
						endif
						convexhull%int => curr
						curr%flagcon = .true.
							currch => convexhull
						nconv = 1
						counter = counter+1
					else
						allocate(currch%next, stat = sv)
						if(sv.ne.0) then
							memerror = .true.
							return
						endif
						currch%next%int => curr
						curr%flagcon = .true.
						currch => currch%next
						nconv  = nconv + 1
						counter = counter+1
					endif
					if (iprint > 0)	then
						write(*,*) currch%int%id, currch%int%diam, currch%int%fobs
					endif
				endif
				curr => curr%next
			enddo
		endif

		if(iprint > 0) write(*,*) 'i1 = ',i1,' nconv = ',nconv
		if((i1==0).and.(nconv==0)) then
			deallocate(convexhull)
			nullify(convexhull)
		else
			nullify(currch%next)
		endif

		if(iprint > 0) then
			write(*,*) ' nconv = ',nconv
			write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			!pause
		endif

		return
	end subroutine ricintervallo_pareto

	subroutine print_filter(n,q)
		implicit none
		integer n,q,i,j

		!open(33,file='fort.44',status='replace')
		do i = 1,maxdim
			if(Ldir(i)%flag_busy) then
				write(*,100) i
				write(*,200)
				do j = 1,n
					!write(33,110) Ldir(i)%x(j)
					write(* ,110) Ldir(i)%x(j)
				enddo
				write(*,300)
				do j = 1,q
					!write(33,110) Ldir(i)%f(j)
					write(* ,110) Ldir(i)%f(j)
				enddo
				write(*,*)
				!write(33,*)
			endif
		enddo
		!close(33)
		return

	100 format(1x,i7,$)
        200 format(1x,'x=',$)
        300 format(1x,'f=',$)
	110 format(1x,es20.13,$)
	end subroutine print_filter
end module paretoset

module paretodim
	use mod_type
	use mod_mem
	use mod_multiobj
	integer, parameter :: maxdim = 20000
	!--------------------------------------------------------------------------------------------------
	! il filtro e' un array la cui generica riga di dimensione n+q+1+1+1 
	! e' intesa come segue:
	!
	!   ------------------------------------------------------------------------------------
	!   | x(1)  ....   x(n) | f(1) .... f(q) f(q+1) | flag_busy | puntatore_all_intervallo |
	!   ------------------------------------------------------------------------------------
	! maxdim indica il numero di righe nella matrice
	! flag_busy = 100 se la riga corrispondente e' occupata da un non-dominato
	!           = 0   se la riga corrispondente e' vuota/libera
	! puntatore_all_intervallo è il puntatore all'intervallo che nella
	!	struttura dati di MODIRECT *contiene* il punto x
	! N.B. il punto x potrebbe ed anzi spesso NON è il centroide dell'intervallo
	!	puntato da puntatore_all_intervallo. Questo perché i non-dominati
	!	possono essere stati ereditati dalla lista su cui opera DFMO
	! 	cioè l'algoritmo LOCALE multiobiettivo
	!
	! x è il punto non-dominato
	! f è il vettore di funzioni obiettivo e risulta: f(1:q) = F(x)
	! NO    f(q+1) = ampiezza dell'intervallo puntato da puntatore_all_intervallo quindi
	! NO	f(q+1) = punt%diam
	! ATTENZIONE: quando si fanno le divisioni e quindi quando cambiano i diametri degli intervalli
	!	è necessario ricontrollare la dominanza dei punti dentro Ldir !!!!
	!
	!--------------------------------------------------------------------------------------------------
	type paretodim_elem
		real*8, allocatable		:: x(:), f(:)
		logical				:: flag_busy
		type(intervallo), pointer	:: punt
	end type paretodim_elem

	type(paretodim_elem) 			:: Ldim(maxdim)

	contains 

	subroutine init_paretodim(n,q)
		implicit none
		integer				:: n, q
		integer				:: i

		do i = 1,maxdim
			allocate(Ldim(i)%x(n))
			allocate(Ldim(i)%f(q+1))
			Ldim(i)%flag_busy = .false.
			nullify(Ldim(i)%punt)
		enddo

		return
	end subroutine init_paretodim

	subroutine canc_paretodim()
		implicit none
		integer				:: i

		do i = 1,maxdim
			deallocate(Ldim(i)%x)
			deallocate(Ldim(i)%f)
		enddo

		return
	end subroutine canc_paretodim

	logical function intorno_dim(n,q,x,f,diam)
		implicit none
		integer n,q,i
		real*8 x(n), f(q), delta, diam
		real*8 theta

		theta = 1.d-1
		intorno_dim = .true.
		do i = 1,maxdim
			if(Ldim(i)%flag_busy) then
				if( maxval(abs(f-theta*diam-Ldim(i)%f)) < 1.d-4*maxval(abs(Ldim(i)%f)) ) then
					intorno_dim = .false.
					exit
				endif
			endif
		enddo

		return
	end function intorno_dim

	subroutine aggiorna_dim(n,q,punt)
		implicit none
		integer						:: n, q
		type(intervallo), pointer	:: punt
		integer						:: i
		logical						:: flag
		double precision			:: xtemp(n), fdim(q+1)

		do i = 1,maxdim
			if(Ldim(i)%flag_busy) then
				if(Ldim(i)%punt%id == punt%id) then
					Ldim(i)%flag_busy = .false.
					xtemp = Ldim(i)%x
					fdim  = Ldim(i)%f
					fdim(q+1) = -punt%diam
					call update_point_dim(n,q,xtemp,fdim,punt,flag)
				endif
			endif
		enddo

	end subroutine aggiorna_dim

	logical function domina_dim(q,f1,f2)
		implicit none
		integer q,i
		real*8 f1(q+1),f2(q+1)
		!---------------------------------------------------
		! return TRUE   if f1 domina (sec. Pareto, <=_P) f2
		!        FALSE  otherwise
		!---------------------------------------------------
		domina_dim = .true.
		do i = 1,q+1
			if(f1(i) > f2(i)) then
				domina_dim = .false.
				exit
			endif
		enddo
		if(.not.domina_dim) then
			!write(*,*) 'domina: f1 = ',f1
			!write(*,*) 'domina: f2 = ',f2
			!write(*,*) 'domina = ',domina
			return
		endif
		domina_dim = .false.
		do i = 1,q+1
			if(f1(i) < f2(i)) then
				domina_dim = .true.
				exit
			endif
		enddo

		!write(*,*) 'domina: f1 = ',f1
		!write(*,*) 'domina: f2 = ',f2
		!write(*,*) 'domina = ',domina

		return
	end function domina_dim

	logical function minore_uguale_dim(q,f1,f2)
		implicit none
		integer q,i
		real*8 f1(q+1),f2(q+1)
		!-----------------------------------
		! return TRUE   if f1 <= f2
		!        FALSE  otherwise
		!-----------------------------------
		minore_uguale_dim = .true.
		do i = 1,q+1
			if(f1(i) > f2(i)) then
				minore_uguale_dim = .false.
				exit
			endif
		enddo

		return
	end function minore_uguale_dim

	subroutine update_point_dim(n,q,x,f,punt,flag_change)
		implicit none
		integer 			:: n,q
		real*8 				:: x(n), f(q+1)
		type(intervallo), pointer 	:: punt
		integer 			:: i, j, indfree
		integer 			:: flag
		logical 			:: flag_change

		!--------------------------------------------------------------
		! flag assume tre valori: 1,2,3
		! flag = 1 : Ldim(i,:) < (x,f) 
		! flag = 2 : (x,f) < Ldim(i,:)
		! flag = 3 : (x,f) non domina e non e' dominato Ldim(i,:)
		!--------------------------------------------------------------
		flag_change = .false.
		indfree = -1
		flag = 0
		!write(*,*) 'update_point: f = ',f
		do i = 1,maxdim
			flag = 0
			if((.not.Ldim(i)%flag_busy).and.(indfree == -1)) then
				indfree = i
			elseif(Ldim(i)%flag_busy) then
				if( (minore_uguale_dim(q,Ldim(i)%f,f)) .and. &
				    (minore_uguale_dim(q,f,Ldim(i)%f)) ) then
					flag = 1
					exit
				endif
				!write(*,*) '    i = ',i,' flag_busy = ',Ldim(i)%flag_busy,' f = ',Ldim(i)%f
				if(domina_dim(q,Ldim(i)%f,f)) then
					flag = 1
				elseif(domina_dim(q,f,Ldim(i)%f)) then
					flag = 2
				else
					flag = 3
				endif
				!write(*,*) '    flag = ',flag
				select case (flag)
				case(1)
					!------------------------------------------------
					!il punto passato e' dominato da uno nel filtro
					!------------------------------------------------
					exit
				case(2)
					!------------------------------------------------
					!il punto passato domina uno nel filtro
					!------------------------------------------------
					Ldim(i)%flag_busy = .false.
					if(indfree == -1) indfree = i
				case(3)
				end select
			endif
		enddo
		!write(*,*) 'update_point:    flag = ',flag,' indfree = ',indfree
		if(flag /= 1) then
			flag_change = .true.
			if(indfree /= -1) then
				Ldim(indfree)%x = x
				Ldim(indfree)%f = f
				Ldim(indfree)%flag_busy = .true.
				Ldim(indfree)%punt => punt
				!write(*,*) 'update_point: Ldim%f = ',Ldim(indfree)%f
			else
				write(*,*) 'WARNING: Filtro pieno'
				!pause
			endif
		endif
		!write(*,*) 'update_point: -----------------------------------'
		return
	end subroutine update_point_dim

	subroutine ricintervallo_paretodim(n,q,root,convexhull,nconv,iprint,maxL,maxdiam)
		implicit none

		type(colonna),pointer		:: root, tempcol
		type(intervallo), pointer	:: curr
		type(vertice),pointer		:: convexhull, currch
		integer				:: n, q
		integer				:: nconv, iprint, i, j, i1, i2, sv
		integer				:: counter
		real*8				:: maxL, maxdiam
		real*8				:: app(maxdim)
		integer				:: ipt(maxdim)
		logical				:: flag

		allocate(convexhull,stat=sv)
		if(sv.ne.0) then
			memerror = .true.
			return
		endif

		call spread_ord(n,q,app,ipt,flag)

		! cerca il primo non-dominato
		i1 = 0
		nconv = 0
		do j = 1,maxdim
			i = ipt(j)
			if((Ldim(i)%flag_busy).and.(Ldim(i)%punt%diam > 1.d-1*maxdiam).and.(.not.Ldim(i)%punt%flagcon)) then
			!if((Ldim(i)%flag_busy).and.(.not.Ldim(i)%punt%flagcon)) then
				i1 = i
				i2 = j
				exit
			endif
		enddo

		!write(*,*) 'i1 = ',i1,ipt(1),ipt(2),app(1),app(2)
		!pause

		if(i1 > 0) then
			convexhull%int => Ldim(i1)%punt
			Ldim(i1)%punt%flagcon = .true.
		    	currch => convexhull

			nconv = 1

			if(iprint > 0) then
				write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
				write(*,*) '         i_obj = ', i_obj
				write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
				write(*,*) '  id          diam           fint'
				write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fobs
			endif

			do j = i2+1,maxdim
				i = ipt(j)
				if((Ldim(i)%flag_busy).and.(Ldim(i)%punt%diam >  1.d-1*maxdiam).and.(.not.Ldim(i)%punt%flagcon)) then
				!if((Ldim(i)%flag_busy).and.(.not.Ldim(i)%punt%flagcon)) then
					allocate(currch%next, stat = sv)
					if(sv.ne.0) then
						memerror = .true.
						return
					endif
					currch%next%int => Ldim(i)%punt
					Ldim(i)%punt%flagcon = .true.
					currch => currch%next
					nconv  = nconv + 1
					if (iprint > 0)	then
						write(*,*) currch%int%id, currch%int%diam, currch%int%fobs
					endif
				endif
			enddo
		endif

		!write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ', nconv
		tempcol => root

		do while(associated(tempcol%next))
			tempcol => tempcol%next
		enddo
		
		if(nconv <= 0) then
			counter = 0
			curr => tempcol%int
			!write(*,*) '--------- QUI ----------'
			!do while ((counter <= 10*q).and.(associated(curr)))
			do while ((counter <= q).and.(associated(curr)))
				!write(*,*) counter, q
				if(.not.curr%flagcon) then
					if((counter == 0).and.(i1 == 0)) then
						allocate(convexhull,stat=sv)
						if(sv.ne.0) then
							memerror = .true.
							return
						endif
						convexhull%int => curr
						curr%flagcon = .true.
							currch => convexhull
						nconv = 1
						counter = counter+1
					else
						allocate(currch%next, stat = sv)
						if(sv.ne.0) then
							memerror = .true.
							return
						endif
						currch%next%int => curr
						curr%flagcon = .true.
						currch => currch%next
						nconv  = nconv + 1
						counter = counter+1
					endif
					if (iprint > 0)	then
						write(*,*) currch%int%id, currch%int%diam, currch%int%fobs
					endif
				endif
				curr => curr%next
			enddo
		endif

		if(iprint > 0) write(*,*) 'i1 = ',i1,' nconv = ',nconv
		if((i1==0).and.(nconv==0)) then
			deallocate(convexhull)
			nullify(convexhull)
		else
			nullify(currch%next)
		endif

		if(iprint > 0) then
			write(*,*) ' nconv = ',nconv
			write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			!pause
		endif

		return
	end subroutine ricintervallo_paretodim

	subroutine print_filterdim(n,q)
		implicit none
		integer n,q,i,j

		open(33,file='fort.44',status='replace')
		do i = 1,maxdim
			if(Ldim(i)%flag_busy) then
				write(*,100) i
				write(*,200)
				do j = 1,n
					!write(33,110) Ldir(i)%x(j)
					write(* ,110) Ldim(i)%x(j)
				enddo
				write(*,300)
				do j = 1,q+1
					write(33,110) Ldim(i)%f(j)
					write(* ,110) Ldim(i)%f(j)
				enddo
				write(*,*)
				write(33,*)
			endif
		enddo
		close(33)
		return

	100 format(1x,i7,$)
        200 format(1x,'x=',$)
        300 format(1x,'f=',$)
	110 format(1x,es20.13,$)
	end subroutine print_filterdim

	subroutine spread_ord(n,q,app,ipt,flag)
		implicit none
		integer	:: n, q, i, j
		real*8		:: alfas(maxdim), fobj(maxdim,q), delta(maxdim), app(maxdim)	
		real*8		:: dummy(n+q+2+n), dummyr
		real*8		:: Lktemp(maxdim,n+q+2+n)
		integer	:: dummyi
		integer	:: ipt(maxdim), jlast, indf(maxdim), ncomp, iapp(maxdim)
		logical	:: primo,flag

!	type paretoset_elem
!		real*8, allocatable		:: x(:), f(:)
!		logical				:: flag_busy
!		type(intervallo), pointer	:: punt
!	end type paretoset_elem

		!write(*,*) 'in SPREAD'
		flag = .true.
		delta = 0.d0
		ncomp = 0
		do i = 1,maxdim
			!write(*,*) 'in SPREAD',i,ncomp
			if(Ldim(i)%flag_busy) then
				ncomp = ncomp+1
				fobj(ncomp,:) = Ldim(i)%f
				indf(ncomp) = i
				!write(*,*) fobj(ncomp,:)
			endif
		enddo

		!do i = 1,ncomp
		!	write(*,*) fobj(i,:),indf(i)
		!enddo
		!write(*,*) 'in SPREAD', ncomp
		if(ncomp.le.1) then
		!	pause
			flag = .false.
			ipt(1) = 1
			return
		endif
	
		do j = 1,q

		!	write(*,*) 'SPREAD: ordering obj ',j
			alfas = fobj(:,j)
		!	write(*,*) 'SPREAD: alfas ',alfas(1:ncomp)
		!	write(*,*) 'SPREAD:  indf ',indf(1:ncomp)
			call qsortd(alfas(1:ncomp),ipt(1:ncomp),ncomp)
		!	write(*,*) alfas(1:ncomp)
		!	pause
		!	write(*,*) 'SPREAD: ... done'
		!	write(*,*) 'SPREAD: ordering indf '
			do i=1, ncomp
				app(ncomp+1-i)=alfas(ipt(i))
				iapp(ncomp+1-i)=indf(ipt(i))
			end do
			alfas(1:ncomp) = app(1:ncomp)
			indf(1:ncomp) = iapp(1:ncomp)
		!	write(*,*) 'SPREAD: ... done'
		
			delta(indf(1)) = delta(indf(1))  + 1.d+10
			delta(indf(ncomp)) = delta(indf(ncomp)) + 1.d+10
			!delta(indf(1)) = delta(indf(1)) + abs(alfas(1)-alfas(2))
			do i = 2,ncomp-1
				delta(indf(i)) = delta(indf(i)) + abs(alfas(i-1)-alfas(i+1))/(alfas(1)-alfas(ncomp))
			enddo
			!delta(indf(ncomp)) = delta(indf(ncomp)) + abs(alfas(ncomp-1)-alfas(ncomp))
		enddo

		delta = delta/dble(ncomp)
		do i = 1,maxdim
			!write(*,*) 'in SPREAD',i,ncomp
			if(Ldim(i)%flag_busy) then
				delta(i) = -100.d0
			endif
		enddo

		call qsortd(delta,ipt,maxdim)
		do i=1,maxdim
			app(maxdim-i+1) = delta(ipt(i))
		!	Lktemp(ndim-i+1,:) = Lkappa(ipt(i),:)
		end do
		!Lkappa = Lktemp
	
		return
	end subroutine spread_ord

end module paretodim
