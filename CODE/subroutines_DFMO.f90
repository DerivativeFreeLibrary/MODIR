subroutine insert_row(n,q,ndim,row,Lout,flag_change)
	implicit none
	integer n,q,ndim
	real*8 row(n+q+2+n), Lout(ndim,n+q+2+n)
	integer i, j, indfree
	integer flag
	logical flag_change, domina, minore_uguale
	!--------------------------------------------------------------
	! flag assume tre valori: 1,2,3
	! flag = 1 : Lout(i,n+1:n+q) < row(n+1:n+q) 
	! flag = 2 : row(n+1:n+q) < Lkappa(i,n+1:n+q)
	! flag = 3 : row(n+1:n+q) non domina e non e' dominato da Lkappa(i,n+1:n+q)
	!--------------------------------------------------------------
	flag_change = .false.
        flag=3
	indfree = -1
	do i = 1,ndim
		if((Lout(i,n+q+1) < 50.d0).and.(indfree == -1)) then
			indfree = i
		elseif(Lout(i,n+q+1) > 50.d0) then
			if( (minore_uguale(q,Lout(i,n+1:n+q),row(n+1:n+q))) .and. &
		            (minore_uguale(q,row(n+1:n+q),Lout(i,n+1:n+q))) ) then
				flag = 1
				exit
			endif
			if(domina(q,Lout(i,n+1:n+q),row(n+1:n+q))) then
				flag = 1
			elseif(domina(q,row(n+1:n+q),Lout(i,n+1:n+q))) then
				flag = 2
			else
				flag = 3
			endif
			!flag = 0
			!do j = 1,q
			!	if(Lkappa(i,n+j) >= f(j)) then
			!		flag = .false.
			!		exit
			!	endif
			!enddo
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
				Lout(i,n+q+1) = 0.d0
				if(indfree == -1) indfree = i
			case(3)
			end select
		endif
	enddo
	!write(*,*) 'insert row: flag = ',flag,indfree
	if(flag /= 1) then
		flag_change = .true.
		if(indfree /= -1) then
			Lout(indfree,1:n) = row(1:n)
			Lout(indfree,n+1:n+q) = row(n+1:n+q)
			Lout(indfree,n+q+1) = 100.d0
			Lout(indfree,n+q+2) = row(n+q+2)
			Lout(indfree,n+q+2+1:n+q+2+n) = row(n+q+2+1:n+q+2+n)
			!write(*,*) 'insert row: ', indfree, Lout(indfree,n+q+2+1:n+q+2+n)
		else
			write(*,*) 'WARNING: insert_row'
			!pause
		endif
	endif
	return
end subroutine insert_row

subroutine merge(n,q,ndim,Lin,Lout)
	implicit none
	integer, intent( IN ) 	:: n, q, ndim
	real*8, intent( IN )	:: Lin(ndim,n+q+2+n)
	real*8, intent( INOUT )	:: Lout(ndim,n+q+2+n)
	integer			:: i
	logical 		:: flag_change, flag
	integer, parameter	:: iprint = 0

	flag_change = .false.
	if(iprint > 0) then
		write(*,*) '---------- MERGE START ----------',n,q
		do i = 1,ndim
			if(Lin(i,n+q+1) > 50.d0) then
				write(*,*) 'Lin:',i,Lin(i,1:n+q)
			endif
			if(Lout(i,n+q+1) > 50.d0) then
				write(*,*) 'Lou:',i,Lout(i,n+1:n+q)
			endif
		enddo
	endif
	do i = 1,ndim
		if(Lin(i,n+q+1) < 50.d0) cycle
		call insert_row(n,q,ndim,Lin(i,:),Lout,flag)

		flag_change = flag_change.or.flag
	enddo
	if(iprint > 0) then
		do i = 1,ndim
			if(Lout(i,n+q+1) > 50.d0) then
				write(*,*) 'Lou:',i,Lout(i,1:n+q)
			endif
		enddo
		write(*,*) '---------- MERGE END ----------'
	endif
	!pause
	return
end subroutine merge

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

logical function Lkappa_domina(n,q,f,rho)
	use filtro
	implicit none
	integer		:: n,q
	real*8		:: f(q),rho
	integer		:: i
	logical 	:: domina

	Lkappa_domina = .false.

	do i = 1,ndim
		if(Lkappa(i,n+q+1) > 50.d0) then
			if(domina(q,Lkappa(i,n+1:n+q)-rho,f)) then
				Lkappa_domina = .true.
				return
			endif
		endif
	enddo
	return
end function Lkappa_domina

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
	if(.not.domina) return
	domina = .false.
	do i = 1,q
		if(f1(i) < f2(i)) then
			domina = .true.
			exit
		endif
	enddo

	return
end function domina


logical function maggiore_stretto(q,f1,f2)
	implicit none
	integer q,i
	real*8 f1(q),f2(q)
	!---------------------------------------------------
	! return TRUE   if f1 maggiore stretto f2
	!        FALSE  otherwise
	!---------------------------------------------------
	maggiore_stretto = .true.
	do i = 1,q
		if(f1(i) <= f2(i)) then
			maggiore_stretto = .false.
			exit
		endif
	enddo
	
	return
end function maggiore_stretto
