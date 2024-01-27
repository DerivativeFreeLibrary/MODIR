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
subroutine riduciconvexhull(convexhull,nelim,eps,toldiam,fmin)
	use mod_type
	implicit none

	type(vertice),pointer		:: convexhull
	type(vertice),pointer		:: currch
	integer						:: nelim, nugua
	real*8						:: eps, toldiam, fmin, L
	logical						:: halt
	
	nelim = 0
	nugua = 0
	halt = .false.
	currch => convexhull
	do while (.not.halt)
		if(associated(currch)) then
			if(currch%int%diam < toldiam) then
				nelim = nelim + 1
				currch => currch%next
				cycle
			endif
		else
			halt = .true.
			exit
		endif

		if(associated(currch%next)) then
			if((currch%next%int%diam - currch%int%diam) > 0.d0) then
			  
				L = (currch%next%int%fint - currch%int%fint) / (currch%next%int%diam - currch%int%diam)

				if( currch%int%fint - L*currch%int%diam >  fmin - eps*max(abs(fmin),1.d-6)  ) then
					
					nelim = nelim + 1 + nugua
					nugua = 0
					currch => currch%next

				else
					halt = .true.
				endif
			else
				nugua = nugua + 1
				currch => currch%next
			endif
		else
			halt = .true.
		endif
	enddo

	nullify(currch)

	return
end subroutine riduciconvexhull

subroutine riduciconvexhull_altre(convexhull,nelim,eps,toldiam,fmin)
	use mod_type
	use mod_multiobj
	implicit none

	type(vertice),pointer		:: convexhull
	type(vertice),pointer		:: currch
	integer						:: nelim, nugua
	real*8						:: eps, toldiam, fmin, L
	logical						:: halt
	
	nelim = 0
	nugua = 0
	halt = .false.
	currch => convexhull
	do while (.not.halt)
		if(associated(currch)) then
			if(currch%int%diam < toldiam) then
				nelim = nelim + 1
				currch => currch%next
				cycle
			endif
		else
			halt = .true.
			exit
		endif

		if(associated(currch%next)) then
			if((currch%next%int%diam - currch%int%diam) > 0.d0) then
			  
				L = (currch%next%int%fobs(i_obj) - currch%int%fobs(i_obj)) / (currch%next%int%diam - currch%int%diam)

				if( currch%int%fobs(i_obj) - L*currch%int%diam >  fmin - eps*max(abs(fmin),1.d-6)  ) then
					
					nelim = nelim + 1 + nugua
					nugua = 0
					currch => currch%next

				else
					halt = .true.
				endif
			else
				nugua = nugua + 1
				currch => currch%next
			endif
		else
			halt = .true.
		endif
	enddo

	nullify(currch)

	return
end subroutine riduciconvexhull_altre

subroutine ricintervallo_dx(root,convexhull,nconv,iprint,Ltilde,eps,fmin)
	use mod_type
	use mod_mem
	implicit none
 
	type(colonna),pointer		:: root, currcol, tempcol, ultimacol
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo
	real*8				:: Ltilde, stimaL, fh, dh, eps, maxL, maxLtemp
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam, fmin
	logical				:: halt
	integer				:: nconv, iprint, sv

	!search the column with max diameter
	ultimacol => root
	do while (associated(ultimacol%next))
		ultimacol => ultimacol%next
	enddo

	!ultimacol points to the column with maximum diameter
	!record the first interval (top, right in the CVX HULL)
	primo => ultimacol%int

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	convexhull%int => primo
	primo%flagcon = .true.
   	currch => convexhull
	currcol => root

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '  id          diam           fint'
		write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fint
		write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*) '     diam           fint',fmin
		write(24,*) convexhull%int%diam, convexhull%int%fint
	endif

	nconv = 1

	maxL = -1.d+10
	currcol => ultimacol

	do while (associated(currcol%pred))
		if(maxL < (ultimacol%int%fint - currcol%pred%int%fint)/(ultimacol%diam - currcol%pred%diam)) then
		  maxL = (ultimacol%int%fint - currcol%pred%int%fint)/(ultimacol%diam - currcol%pred%diam)
		endif
		currcol => currcol%pred
	enddo


	halt = .false.
	if (associated(ultimacol%pred)) then
		currcol => ultimacol%pred
	else
		nullify(currch%next)
		if(iprint > 0) then
			write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			write(24,*)
		endif
		return
	endif
	
	do while (associated(currcol%pred))
		tempcol => currcol%pred
		maxLtemp = -1.d+10
		!do while ((maxLtemp <= maxL).and.(associated(tempcol)))
		do while (associated(tempcol))
			!write(*,*) '=====',maxLtemp,(currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam),maxL,tempcol%int%fint
			if(maxLtemp < (currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam)) then
				maxLtemp = (currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam)
			endif
			tempcol => tempcol%pred
		enddo
		!write(*,*) '[[[[[[[[[[[[',currcol%diam,currcol%int%fint,maxLtemp

		if(maxLtemp <= maxL) then

			!write(*,*) currcol%int%fint,fmin,eps,currcol%diam,maxL
			!write(*,*) (currcol%int%fint - fmin + eps*abs(fmin))/currcol%diam, maxL
			!write(*,*)
			if((currcol%int%fint - fmin + eps*abs(fmin))/currcol%diam <= maxL) then
				allocate(currch%next, stat = sv)
				if(sv.ne.0) then
					memerror = .true.
					return
				endif
				currch%next%int => currcol%int
				currch => currch%next
				if(iprint > 0) then
					write(*,*) currch%int%id, currch%int%diam, currch%int%fint
					write(24,*) currch%int%diam, currch%int%fint
				endif
				nconv  = nconv + 1
				maxL = maxLtemp
			else
				exit
			endif

		endif
		currcol => currcol%pred
		
	enddo

	!write(*,*) root%int%fint,fmin,eps,root%diam,maxL
	!write(*,*) (root%int%fint - fmin + eps*abs(fmin))/root%diam, maxL
	!write(*,*)
	if((root%int%fint - fmin + eps*abs(fmin))/root%diam <= maxL) then
		allocate(currch%next, stat = sv)
		if(sv.ne.0) then
			memerror = .true.
			return
		endif
		currch%next%int => root%int
		currch => currch%next
		if(iprint > 0) then
			write(*,*) currch%int%id, currch%int%diam, currch%int%fint
			write(24,*) currch%int%diam, currch%int%fint
		endif
		nconv  = nconv + 1
	endif

	nullify(currch%next)

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*)
	endif

	return
end subroutine ricintervallo_dx

subroutine ricintervallo(root,convexhull,nconv,iprint,Ltilde)
	use mod_type
	use mod_mem
	use mod_multiobj
	implicit none
 
	type(colonna),pointer		:: root, currcol, aux
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo
	real*8				:: Ltilde, stimaL, fh, dh
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam
	logical				:: halt
	integer				:: nconv, iprint, sv

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	currcol  => root
	aux      => root
	minfunc  = root%int%fint
	maxdiam  = root%diam
	primo    => root%int

	!----------------------------------------
	! search for the interval with max. diam.
	! and min. objective function
	!----------------------------------------
	do while (associated(currcol%next))
		if( ((currcol%next%int%fint < minfunc).or.		&
		    ((currcol%next%diam > maxdiam).and.(currcol%next%int%fint - minfunc <= 1.d-9))).and. &
			(currcol%next%diam > 1.d-8) ) then
			primo   =>  currcol%next%int
			aux     =>  currcol%next	
			maxdiam =  primo%diam
			minfunc =  primo%fint
		endif
		currcol => currcol%next
	enddo

	currcol => aux%next

	!--------------------------------------
	! record the first interval belonging
	! to the convex hull so far identified
	!--------------------------------------
	convexhull%int => primo
	primo%flagcon = .true.
    	currch => convexhull

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '         i_obj = ', i_obj
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '  id          diam           fint'
		write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fint
		write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*) '      diam           fint'
		write(24,*) convexhull%int%diam, convexhull%int%fint
	endif

	nconv = 1

	halt = .false.
	
	do while (.not.halt)
		!-------------------------------------
		! among those intervals in the upper
		! right region, find the one with
		! maximum cosine with vector (1,0)
		!-------------------------------------
		maxcos = -1.d0
		stimaL = 0.d0
		do while (associated(currcol))
			norma = dsqrt((currcol%diam - currch%int%diam)**2.d0+(currcol%int%fint - currch%int%fint)**2.d0)
			coseno = (currcol%diam - currch%int%diam) / norma			
			if(coseno > maxcos) then
				!write(*,*) 'coseno ',coseno,maxcos
				stimaL = (currcol%int%fint - currch%int%fint)/(currcol%diam - currch%int%diam)
				maxcos = coseno
				primo => currcol%int

				aux   => currcol
			endif
			currcol => currcol%next
		enddo
		currcol => aux%next
		if(stimaL > Ltilde) exit
		if(maxcos > 0.d0) then
			allocate(currch%next, stat = sv)
			if(sv.ne.0) then
				memerror = .true.
				return
			endif
			currch%next%int => primo
			currch => currch%next
			nconv  = nconv + 1
			primo%flagcon = .true.
			if (iprint > 0)	then
				write(*,*) currch%int%id, currch%int%diam, currch%int%fint
				write(24,*) currch%int%diam, currch%int%fint
			endif
		else
			halt = .true.
		endif
	enddo

	nullify(currch%next)
	
	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*)
	endif

	return
end subroutine ricintervallo

subroutine ricintervallo_altre(root,convexhull,nconv,iprint,Ltilde)
	use mod_type
	use mod_mem
	use mod_multiobj
	implicit none
 
	interface
		subroutine search_min(col,i_obj,punt,minfunc)
			use mod_type
			implicit none
			type(colonna), pointer		:: col
			type(intervallo), pointer	:: punt
			integer				:: i_obj 
			real*8				:: minfunc
		end subroutine search_min
	end interface

	type(colonna),pointer		:: root, currcol, aux
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo, temp
	real*8				:: Ltilde, stimaL, fh, dh, tempmin
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam
	logical				:: halt
	integer				:: nconv, iprint, sv

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	currcol  => root
	aux      => root

	maxdiam  = root%diam
	call search_min(root,i_obj,primo,minfunc)

	!minfunc  = root%int%fint
	!primo    => root%int

	!----------------------------------------
	! search for the interval with max. diam.
	! and min. objective function
	!----------------------------------------
	do while (associated(currcol%next))

		call search_min(currcol%next,i_obj,temp,tempmin)

		if( ((tempmin < minfunc).or.		&
		    ((currcol%next%diam > maxdiam).and.(tempmin - minfunc <= 1.d-9))).and. &
			(currcol%next%diam > 1.d-8) ) then
			primo   =>  temp
			aux     =>  currcol%next	
			maxdiam =  primo%diam
			minfunc =  tempmin
		endif
		currcol => currcol%next
	enddo

	currcol => aux%next

	!--------------------------------------
	! record the first interval belonging
	! to the convex hull so far identified
	!--------------------------------------
	convexhull%int => primo
	primo%flagcon = .true.
    	currch => convexhull

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '         i_obj = ', i_obj
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '  id          diam           fint'
		write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fobs(i_obj)
		write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*) '      diam           fint'
		write(24,*) convexhull%int%diam, convexhull%int%fobs(i_obj)
	endif

	nconv = 1

	halt = .false.
	
	do while (.not.halt)
		!-------------------------------------
		! among those intervals in the upper
		! right region, find the one with
		! maximum cosine with vector (1,0)
		!-------------------------------------
		maxcos = -1.d0
		stimaL = 0.d0
		do while (associated(currcol))
			call search_min(currcol,i_obj,temp,tempmin)
			norma = dsqrt((currcol%diam - currch%int%diam)**2.d0+(tempmin - currch%int%fobs(i_obj))**2.d0)
			coseno = (currcol%diam - currch%int%diam) / norma			
			if(coseno > maxcos) then
				!write(*,*) 'coseno ',coseno,maxcos
				stimaL = (tempmin - currch%int%fobs(i_obj))/(currcol%diam - currch%int%diam)
				maxcos = coseno
				primo => temp
				aux   => currcol
			endif
			currcol => currcol%next
		enddo
		currcol => aux%next
		if(stimaL > Ltilde) exit
		if(maxcos > 0.d0) then
			allocate(currch%next, stat = sv)
			if(sv.ne.0) then
				memerror = .true.
				return
			endif
			currch%next%int => primo
			currch => currch%next
			nconv  = nconv + 1
			primo%flagcon = .true.
			if (iprint > 0)	then
				write(*,*) currch%int%id, currch%int%diam, currch%int%fobs(i_obj)
				write(24,*) currch%int%diam, currch%int%fobs(i_obj)
			endif
		else
			halt = .true.
		endif
	enddo

	nullify(currch%next)
	
	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*)
	endif

	return
end subroutine ricintervallo_altre

subroutine ricintervallo_maxnnondom(root,convexhull,nconv,iprint,Ltilde)
	use mod_type
	use mod_mem
	use mod_multiobj
	implicit none
 
	type(colonna),pointer		:: root, currcol, aux
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo, temp
	real*8				:: Ltilde, stimaL, fh, dh, tempmin
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam
	logical				:: halt
	integer				:: nconv, iprint, sv, maxnnondom
	real*8, parameter		:: maxfnondom = 1.d+4

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	currcol  => root
	convexhull%int => root%int
	maxnnondom = root%int%nnondom

	do while (associated(currcol))
		temp => currcol%int
		do while (associated(temp))
			if(temp%nnondom > maxnnondom) then
				maxnnondom = temp%nnondom
				convexhull%int => temp
			endif
			temp => temp%next
		enddo
		currcol => currcol%next
	enddo
	nullify(convexhull%next)

	return
end subroutine ricintervallo_maxnnondom

subroutine ricintervallo_nnondom(root,convexhull,nconv,iprint,Ltilde)
	use mod_type
	use mod_mem
	use mod_multiobj
	implicit none
 
	interface
		subroutine search_min_nnondom(col,punt,minfunc)
			use mod_type
			implicit none
			type(colonna), pointer		:: col
			type(intervallo), pointer	:: punt
			real*8				:: minfunc
		end subroutine search_min_nnondom
	end interface

	type(colonna),pointer		:: root, currcol, aux
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo, temp
	real*8				:: Ltilde, stimaL, fh, dh, tempmin
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam
	logical				:: halt
	integer				:: nconv, iprint, sv
	real*8, parameter		:: maxfnondom = 1.d+4

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	currcol  => root
	aux      => root

	maxdiam  = root%diam
	call search_min_nnondom(root,primo,minfunc)

	!minfunc  = root%int%fint
	!primo    => root%int

	!----------------------------------------
	! search for the interval with max. diam.
	! and min. objective function
	!----------------------------------------
	do while (associated(currcol%next))

		call search_min_nnondom(currcol%next,temp,tempmin)

		if( ((tempmin < minfunc).or.		&
		    ((currcol%next%diam > maxdiam).and.(tempmin - minfunc <= 1.d-9))).and. &
			(currcol%next%diam > 1.d-8) ) then
			primo   =>  temp
			aux     =>  currcol%next	
			maxdiam =  primo%diam
			minfunc =  tempmin
		endif
		currcol => currcol%next
	enddo

	currcol => aux%next

	!--------------------------------------
	! record the first interval belonging
	! to the convex hull so far identified
	!--------------------------------------
	convexhull%int => primo
	primo%flagcon = .true.
    	currch => convexhull

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '         ricerca su 1/nnondom  '
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '  id          diam           fint'
		write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*)'  id          diam           fint'
		if(convexhull%int%nnondom > 0) then
			write(*,*) convexhull%int%id, convexhull%int%diam, 1.d0/convexhull%int%nnondom
			write(24,*)convexhull%int%id, convexhull%int%diam, 1.d0/convexhull%int%nnondom
		else
			write(*,*) convexhull%int%id, convexhull%int%diam, maxfnondom
			write(24,*)convexhull%int%id, convexhull%int%diam, maxfnondom
		endif
	endif

	nconv = 1

	halt = .false.
	
	do while (.not.halt)
		!-------------------------------------
		! among those intervals in the upper
		! right region, find the one with
		! maximum cosine with vector (1,0)
		!-------------------------------------
		maxcos = -1.d0
		stimaL = 0.d0
		do while (associated(currcol))
			call search_min_nnondom(currcol,temp,tempmin)
			if(currch%int%nnondom > 0) then
				norma = dsqrt((currcol%diam - currch%int%diam)**2.d0+(tempmin - 1.d0/currch%int%nnondom)**2.d0)
			else
				norma = dsqrt((currcol%diam - currch%int%diam)**2.d0+(tempmin - maxfnondom)**2.d0)
			endif
			coseno = (currcol%diam - currch%int%diam) / norma			
			if(coseno > maxcos) then
				!write(*,*) 'coseno ',coseno,maxcos
				if(currch%int%nnondom > 0) then
					stimaL = abs(tempmin - 1.d0/currch%int%nnondom)/abs(currcol%diam - currch%int%diam)
				else
					stimaL = abs(tempmin - maxfnondom)/abs(currcol%diam - currch%int%diam)
				endif
				maxcos = coseno
				primo => temp
				aux   => currcol
			endif
			currcol => currcol%next
		enddo
		currcol => aux%next
		if(stimaL > Ltilde) exit
		if(maxcos > 0.d0) then
			allocate(currch%next, stat = sv)
			if(sv.ne.0) then
				memerror = .true.
				return
			endif
			currch%next%int => primo
			currch => currch%next
			nconv  = nconv + 1
			primo%flagcon = .true.
			if (iprint > 0)	then
				if(currch%int%nnondom > 0) then 
					write(*,*) currch%int%id, currch%int%diam, 1.d0/currch%int%nnondom
					write(24,*)currch%int%id, currch%int%diam, 1.d0/currch%int%nnondom
				else
					write(*,*) currch%int%id, currch%int%diam, maxfnondom
					write(24,*)currch%int%id, currch%int%diam, maxfnondom
				endif
			endif
		else
			halt = .true.
		endif
	enddo

	nullify(currch%next)
	
	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(24,*)
	endif

	return
end subroutine ricintervallo_nnondom

subroutine suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
	use mod_type
	use mod_globale
	use mod_suddividi
	use mod_mem
	use mod_box
	use mod_multiobj
	implicit none

	type(vertice), pointer		:: currch
	type(intervallo), pointer	:: curr
	type(colonna),    pointer	:: root, temp
	integer						:: n, nf, nint, cont
	real*8						:: xdir(n), fdir(q), xdir_unscaled(n)

	integer						:: i
	integer						:: numtrue, ind1(1), ind2(1)
	real*8						:: maxL,Ltilde,tol
	real*8						:: fobs(q)
	logical						:: flder

	interface
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8						:: tol
		end subroutine find_colonna
		subroutine insert_intervallo(curcol,int)
			use mod_type
			implicit none

			type(colonna),    pointer	:: curcol
			type(intervallo), pointer	:: int
		end subroutine insert_intervallo
		subroutine triplica(primo,root,n,ind,f1,f2,fobs1,fobs2,nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flag,cont,tol)
			use mod_type
			use mod_box
			use mod_mem
			use mod_multiobj
			implicit none

			type(intervallo), pointer	:: primo
			type(colonna),    pointer	:: root
			integer						:: n, ind, nint, cont
			real*8						:: xdir(n), fdir(q),xdir_unscaled(n)
			real*8						:: f1, f2, maxL, Ltilde, tol
			real*8						:: fobs1(q), fobs2(q)
			logical						:: flag
		end subroutine triplica
	end interface

	curr => currch%int
	numtrue = 0

	do i = 1,n
		if(curr%maxdim == curr%dimen(i)) then
			ysud = curr%cent 
			ysud(i) = curr%cent(i) + 1.d0*curr%dimen(i)/3.d0
			call unscalevars(n,ysud,ytemp,xbar,lbs,ubs)
			call unscalevars_direct(n,ytemp,xsud)
			call funct(n,xsud,vetf1(i),matf1(:,i))
		      if((vetf1(i) < globfbest).and.(nf <= globmaxnf)) then
		 	globxbest = xsud
			globfbest = vetf1(i)
			globnftot = nf
			globnf    = nf
		      endif

			ysud(i) = curr%cent(i) - 1.d0*curr%dimen(i)/3.d0
			call unscalevars(n,ysud,ytemp,xbar,lbs,ubs)
			call unscalevars_direct(n,ytemp,xsud)
			call funct(n,xsud,vetf2(i),matf2(:,i))
		      if((vetf2(i) < globfbest).and.(nf <= globmaxnf)) then
		 	globxbest = xsud
			globfbest = vetf2(i)
			globnftot = nf
			globnf    = nf
		      endif
			mask(i) = .true.
			numtrue = numtrue + 1

			nf = nf+2
		else
			vetf1(i) = 1.d+30
			vetf2(i) = 1.d+30
			matf1(:,i) = 1.d+30
			matf2(:,i) = 1.d+30
			mask(i)  = .false.
		endif
	enddo

	curr%der = 0.d0
	flder    = .false.

	do i = 1,numtrue
		ind1 = minloc(matf1(i_obj,:),mask)
		ind2 = minloc(matf2(i_obj,:),mask)
		if(matf1(i_obj,ind1(1)) < matf2(i_obj,ind2(1))) then
			mask(ind1(1)) = .false.
			!write(*,*) '======= suddividi: chiamo triplica', numtrue
			call triplica(curr,root,n,ind1(1),vetf1(ind1(1)),vetf2(ind1(1)),matf1(:,ind1(1)), matf2(:,ind1(1)), &
				      nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flder,cont,tol)
			if (memerror) then
				write(*,*)'fine memoria disponibile mentre triplico'
				return
			endif
		else
			mask(ind2(1)) = .false.
			!write(*,*) '======= suddividi: chiamo triplica', numtrue
			call triplica(curr,root,n,ind2(1),vetf1(ind2(1)),vetf2(ind2(1)),matf1(:,ind2(1)), matf2(:,ind2(1)), &
				      nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flder,cont,tol)
			if (memerror) then
				write(*,*)'fine memoria disponibile mentre triplico'
				return
			endif
		endif
		nint = nint + 2
	enddo

	!write(*,*) 'controllo il centrale',associated(root),curr%id,curr%diam
	if(curr%fint - Ltilde*curr%diam/2.d0 <= fdir(1)) then
	!if(.true.) then
		curr%flagcon = .false.
		!write(*,*) '-------insert   primo--------',associated(root),curr%id,curr%diam
		!if(associated(root)) write(*,*) '-------insert   primo--------',associated(root),associated(root%pred),associated(root%next),root%diam
		call find_colonna(root,curr,tol,temp)
		call insert_intervallo(temp,curr)
		!nullify(curr)
	else
		!write(*,*) '-------elimino  primo--------',curr%id

		deallocate(curr%cent,curr%dimen)
		nullify(curr%next)
		nullify(curr%pred)
		deallocate(curr)
		!deallocate(currch%int)
		!nullify(currch%int)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif

	return
end subroutine suddividi

subroutine triplica(primo,root,n,ind,f1,f2,fobs1,fobs2,nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flag,cont,tol)
	use mod_type
	use mod_box
	use mod_mem
	use mod_multiobj
	use paretoset
	use paretodim
	implicit none

	interface
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8						:: tol
		end subroutine find_colonna
		subroutine insert_intervallo(curcol,int)
			use mod_type
			implicit none

			type(colonna),    pointer	:: curcol
			type(intervallo), pointer	:: int
		end subroutine insert_intervallo
	end interface

	type(intervallo), pointer	:: primo
	type(colonna),    pointer	:: root, temp
	type(typ_lista), pointer	:: templist, primolist

	integer						:: n, ind, nint, cont, sv, i, j
	real*8						:: xdir(n), fdir(q), xdir_unscaled(n)
	real*8						:: f1, f2, norma, deltax, maxL, Ltilde, tol
	real*8						:: fobs1(q), fobs2(q), fdim(q+1)
	real*8						:: g(n), lb1(n), ub1(n), lb2(n), ub2(n), lb3(n), ub3(n)
	real*8						:: newder
	logical						:: flag, soddisfa, e_il_primo

	type(intervallo), pointer	:: secondo
	type(intervallo), pointer	:: terzo

	allocate(secondo, stat = sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif
	allocate(terzo, stat = sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	call alloca_intervallo(secondo,n,q)
	if (memerror) then
		write(*,*)'fine memoria disponibile mentre triplico'
		return
	endif
	
	call alloca_intervallo(terzo,n,q)
	if (memerror) then
		write(*,*)'fine memoria disponibile mentre triplico'
		return
	endif
 
	!allocate(secondo%cent(n),secondo%dimen(n))
	!allocate(terzo%cent(n),  terzo%dimen(n))

	secondo%cent = primo%cent
	terzo%cent   = primo%cent

	secondo%cent(ind) = secondo%cent(ind) + 1.d0*primo%dimen(ind)/3.d0
	terzo%cent(ind)   = terzo%cent(ind)   - 1.d0*primo%dimen(ind)/3.d0

	secondo%dimen = primo%dimen
	terzo%dimen   = primo%dimen

	primo%dimen(ind)   = primo%dimen(ind)/3.d0
	secondo%dimen(ind) = secondo%dimen(ind)/3.d0
	terzo%dimen(ind)   = terzo%dimen(ind)/3.d0

	primo%maxdim = maxval(primo%dimen)
	primo%diam   = norma(n,primo%dimen)/2.d0

	call aggiorna_dim(n,q,primo)

	secondo%maxdim = maxval(secondo%dimen)
	secondo%diam   = norma(n,secondo%dimen)/2.d0
	secondo%xbars  = primo%xbars
	secondo%lbs    = primo%lbs
	secondo%ubs    = primo%ubs

	terzo%maxdim = maxval(terzo%dimen)
	terzo%diam   = norma(n,terzo%dimen)/2.d0
	terzo%xbars  = primo%xbars
	terzo%lbs    = primo%lbs
	terzo%ubs    = primo%ubs

	secondo%fint = f1
	secondo%fobs = fobs1

	if(f1 < fdir(1)) then
		fdir(1) = f1
		xdir = secondo%cent
		call unscalevars(n,xdir,ytemp,xbar,lbs,ubs)
		call unscalevars_direct(n,ytemp,xdir_unscaled)
		cont = 0
	endif
	xtemp = secondo%cent
	call unscalevars(n,xtemp,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	call update_point(n,q,xtemp,fobs1,secondo,flag)
	fdim(1:q) = fobs1
	fdim(q+1) = -secondo%diam
	call update_point_dim(n,q,xtemp,fdim,secondo,flag)


	terzo%fint = f2
	terzo%fobs = fobs2

	if(f2 < fdir(1)) then
		fdir(1) = f2
		xdir = terzo%cent
		call unscalevars(n,xdir,ytemp,xbar,lbs,ubs)
		call unscalevars_direct(n,ytemp,xdir_unscaled)
		cont = 0
	endif
	xtemp = terzo%cent
	call unscalevars(n,xtemp,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	call update_point(n,q,xtemp,fobs2,terzo,flag)
	fdim(1:q) = fobs2
	fdim(q+1) = -terzo%diam
	call update_point_dim(n,q,xtemp,fdim,secondo,flag)

	do i = 1,q
		if(fobs1(i) < fdir(i)) then
			fdir(i) = fobs1(i)
		endif
		if(fobs2(i) < fdir(i)) then
			fdir(i) = fobs2(i)
		endif
	enddo

	secondo%flagopt = .false.
	terzo%flagopt   = .false.
	primo%flagopt   = .false.

	secondo%flagloc = .false.
	terzo%flagloc   = .false.

	secondo%flagdiv = .true.
	terzo%flagdiv   = .true.

	secondo%flagcon = .false.
	terzo%flagcon   = .false.

	secondo%id      = nint+1
	secondo%nnondom = 0
	nullify(secondo%list_nondom)
	terzo%id        = nint+2
	terzo%nnondom   = 0
	nullify(terzo%list_nondom)

	if(primo%nnondom > 0) then
		!-----------------------------------------------
		! divide la lista dei nondominati del primo 
		! tra primo secondo e terzo
		!-----------------------------------------------
		call unscalevars(n,primo%cent - primo%dimen/2.d0,xtemp, primo%xbars,primo%lbs,primo%ubs)
		call unscalevars_direct(n,xtemp,lb1)
		call unscalevars(n,primo%cent + primo%dimen/2.d0,xtemp, primo%xbars,primo%lbs,primo%ubs)
		call unscalevars_direct(n,xtemp,ub1)

		call unscalevars(n,secondo%cent - secondo%dimen/2.d0,xtemp, secondo%xbars,secondo%lbs,secondo%ubs)
		call unscalevars_direct(n,xtemp,lb2)
		call unscalevars(n,secondo%cent + secondo%dimen/2.d0,xtemp, secondo%xbars,secondo%lbs,secondo%ubs)
		call unscalevars_direct(n,xtemp,ub2)

		!call unscalevars(n,terzo%cent - terzo%dimen/2.d0,xtemp, terzo%xbars,terzo%lbs,terzo%ubs)
		!call unscalevars_direct(n,xtemp,lb3)
		!call unscalevars(n,terzo%cent + terzo%dimen/2.d0,xtemp, terzo%xbars,terzo%lbs,terzo%ubs)
		!call unscalevars_direct(n,xtemp,ub3)

		write(*,*) 'triplica: primo%nnondom=',primo%nnondom
		primolist => primo%list_nondom
		i = 0
		do while (associated(primolist))
			i = i+1
			primolist => primolist%next
		enddo		
		write(*,*) 'triplica: check nnondom di primo e'' ',i
		primo%nnondom = 0
		i = 0
		nullify(primolist)
		primolist => primo%list_nondom

		if(associated(primolist)) then
			write(*,*) 'la lista di nondom di primo non e'' vuota'
			do while (associated(primolist))
				i = i+1
				write(*,*) 'triplica: check 1',primo%nnondom
				write(*,*) 'triplica: check 2',secondo%nnondom
				write(*,*) 'triplica: check 3',terzo%nnondom
				soddisfa = .true.
				do j = 1,n
					if( (primolist%x(j) > ub1(j)).or.(primolist%x(j) < lb1(j)) ) then
						soddisfa = .false.
						exit
					endif
				enddo
				if(soddisfa) then
				!point belongs to primo
					primo%nnondom = primo%nnondom+1
					primolist => primolist%next
				else
					!per prima cosa sgancia l'intervallo dalla lista di primo
					templist => primolist
					primolist => primolist%next
					if(associated(templist%prev) .and. associated(templist%next)) then
						! sta in mezzo
						templist%prev%next => templist%next
						templist%next%prev => templist%prev
					elseif(associated(templist%prev)) then
						! e' l'ultimo
						nullify(templist%prev%next)
					elseif(associated(templist%next)) then
						! e' il primo
						nullify(templist%next%prev)
						primo%list_nondom => templist%next
					else
						! e' l'unico
						nullify(primo%list_nondom)
					endif
					nullify(templist%prev)
					nullify(templist%next)
					!qui ho sganciato l'intervallo dalla lista di primo
					soddisfa = .true.
					do j = 1,n
						if( (templist%x(j) > ub2(j)).or.(templist%x(j) < lb2(j)) ) then
							soddisfa = .false.
							exit
						endif
					enddo
					if(soddisfa) then
					! point belongs to secondo
						if(associated(secondo%list_nondom)) then
							templist%next => secondo%list_nondom
							secondo%list_nondom => templist
							templist%next%prev => templist
						else
							secondo%list_nondom => templist
						endif
						secondo%nnondom = secondo%nnondom+1
					else
					! point belongs to terzo
						if(associated(terzo%list_nondom)) then
							templist%next => terzo%list_nondom
							terzo%list_nondom => templist
							templist%next%prev => templist
						else
							terzo%list_nondom => templist
						endif
						terzo%nnondom = terzo%nnondom+1
					endif
				endif
			enddo
		else
			write(*,*) 'la lista di nondom di primo si e'' svuotata'
		endif
		if(primo%nnondom == 0) then
			nullify(primo%list_nondom)
		endif

		write(*,*) 'triplica:   primo%nnondom=',primo%nnondom,i
		write(*,*) 'triplica: secondo%nnondom=',secondo%nnondom
		write(*,*) 'triplica:   terzo%nnondom=',terzo%nnondom
	endif

	call unscalevars(n,secondo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = xtemp(ind)

	call unscalevars(n,primo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = deltax - xtemp(ind)
	secondo%der		= abs(f1 - primo%fint)/abs(deltax)
!	call grad(xtemp,n,g)
!	secondo%der = norma(n,g)

	if (maxL < secondo%der	) then
		maxL = secondo%der
		cont = 0
	endif

	call unscalevars(n,terzo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = xtemp(ind)

	call unscalevars(n,primo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = deltax - xtemp(ind)
	terzo%der		= abs(primo%fint - f2)/abs(deltax)

	if (maxL < terzo%der ) then
		maxL = terzo%der
		cont = 0
	endif

	primo%der = max(primo%der,abs(secondo%der),abs(terzo%der))

	if(secondo%fint - Ltilde*secondo%diam/2.d0 <= fdir(1)) then
		!write(*,*) '-------insert  secondo--------',associated(root)
		call find_colonna(root,secondo,tol,temp)
		call insert_intervallo(temp,secondo)
		!write(*,*) '-------insert  secondo--------',associated(root),associated(root%pred),associated(root%next)
	else
		!write(*,*) '-------elimino secondo--------'
		deallocate(secondo%cent,secondo%dimen)
		deallocate(secondo)
		nullify(secondo)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif

	if(terzo%fint - Ltilde*terzo%diam/2.d0 <= fdir(1)) then
		!write(*,*) '-------insert   terzo--------',associated(root)
		call find_colonna(root,terzo,tol,temp)
		call insert_intervallo(temp,terzo)
		!write(*,*) '-------insert   terzo--------',associated(root),associated(root%pred),associated(root%next)
	else
		!write(*,*) '-------elimino  terzo--------'
		deallocate(terzo%cent,terzo%dimen)
		deallocate(terzo)
		nullify(terzo)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif
		


	return
end subroutine triplica

subroutine scalevars(n,x,y,xbar,lbs,ubs)
	implicit none

	!-----------------------------------------------------------------
	! dato [lbs , ubs] subset [0 , 1] e xbar in [lbs , ubs]
	! dato x in [lbs , ubs]
	! restituisce y trasformato ricentrando xbar in (lbs+ubs)/2
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n), xbar(n), lbs(n), ubs(n)
	real*8						:: cent(n)
	integer						:: i
	
	cent = (lbs+ubs) / 2.d0

	do i = 1,n
		if(x(i) <= xbar(i)) then
			y(i) = ( (x(i)-lbs(i)) / (xbar(i) - lbs(i)) ) * ( cent(i) -lbs(i) ) + lbs(i)
		else
			y(i) = ( (x(i) - xbar(i)) / (ubs(i) - xbar(i)) ) * ( ubs(i) - cent(i) ) + cent(i)
		endif
	enddo

	return

end subroutine scalevars

subroutine unscalevars(n,y,x,xbar,lbs,ubs)
	implicit none

	!-----------------------------------------------------------------
	! dato [lbs , ubs] subset [0 , 1] e xbar in [lbs , ubs]
	! dato y in [lbs , ubs] trasformato
	! restituisce x in modo che (lbs+ubs)/2 va in xbar
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n), xbar(n), lbs(n), ubs(n)
	real*8						:: cent(n)
	integer						:: i
	
	cent = (lbs+ubs) / 2.d0

	do i = 1,n
		if(y(i) <= cent(i)) then
			x(i) = ( (y(i)-lbs(i)) / (cent(i) - lbs(i)) ) * ( xbar(i) -lbs(i) ) + lbs(i)
		else
			x(i) = ( (y(i) - cent(i)) / (ubs(i) - cent(i)) ) * ( ubs(i) - xbar(i) ) + xbar(i)
		endif
	enddo

	return

end subroutine unscalevars

subroutine scalevars_direct(n,x,y)
	use mod_box
	implicit none

	!-----------------------------------------------------------------
	! dato [lb , ub] e x in [lb , ub]
	! restituisce y in [0 , 1]
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n)
	integer						:: i

	y = (x-lb)/(ub-lb)

	return

end subroutine scalevars_direct

subroutine unscalevars_direct(n,y,x)
	use mod_box
	implicit none

	!-----------------------------------------------------------------
	! dato [lb , ub] e y in [0 , 1]
	! restituisce x in [lb , ub]
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n)
	integer						:: i

	x = lb + y*(ub-lb)
	
	return
end subroutine unscalevars_direct

real*8 function norma(n,x)
	implicit none
	
	integer						:: n
	real*8						:: x(n)
	
	norma = dsqrt(dot_product(x,x))

	return 
end function norma

