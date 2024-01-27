subroutine main_MODIRECT(n,nobj, lb0, ub0, nftot, root,first_time_call,maxnf_MODIRECT)
	use mod_type
	use mod_box
	use mod_suddividi
	use mod_mem
	use mod_globale
	use mod_multiobj
	use vincoli
	use alfa_mod
	use paretoset
	
	implicit none

	interface
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
			implicit none


			type(colonna),    pointer		:: root

			integer				:: n, iprint
			real*8				:: tolglob, trigLS
			real*8				:: lbb(n),ubb(n)
		 	real*8				:: xbest(n), fbest, fglob, maxdiam, mindiam, maxL
			real*8				:: xmin(n,100), fmin(100)
			integer				:: imin, jmin
			integer				:: nf, nminloc, maxnf, maxint, nint
			logical				:: trovato,first_time_call
		end subroutine nuovoglobale
	end interface
	type(colonna), pointer	:: root
	integer			:: n, nobj, nftot, maxnf_MODIRECT
	real*8			:: lb0(n), ub0(n)
	logical			:: first_time_call

	character*40		:: nomefun
	character*10		:: fun_name
	character*11		:: tempstr
	character ( len = 64)	:: file_input
	real*8			:: lbb(n), ubb(n), xaux(n), yaux(n), xbest(n)
	real*8			:: fbest, ftemp, mindist, maxdist
	real*8			:: mindiam,maxdiam,maxL, c, fglob
	real			:: rr,time_begin, time_end, timetot, f_max, f_min
	integer			:: imin, jmin
	integer			:: seed(2), iprint, i, nf, numnf, numng, k, ktot, j, ng
	integer			:: maxint, nint, maxnf
	integer*8		:: ninttot, nminloc
	logical			:: trovato
	real*8			:: alfa_stop, tolglob, trigLS
	real*8, allocatable	:: xmin(:,:), fmin(:), fobs(:)
	real, allocatable	:: lbins(:), ubins(:), xr(:)
	external		:: funct
 
	memerror = .false.

	write(*,*) 'MODIRECT: n=',n,' m=',m,' mm=',mm

	!----------------------------------------------------
	! Set number of obj. functions into mod_multiobj
	!----------------------------------------------------
	q = nobj

	!----------------------------------------
	! Allocate storage arrays
	!----------------------------------------
	allocate(vetf1(n),vetf2(n),xsud(n),ysud(n),mask(n))
	allocate(matf1(q,n),matf2(q,n))
	allocate(xmin(n,100),fmin(100))
	allocate(lbins(n),ubins(n),xr(n))
	allocate(globxbest(n),fobs(q))

	!----------------------------------------
	! Set bounds, nomefun and other problem
	! dependent parameters
	!----------------------------------------
	nomefun = 'insean multiobj'
	fglob   = -1d+30
	do i = 1,n
		lbb(i)=lb0(i)
		ubb(i)=ub0(i)
	enddo

	xbest    = (ubb+lbb)/2.d0

	call funct(n,xbest,fbest,fobs)
	write(*,*) 'fbest = ',fbest
	write(*,*) 'xbest = ',xbest
	write(*,*) '  ubb = ',ubb
	write(*,*) '  lbb = ',lbb
    
	!-----------------------------------------------------------------------------
	! Read customizing parameters from file: DIRECT_params.txt
	!-----------------------------------------------------------------------------
	!
	! iprint    (MUST BE an integer) = the printing level
	!
	! rmaxnf    (MUST BE an integer) = relative maximum num. of fun. evaluations
	!           the maximum number of function evaluations is:  
	!              maxnf = rmaxnf*N  where N is the problem dimension
	!
	! maxint    (MUST BE an integer) = maximum number of intervals.
	!           if maxint < 0, then it is set to 1000*maxnf
	!
	! tolglob   (MUST BE a  real)    = tolerance in the stopping condition
	!           code stops when:
	!             abs(fbest-fglob)/max(1.0d0,abs(fglob)) < tolglob
	!
	! alfa_stop (MUST BE a  real)    = required accuracy for the LineSearches
	!           if (alfa_stop < 1.e-6) then alfa_stop is set to 1.e-6 
	!
	! trigLS    (MUST BE a  real)    = percentage of maxnf function evaluations
	!           performed BEFORE first Linesearch is started
	!           trigLS must be between 0.0 and 1.0
	! 
	!-----------------------------------------------------------------------------
	iprint = 0
	maxnf  = maxnf_MODIRECT !500*n
	maxint = -1
	if(maxint < 0) then
		maxint = min(1000000,abs(1000*maxnf))
	endif
	tolglob = -1.d0
	alfa_stop = 1.d-3
	if(alfa_stop < 1.d-6) alfa_stop = 1.d-6
	trigLS = 0.1d0
	if(trigLS < 0.d0) trigLS = 0.d0
	if(trigLS > 1.d0) trigLS = 1.d0
	
	nint      = 0
	nftot     = 1
	ninttot   = 0
	nminloc   = 0
	imin      = 0
	jmin      = 0

	globxbest = xbest
	globfbest = fbest
	globnf    = 1
	globnftot = 1
	globmaxnf = maxnf

	call cpu_time(time_begin)

	call nuovoglobale(root,n,lbb,ubb,xbest,fbest,nf,ng,fglob,maxint,maxnf,iprint, &
			  trovato, nint,mindiam,maxdiam,maxL,xmin,fmin,imin,jmin,tolglob,trigLS,first_time_call)

	call cpu_time(time_end)
	timetot = time_end - time_begin

	nftot   = nftot + nf
	ninttot = ninttot + nint
	nminloc = nminloc + ng

	write(*,*)
	write(*,*) '---------- sommario risultati -----------'
	write(*,*) '  cpu_time = ',timetot, ' secondi'
	write(*,*) '        nf = ',nftot
	write(*,*) '        ng = ',ng
	write(*,*) '     fbest = ',fbest
	write(*,*) '     xbest = ',xbest
	write(*,*) '-----------------------------------------'

	write(1,*)
	write(1,*) '---------- sommario risultati -----------'
	write(1,*) '  cpu_time = ',timetot, ' secondi'
	write(1,*) '        nf = ',nftot
	write(1,*) '        ng = ',ng
	write(1,*) '     fbest = ',fbest
	write(1,*) '     xbest = ',xbest
	write(1,*) '-----------------------------------------'


	!------------------ DEBUG PARETOSET -------------------------
	call print_filter(n,q)
	!------------------ DEBUG PARETOSET -------------------------

	deallocate(vetf1,vetf2,xsud,ysud,mask)
	deallocate(matf1,matf2)
	deallocate(xmin,fmin)
	deallocate(lbins,ubins,xr)
	deallocate(globxbest,fobs)

100 format(a40)

800 FORMAT(a20,' & ',i4,' &  \bf ', es16.8,4(' & ',i11),' & ', es10.2,' & ', es10.2,' & ', es10.2,' \\')
801 FORMAT(a20,' & ',i4,' &      ', es16.8,4(' & ',i11),' & ', es10.2,' & ', es10.2,' & ', es10.2,' \\')
802 FORMAT(a20,' & ',i4,' &   *  ', es16.8,4(' & ',i11),' & ', es10.2,' & ', es10.2,' & ', es10.2,' \\')

900 FORMAT(a20,' & ',i4,' &  \bf ', es16.8,3(' & ',i11),' \\')
901 FORMAT(a20,' & ',i4,' &      ', es16.8,3(' & ',i11),' \\')

end subroutine main_MODIRECT

subroutine funct(n,x,f1,fob)
	use mod_multiobj
	use alfa_mod
	use eps_mod
	implicit none
	integer		:: n
	real*8		:: x(n), f1, fob(q)

	!--------------------------------------------
	! computes the vector of objective functions
	!--------------------------------------------
	call functs(n,x,q,fob)
	if(mm.ge.1) call fconstriq(n,mm,x,ciq_m)
	call functs_pen(n, q, mm, x, finiz_m, fob, ciq_m)	
	fob = finiz_m

	!--------------------------------------------
	! assigns f(1) to f1 and returns
	!--------------------------------------------
	f1 = fob(1)

	return
end subroutine funct
