program main
	implicit none
	integer :: maxnf_MODIR
	integer :: maxnf_DFMO
	integer :: n, m, q

	call setdim(n,m,q)

	maxnf_MODIR = 0 !500*n
	maxnf_DFMO  = 20000 - maxnf_MODIR

	call wrapper(maxnf_MODIR,maxnf_DFMO)

end program main
