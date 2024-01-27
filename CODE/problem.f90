!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#
!#   Example BK1, see the previous cited paper for the original reference.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#   Written by the authors in June 1, 2010.
!#
!###############################################################################
!
!# x variable
!var x{1..2} >=-5, <=10;
!
!minimize f1:
!	x[1]^2+x[2]^2;
!minimize f2:
!    (x[1]-5)^2+(x[2]-5)^2;

subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 2
	m = 1
	q = 2

	return
end subroutine setdim

subroutine startp(n,x)
	implicit none
	integer	:: n
	real*8		:: x(n)

	x = 0.d0
	
	return
end subroutine startp

subroutine functs(n,x,q,f)
	implicit none
	integer	:: n,q
	real*8		:: x(n), f(q)

	f(1) = x(1)**2.d0 + x(2)**2.d0
	f(2) = (x(1)-5.d0)**2.d0 + (x(2)-5.d0)**2.d0

	return
end subroutine functs

subroutine setbounds(n,l,u)
	implicit none
	integer	:: n
	real*8		:: l(n), u(n)

	l = -5.d0
	u = 10.d0

	return
end subroutine setbounds

subroutine fconstriq(n,m,x,ciq)
	implicit none
	integer n,m,I
	real*8 x(n), ciq(m)

	do I = 1,m
		ciq(I) = x(I)**2 + x(I+1)**2 + &
			   x(I)*x(I+1)-2.D0*x(I)-2.D0*x(I+1)+1.D0
	enddo

	return
end subroutine fconstriq

