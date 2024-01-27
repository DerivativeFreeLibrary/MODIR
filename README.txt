-----------------------------------------------------------
 How to run the derivative-free optimizer MODIR (and, 
 optionally, the hybrid alg. MODR+DFMO) for inequality 
 constrained multiobjective global optimization problems
-----------------------------------------------------------
 The package provides a Fortran90 version of the MODIR algorithm.

0- Gunzip and untar the archive in a folder on your computer by
   issuing in a directory of your choice (ex. curdir) the command

   $> tar -zxvf MODIR.tar.gz

1- Edit the file main.f90:

	maxnf_MODIR		maximum number of function evaluations for MODIR
	
	maxnf_DFMO		maximum number of function evaluations for DFMO

2- Edit the file problem.f90 providing, in particular:

	subroutine setdim(n,m,q) : problem dimensions
		n = number of variables
		m = number of inequality constraints (g(x) <= 0)
		q = number of objective functions to be (simoultaneously) minimized

	subroutine startp(n,x) : initial point

	subroutine functs(n,x,q,f) 
		subroutine that computes the objective function values (f) on a given point (x)

	subroutine setbounds(n,l,u)
		subroutine that defines lower and upper bounds on the variables

	subroutine fconstriq(n,m,x,ciq)
		subroutine that computes contraint values on a given point (x)

3- At command prompt in curdir execute 

     $> make 
 
   which will create the executable 'modir_dfmo'

4- execute

     $> ./modir_dfmo

   and you are done!
