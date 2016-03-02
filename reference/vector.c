static const char help[] = "A help string, printed on -help\n\n";

/*
   Include "petscdmda.h" so that we can use distributed arrays (DMDAs).
   Include "petscsnes.h" so that we can use SNES solvers.  Note that this
   file automatically includes:
     petsc.h       - base PETSc routines   petscvec.h - vectors
     petscsys.h    - system routines       petscmat.h - matrices
     petscis.h     - index sets            petscksp.h - Krylov subspace methods
     petscviewer.h - viewers               petscpc.h  - preconditioners
     petscksp.h    - linear solvers        petscdm.h  - generic grid/physics management
*/
#include "petscdmda.h"
#include "petscsnes.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Vec                x,y,w;               /* vectors */
  PetscErrorCode     ierr;
  PetscMPIInt        rank;
  PetscScalar        norm;
  PetscScalar        one = 1.0,two = 2.0;
  PetscInt           i,N;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInitialize(&argc,&argv,0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  /*
     Create a parallel vector.
      - In this case, we specify the size of each processor's local
        portion, and PETSc computes the global size.  Alternatively,
        if we pass the global size and use PETSC_DECIDE for the
        local size PETSc will choose a reasonable partition trying
        to put nearly an equal number of elements on each processor.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,rank+1,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecGetSize(x,&N);CHKERRQ(ierr);

  /*
     Duplicate some work vectors (of the same format and
     partitioning as the initial vector).
  */
  ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&w);CHKERRQ(ierr);


  ierr = VecSet(x,one);CHKERRQ(ierr);
  ierr = VecSet(y,two);CHKERRQ(ierr);


  /* Vector norm */
  ierr = VecNorm(y,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of y: %g\n",(double)norm);CHKERRQ(ierr);



  /*
     Set the vector elements.
      - Always specify global locations of vector entries.
      - Each processor can contribute any vector entries,
        regardless of which processor "owns" them; any nonlocal
        contributions will be transferred to the appropriate processor
        during the assembly process.
      - In this example, the flag ADD_VALUES indicates that all
        contributions will be added together. Use INSERT_VALUES to
        replace existing contributins with new values.
  */
  for (i=0; i<N-rank; i++) {
    ierr = VecSetValues(x,1,&i,&one,ADD_VALUES);CHKERRQ(ierr);
  }


  /*
     Assemble vector, using the 2-step process:
       VecAssemblyBegin(), VecAssemblyEnd()
     Computations can be done while messages are in transition
     by placing code between these two statements.
  */
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);


  /*
      View the vector; then destroy it.
  */
  ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Finalize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
