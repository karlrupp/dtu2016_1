static const char help[] = "A help string, printed on -help\n\n";

#include "petscdmda.h"
#include "petscsnes.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Vec                x;
  PetscErrorCode     ierr;
  PetscMPIInt        rank;
  PetscScalar        one = 1.0;
  PetscScalar        *local_x;
  PetscInt           i,n,N;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Initialize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  PetscInitialize(&argc,&argv,0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  /*
     Create a parallel vector.
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x, rank+1, PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecGetSize(x,&N);CHKERRQ(ierr);

  ierr = VecSet(x,one);CHKERRQ(ierr);

  /*
     Set the vector elements.
  */
  ierr = VecGetArray(x, &local_x);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x, &n);CHKERRQ(ierr);
  for (i=0; i<=n; i++) {
    local_x[i] = i;
  }
  ierr = VecRestoreArray(x, &local_x);CHKERRQ(ierr);


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
