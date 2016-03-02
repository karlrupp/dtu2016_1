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
  VecCreate(PETSC_COMM_WORLD,&x);
  VecGetSize(x,&N);

  VecSet(x,one);

  /*
     Set the vector elements.
  */
  VecGetArray(x, &local_x);
  VecGetLocalSize(x, &n);
  for (i=0; i<=n; i++) {
    local_x[i] = i;
  }
  VecRestoreArray(x, &local_x);


  /*
      View the vector; then destroy it.
  */
  VecView(x,PETSC_VIEWER_STDOUT_WORLD);
  VecDestroy(&x);


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Finalize program
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscFinalize();

  return 0;
}
