
/*
Laplacian in 3D. Modeled by the partial differential equation

   - Laplacian u = 1,   0 < x,y,z < 1,

with boundary conditions

   u = 1 for x = 0, x = 1, y = 0, y = 1, z = 0, z = 1.

   This is able to use multigrid to solve the linear system

*/

static char help[] = "Solves 3D Laplacian using multigrid.\n\n";

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>

extern PetscErrorCode ComputeMatrix(KSP,Mat,Mat,void*);
extern PetscErrorCode ComputeRHS(KSP,Vec,void*);
extern PetscErrorCode ComputeInitialGuess(KSP,Vec,void*);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  PetscErrorCode ierr;
  KSP            ksp;
  PetscReal      norm;
  DM             da;
  Vec            x,b,r;
  Mat            A;
  PetscInt       its;

  PetscInitialize(&argc,&argv,(char*)0,help);

  /* Create DMDA and KSP and attach to each other */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DMDA_STENCIL_STAR,-7,-7,-7,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,0,0,0,&da);CHKERRQ(ierr);
  ierr = KSPSetDM(ksp,da);CHKERRQ(ierr);

  /* Register callback functions for 
       - computing the initial guess
       - computing the right hand side vector
       - computing the system matrix
     on each level. */
  ierr = KSPSetComputeInitialGuess(ksp,ComputeInitialGuess,NULL);CHKERRQ(ierr);
  ierr = KSPSetComputeRHS(ksp,ComputeRHS,NULL);CHKERRQ(ierr);
  ierr = KSPSetComputeOperators(ksp,ComputeMatrix,NULL);CHKERRQ(ierr);
  ierr = DMDestroy(&da);CHKERRQ(ierr);

  /* Solve the system */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,NULL,NULL);CHKERRQ(ierr);

  /* Compute the error */
  ierr = KSPGetSolution(ksp,&x);CHKERRQ(ierr);
  ierr = KSPGetRhs(ksp,&b);CHKERRQ(ierr);
  ierr = VecDuplicate(b,&r);CHKERRQ(ierr);
  ierr = KSPGetOperators(ksp,&A,NULL);CHKERRQ(ierr);
  ierr = MatMult(A,x,r);CHKERRQ(ierr);
  ierr = VecAXPY(r,-1.0,b);CHKERRQ(ierr);
  ierr = VecNorm(r,NORM_2,&norm);CHKERRQ(ierr);

  /*
     Print convergence information.  PetscPrintf() produces a single
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
  */
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %g iterations %D\n",(double)norm,its);CHKERRQ(ierr);

  /* Cleanup */
  ierr = VecDestroy(&r);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = PetscFinalize();

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
PetscErrorCode ComputeRHS(KSP ksp,Vec b,void *ctx)
{
  PetscErrorCode ierr;
  PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
  DM             dm;
  PetscScalar    Hx,Hy,Hz,HxHydHz,HyHzdHx,HxHzdHy;
  PetscScalar    ***barray;

  PetscFunctionBeginUser;
  ierr    = KSPGetDM(ksp,&dm);CHKERRQ(ierr);
  ierr    = DMDAGetInfo(dm,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  Hx      = 1.0 / (PetscReal)(mx-1); Hy = 1.0 / (PetscReal)(my-1); Hz = 1.0 / (PetscReal)(mz-1);
  HxHydHz = Hx*Hy/Hz; HxHzdHy = Hx*Hz/Hy; HyHzdHx = Hy*Hz/Hx;
  ierr    = DMDAGetCorners(dm,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);
  ierr    = DMDAVecGetArray(dm,b,&barray);CHKERRQ(ierr);

  /* TODO: Fill RHS */


  ierr = DMDAVecRestoreArray(dm,b,&barray);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeInitialGuess"
PetscErrorCode ComputeInitialGuess(KSP ksp,Vec b,void *ctx)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = VecSet(b,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMatrix"
PetscErrorCode ComputeMatrix(KSP ksp,Mat jac,Mat B,void *ctx)
{
  DM             da;
  PetscErrorCode ierr;
  PetscInt       i,j,k,mx,my,mz,xm,ym,zm,xs,ys,zs;
  PetscScalar    v[7],Hx,Hy,Hz,HxHydHz,HyHzdHx,HxHzdHy;
  MatStencil     row,col[7];

  PetscFunctionBeginUser;
  ierr    = KSPGetDM(ksp,&da);CHKERRQ(ierr);
  ierr    = DMDAGetInfo(da,0,&mx,&my,&mz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
  Hx      = 1.0 / (PetscReal)(mx-1); Hy = 1.0 / (PetscReal)(my-1); Hz = 1.0 / (PetscReal)(mz-1);
  HxHydHz = Hx*Hy/Hz; HxHzdHy = Hx*Hz/Hy; HyHzdHx = Hy*Hz/Hx;
  ierr    = DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm);CHKERRQ(ierr);

  /* TODO: Fill matrix B */ 

  ierr   = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr   = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
