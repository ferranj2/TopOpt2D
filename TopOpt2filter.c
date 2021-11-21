
static char help[] = "PETSc version of O.Sigmund's 1999 MATLAB code.\n"
"Option prefix = opt_.\n";
#include <petsc.h>
//Executable callout:
//./TopOpt -ksp_view_solution -mat_view draw -draw_pause 100

//Future upgrades: Isolate the matrix builder as its own function.

int main(int argc, char **args){
  //Variable declarations
  PetscErrorCode ierr;//Error tracker
  //PetscViewer viewer;
  Vec F, //Force vector (F)
      U, //Gloval displacement vector (U)
      Ue, //Displacement vector of element_ij #1.
      Uc; //Displacement vector of element_ij #2.
  Mat K,//Global stiffness matrix. (K)
      KE,//Local stiffness matrix. (KE) [Common to ALL elements]
      x,//Elemental density matrix. (x)
      dc//Elemental sensitivity matrix (dx)
      dcn;//Filtered dc.
  KSP ksp;//Linear solver object.
  PetscInt nelx = 5, //Number of elements to distribute along x-direction.
           nely = 5, //Number of elements to distribute along y-direction.
           rmin = 2, //Filter radius.
           nx, //Number of nodes distributed along x-direction.
           ny, //Number of nodes distributed along y-direction.
           size,//Size of global stiffness matrix K (K = size by size)
           sizeKE = 8,//Size of local stiffness matrix KE (KE = sizeKE by sizeKE)
           loop = 0;// TopOpt iteration counter.
  PetscInt i, j, w;//Looping variables
  PetscInt n1,n2,last,edof[8];//Dedicated variables for array indexing.
  PetscInt LIFy[1] = {1};//LI of the dof where Fy is to act.
  PetscReal E = 1, //Young's Modulus
            nu = 0.3, //Poisson's ratio.
            Fy = -1, //Force applied (nondimendional)
            P = 3, //Penalization coefficient.
            volfrac = 0.3;//Volume fraction.
  PetscScalar c = 0,//Compliance of element_ij.
              C = 0;//Compliance of the whole structure.

  ierr = PetscInitialize(&argc,&args,NULL,help); if(ierr) return ierr;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"opt_","options for mesh","");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nelx","Number of x-elements","TopOpt.c",nelx,&nelx,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nely","Number of y-elements","TopOpt.c",nely,&nely,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-E","Young's Modulus","TopOpt.c",E,&E,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-nu","Poisson's Ratio","TopOpt.c",nu,&nu,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-Fy","Force applied.","TopOpt.c",Fy,&Fy,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  //Grid setup variables.
  nx = nelx + 1;//Number of X-nodes.
  ny = nely + 1;//Number of Y-nodes.
  size = 2*nx*ny; //Size of the global stiffness matrix.
  last = size - 1;//LI for the last row or column of the global stiffness matrix.


  //*****************BEGIN ELEMENT STIFFNESS MATRIX (KE)*************************
  //Coefficients of the local stiffness matrix.
  PetscReal mult = E/(1-nu*nu);//A common factor shared by the element stifnesses.
  PetscReal k[8] = {+(3-nu)*mult/6,
          +(1+nu)*mult/8,
          -(3+nu)*mult/12,
          +(3*nu-1)*mult/8,
          +(nu-3)*mult/12,
          -(nu+1)*mult/8,
          +nu*mult/6,
          +(1-3*nu)*mult/8};
  PetscReal coeff[8][8] = {{k[0],k[1],k[2],k[3],k[4],k[5],k[6],k[7]},
                 {k[1],k[0],k[7],k[6],k[5],k[4],k[3],k[2]},
                 {k[2],k[7],k[0],k[5],k[6],k[3],k[4],k[1]},
                 {k[3],k[6],k[5],k[0],k[7],k[2],k[1],k[4]},
                 {k[4],k[5],k[6],k[7],k[0],k[1],k[2],k[3]},
                 {k[5],k[4],k[3],k[2],k[1],k[0],k[7],k[6]},
                 {k[6],k[3],k[4],k[1],k[2],k[7],k[0],k[5]},
                 {k[7],k[2],k[1],k[4],k[3],k[6],k[5],k[0]}};
ierr = MatCreate(PETSC_COMM_WORLD,&KE); CHKERRQ(ierr);
ierr = MatSetSizes(KE,PETSC_DECIDE,PETSC_DECIDE,8,8); CHKERRQ(ierr);
ierr = MatSetType(KE,MATDENSE); CHKERRQ(ierr);
ierr = MatSetUp(KE); CHKERRQ(ierr);
ierr = MatDensePlaceArray(KE,*coeff);CHKERRQ(ierr);
ierr = MatAssemblyBegin(KE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
ierr = MatAssemblyEnd(KE,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
//*****************END ELEMENT STIFFNESS MATRIX (KE)*************************



//*****************BEGIN(OBJECT SETUP)*************************
//Elemental sensitivity matrix (dc)
ierr = MatCreate(PETSC_COMM_WORLD,&dc); CHKERRQ(ierr);
ierr = MatSetSizes(dc,PETSC_DECIDE,PETSC_DECIDE,nelx,nely); CHKERRQ(ierr);
ierr = MatSetType(dc,MATDENSE); CHKERRQ(ierr);
ierr = MatSetUp(dc); CHKERRQ(ierr);

//Filtered dc matrix
ierr = MatDuplicate(dc,MAT_DO_NOT_COPY_VALUES,*dcn); CHKERRQ(ierr);

//Elemental density matrix (x) [Initialize all with the input volfrac]
ierr = MatCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
ierr = MatSetSizes(x,PETSC_DECIDE,PETSC_DECIDE,nelx,nely); CHKERRQ(ierr);
ierr = MatSetType(x,MATDENSE); CHKERRQ(ierr);
ierr = MatSetUp(x); CHKERRQ(ierr);
for(i=0;j<nelx;i++){
  for(j=0;i<nely;j++){
    ierr = MatSetValue(x,i,j,volfrac,INSERT_VALUES); CHKERRQ(ierr);
  }//end "j" loop
}//end "i" loop

//Force vector (RHS in the FEA subroutine).
ierr = VecCreate(PETSC_COMM_WORLD,&F); CHKERRQ(ierr);
ierr = VecSetSizes(F,PETSC_DECIDE,size); CHKERRQ(ierr);
ierr = VecSetFromOptions(F); CHKERRQ(ierr);
ierr = VecSetValues(F,1,LIFy,&Fy,INSERT_VALUES); CHKERRQ(ierr);//Apply the point force on the beam.
ierr = VecAssemblyBegin(F); CHKERRQ(ierr);
ierr = VecAssemblyEnd(F); CHKERRQ(ierr);

//Displacement vector (same dimensions as F vector, entries zeroed later)
ierr = VecDuplicate(F,&U); CHKERRQ(ierr);

//Subset vector (for elementwise compliance computations).
ierr = VecCreate(PETSC_COMM_WORLD,&Ue); CHKERRQ(ierr);
ierr = VecSetSizes(Ue,PETSC_DECIDE,sizeKE); CHKERRQ(ierr);
ierr = VecSetType(Ue,VECSTANDARD); CHKERRQ(ierr);

//Subset vector (for elementwise compliance computations).
ierr = VecCreate(PETSC_COMM_WORLD,&Uc); CHKERRQ(ierr);
ierr = VecSetSizes(Uc,PETSC_DECIDE,sizeKE); CHKERRQ(ierr);
ierr = VecSetType(Uc,VECSTANDARD); CHKERRQ(ierr);

//Global stiffness matrix.
ierr = MatCreate(PETSC_COMM_WORLD,&K); CHKERRQ(ierr);
ierr = MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,size,size); CHKERRQ(ierr);
ierr = MatSetType(K,MATDENSE); CHKERRQ(ierr);
ierr = MatSetUp(K); CHKERRQ(ierr);

//Linear solver setup.
ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);//Create the solver.
//*****************END  (OBJECT SETUP)*************************



//FINITE ELEMENT ANALYSIS
//*****************BEGIN(FEA)*************************
  //Global stiffness matrix (A matrix in the FEA subroutine)
  for(i=0;i<nely;i++){//scan rows of elements.
    for(j=0;j<nelx;j++){//scan columns of elements.
      //LI of nodes in the mesh (dof's are Ux,Uy vertical pairs in column-major order)
      n1 = j + i*ny;//LI of TL corner of element_ij in stencil.
      n2 = n1 + ny;//LI of TR corner of element_ij in stencil.
      edof[0] = 2*n1;//LI of X-displacement of TL corner in U.
      edof[1] = edof[0]+1;//LI of Y-displacement of TL corner in U.
      edof[2] = 2*n2;//LI X-displacement of TR corner in U.
      edof[3] = edof[2]+1;//LI Y-displacement of TR corner in U.
      edof[4] = edof[3]+1;//LI X-displacement of BL corner in U.
      edof[5] = edof[4]+1;//LI Y-displacement of BL corner in U.
      edof[6] = edof[1]+1;//LI X-displacement of BR corner in U.
      edof[7] = edof[6]+1;//LI Y-displacement of BR corner in U.

      for(w=0;w<8;w++){
        ierr = MatSetValues(K,1,&edof[w],8,edof,coeff[w],ADD_VALUES); CHKERRQ(ierr);
      }//end "w" loop.

    }//end "j" loop
  }//end "i" loop
  ierr = MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  //Dirichlet BC: X-displacements clamped on left edge.
  for(w=0;w<ny;w++){//Scan along the y-nodes
    n1 = w*2;//"n1" now stores the row index of the clamped X-nodes.
    ierr = MatZeroRows(K,1,&n1,1,NULL,NULL); CHKERRQ(ierr);
  }//end "w" loop

  //Dirichlet BC: Y-displacement clamped at BR corner
  ierr = MatZeroRows(K,1,&last,1,NULL,NULL); CHKERRQ(ierr);

  //Callout to FEA solver.
  ierr = KSPSetOperators(ksp,K,K); CHKERRQ(ierr);//Set "A" matrix and preconditioner.
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);//Allow for optional user settings.
  ierr = KSPSolve(ksp,F,U); CHKERRQ(ierr);//Solve K*U = F.
  //*****************END  (FEA)*************************





  PetscInt index[8] ={0,1,2,3,4,5,6,7};//indices for "select" below.
  PetscScalar select[8],//Store subset values here.
              x_ij,//Stores elemental density of element_ij.
              dx_ij;//Stores sensitivity of x_ij of element_ij.

  //Elementwise measurement of the compliance.
  for(i=0;i<nely;i++){//Scan elements row-wise
    for(j=0;j<nelx;j++){//At current row scan elements column-wise.
      n1 = j + i*ny;//LI of TL corner of element_ij in stencil.
      n2 = n1 + ny;//LI of TR corner of element_ij in stencil.
      edof[0] = 2*n1;//LI of X-displacement of TL corner in U.
      edof[1] = edof[0]+1;//LI of Y-displacement of TL corner in U.
      edof[2] = 2*n2;//LI X-displacement of TR corner in U.
      edof[3] = edof[2]+1;//LI Y-displacement of TR corner in U.
      edof[4] = edof[3]+1;//LI X-displacement of BL corner in U.
      edof[5] = edof[4]+1;//LI Y-displacement of BL corner in U.
      edof[6] = edof[1]+1;//LI X-displacement of BR corner in U.
      edof[7] = edof[6]+1;//LI Y-displacement of BR corner in U.

      ierr = VecGetValues(U,8,edof,select); CHKERRQ(ierr);//Extract values from the global displacement vector.
      ierr = VecSetValues(Ue,8,index,select,INSERT_VALUES); CHKERRQ(ierr);//Build the displacement vector Ue of element_ij.
      ierr = VecAssemblyBegin(Ue); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(Ue); CHKERRQ(ierr);
      ierr = MatMult(KE,Ue,Uc); CHKERRQ(ierr);
      ierr = VecDot(Ue,Uc,&c); CHKERRQ(ierr);

      ierr = MatGetValues(K,1,&i,1,&j,&x_ij); CHKERRQ(ierr);
      C += c*PetscPowScalar(x_ij,P);//Add the compliance of element_ij to the running total.
      dx_ij = -P*PetscPowScalar(x_ij,P-1)*c;//Compute the "natural" sensitivity of element_ij.
      ierr = MatSetValue(dc,i,j,dx_ij,INSERT_VALUES); CHKERRQ(ierr);
    }//end "j" loop.
  }//end "i" loop.

  PetscPrintf(PETSC_COMM_WORLD,"Compliance C = %.6e",C);
  C = 0; //Reset compliance counter.

  //*****************BEGIN  (FILTERING)*************************
  PetscScalar sum,fac,rad_kl,dcn_ji,x_lk,dc_lk;
  PetscInt k,l, //Looping variables.
           k1,
           k2,
           l1,
           l2;
  for(i=0;i<nelx;i++){
    for(j=0;j<nely;j++){
      sum = 0;

      //NO handy max()/min() functions in PETSc.
      k1 = i-rmin;
      if(k1<1){k1 = 1;}//end if.
      k2 = i + rmin;
      if(k2>nelx){k2 = nelx;}//end if
      l1 = j - rmin;
      if(l1<1){l1 = 1;}//end if
      l2 = j + rmin;
      if(l2>nely){l2 = nely;}

      for(k=k1;k<k2;k++){
        for(l=l1;l<l2;l++){
          rad_kl = PetscPowScalar((i-k)*(i-k),2) + PetscPowScalar((j-l)*(j-l),2);
          fac = rmin - PetscPowScalar(rad_kl,0.5);
          if(fac>0){
            ierr = MatGetValues(x,1,&l,1,&k,&x_lk); CHKERRQ(ierr);
            ierr = MatGetValues(dc,1,&l,1,&k,&dc_lk); CHKERRQ(ierr);
            sum+=fac;
            dcn_ji += x_lk*dc_lk*fac;
          }
        }//end "l" loop.
      }//end "k" loop.
      dcn_ji /= sum*x_ji;
      ierr = MatSetValue(dcn,j,i,dcn_ji,INSERT_VALUES); CHKERRQ(ierr);
    }//end "j" loop.
  }//end "i" loop.
  //*****************END  (FILTERING)*************************

  //*****************BEGIN  (OPTIMALITY CRITERIA)*************************
  //*****************END  (OPTIMALITY CRITERIA)*************************



  //Cleanup (free-up memory)
  KSPDestroy(&ksp);
  MatDestroy(&K); MatDestroy(&KE); MatDestroy(&dc); MatDestroy(&x);
  VecDestroy(&F); VecDestroy(&U); VecDestroy(&Ue); VecDestroy(&Uc);

  return PetscFinalize();
}//end main function
