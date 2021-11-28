
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
      KEc,//Local stiffness matrix of element_ij. (KEc)
      x,//Elemental density matrix. (x)
      xnew,//Updated elemental density matrix (xnew).
      dc,//Elemental sensitivity matrix (dx)
      dcn;//Filtered dc.
  KSP ksp;//Linear solver object.
  PetscInt nelx = 5, //Number of elements to distribute along x-direction.
           nely = 5, //Number of elements to distribute along y-direction.
           rmin = 3, //Filter radius.
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
            volfrac = 0.5,//Volume fraction.
            *coeffKEc;
  PetscReal change = 1,
              x_ij,//Variable to reference element_ij's density.
              c = 0,//Compliance of element_ij.
              C = 0;//Compliance of the whole structure.

  ierr = PetscInitialize(&argc,&args,NULL,help); if(ierr) return ierr;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"opt_","options for mesh","");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nelx","Number of x-elements","TopOpt.c",nelx,&nelx,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-nely","Number of y-elements","TopOpt.c",nely,&nely,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-E","Young's Modulus","TopOpt.c",E,&E,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-nu","Poisson's Ratio","TopOpt.c",nu,&nu,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-Fy","Force applied.","TopOpt.c",Fy,&Fy,NULL); CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  //Grid setup variables.40
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

ierr = MatCreate(PETSC_COMM_WORLD,&KEc); CHKERRQ(ierr);
ierr = MatSetSizes(KEc,PETSC_DECIDE,PETSC_DECIDE,8,8); CHKERRQ(ierr);
ierr = MatSetType(KEc,MATDENSE); CHKERRQ(ierr);



//*****************END ELEMENT STIFFNESS MATRIX (KE)*************************



//*****************BEGIN(OBJECT SETUP)*************************
//Elemental sensitivity matrix (dc)
ierr = MatCreate(PETSC_COMM_WORLD,&dc); CHKERRQ(ierr);
ierr = MatSetSizes(dc,PETSC_DECIDE,PETSC_DECIDE,nelx,nely); CHKERRQ(ierr);
ierr = MatSetType(dc,MATDENSE); CHKERRQ(ierr);
ierr = MatSetUp(dc); CHKERRQ(ierr);

//Filtered dc matrix (dcn)
ierr = MatCreate(PETSC_COMM_WORLD,&dcn); CHKERRQ(ierr);
ierr = MatSetSizes(dcn,PETSC_DECIDE,PETSC_DECIDE,nelx,nely); CHKERRQ(ierr);
ierr = MatSetType(dcn,MATDENSE); CHKERRQ(ierr);
ierr = MatSetUp(dcn); CHKERRQ(ierr);

//Elemental density matrix (x) [Initialize all with the input volfrac]
ierr = MatCreate(PETSC_COMM_WORLD,&x); CHKERRQ(ierr);
ierr = MatSetSizes(x,PETSC_DECIDE,PETSC_DECIDE,nelx,nely); CHKERRQ(ierr);
ierr = MatSetType(x,MATDENSE); CHKERRQ(ierr);
ierr = MatSetUp(x); CHKERRQ(ierr);
for(i=0;i<nely;i++){
  for(j=0;j<nelx;j++){
    ierr = MatSetValue(x,i,j,volfrac,INSERT_VALUES); CHKERRQ(ierr);
  }//end "j" loop
}//end "i" loop
ierr = MatAssemblyBegin(x,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
ierr = MatAssemblyEnd(x,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

//Updated elemental density matrix (xnew).
ierr = MatCreate(PETSC_COMM_WORLD,&xnew); CHKERRQ(ierr);
ierr = MatSetSizes(xnew,PETSC_DECIDE,PETSC_DECIDE,nelx,nely); CHKERRQ(ierr);
ierr = MatSetType(xnew,MATDENSE); CHKERRQ(ierr);
ierr = MatSetUp(xnew); CHKERRQ(ierr);

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
ierr = MatSetType(K,MATAIJ); CHKERRQ(ierr);
ierr = MatSetUp(K); CHKERRQ(ierr);

//Linear solver setup.
ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);//Create the solver.
//*****************END  (OBJECT SETUP)*************************


while(change > 0.01){
loop++;//Update loop iteration.
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
      ierr = MatGetValue(x,i,j,&x_ij); CHKERRQ(ierr);
      ierr = MatDuplicate(KE,MAT_COPY_VALUES,&KEc); CHKERRQ(ierr);
      PetscReal scale =  PetscPowScalar(x_ij,P);
      ierr = MatScale(KEc,scale); CHKERRQ(ierr);
      ierr = MatDenseGetArray(KEc,&coeffKEc);CHKERRQ(ierr);
      ierr = MatSetValues(K,8,edof,8,edof,coeffKEc,ADD_VALUES); CHKERRQ(ierr);

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

  //================BEGIN (COMPLIANCE & SENSITIVITIES)==========================
  PetscInt index[8] ={0,1,2,3,4,5,6,7};//indices for "select" below.
  PetscScalar select[8],//Store subset values here.
              dx_ij;//Stores sensitivity of x_ij of element_ij.
              C = 0; //Reset compliance counter.

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

      ierr = MatGetValue(x,i,j,&x_ij); CHKERRQ(ierr);//Get the density of element_ij.
      C += c*PetscPowScalar(x_ij,P);//Add the compliance of element_ij to the running total.
      dx_ij = -P*PetscPowScalar(x_ij,P-1)*c;//Compute the "natural" sensitivity of element_ij.
      ierr = MatSetValue(dc,j,i,dx_ij,INSERT_VALUES); CHKERRQ(ierr); //FOR WHATEVER REASON IT TRANSPOSES THE MATRIX IF I DO i,j instead of j,i
    }//end "j" loop.
  }//end "i" loop.

  ierr = MatAssemblyBegin(dc,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(dc,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatView(dc,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //ierr = VecView(U,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  //PetscPrintf(PETSC_COMM_WORLD,"Compliance C = %.6e\n",C);
  //================END (COMPLIANCE & SENSITIVITIES)==========================

  //==========================BEGIN  (FILTERING)===============================
  PetscScalar sum,fac,rad_kl,dcn_ji,x_lk,dc_lk,x_ji;
  PetscInt kk,l, //Looping variables.
           k1,
           k2,
           l1,
           l2;
  for(i=0;i<nelx;i++){
    for(j=0;j<nely;j++){
      sum = 0.0;
      dcn_ji = 0;//Reset dcn_ji.

      //NO handy max()/min() functions in PETSc.
      k1 = i-rmin;
      if(k1<0){k1 = 0;}//end if.
      k2 = i + rmin;
      if(k2>(nelx-1)){k2 = nelx-1;}//end if
      l1 = j - rmin;
      if(l1<0){l1 = 0;}//end if
      l2 = j + rmin;
      if(l2>(nely-1)){l2 = nely-1;}

      for(kk=k1;kk<=k2;kk++){
        for(l=l1;l<=l2;l++){
          rad_kl = PetscPowScalar((i-kk),2) + PetscPowScalar((j-l),2);
          fac = rmin - PetscPowScalar(rad_kl,0.5);
          if(fac>0){
            ierr = MatGetValues(x,1,&l,1,&kk,&x_lk); CHKERRQ(ierr);
            ierr = MatGetValues(dc,1,&l,1,&kk,&dc_lk); CHKERRQ(ierr);
            sum += fac;
            dcn_ji += x_lk*dc_lk*fac;
          }//end if
        }//end "l" loop.
      }//end "k" loop.
      ierr = MatGetValues(x,1,&j,1,&i,&x_ji); CHKERRQ(ierr);
      dcn_ji /= (sum*x_ji);
      ierr = MatSetValue(dcn,j,i,dcn_ji,INSERT_VALUES); CHKERRQ(ierr);
    }//end "j" loop.
  }//end "i" loop.
  ierr = MatAssemblyBegin(dcn,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(dcn,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatView(dcn,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); //(DEBUG)
  //============================END  (FILTERING)=========================


  //*****************BEGIN  (OPTIMALITY CRITERIA)*************************
  PetscReal lambda1=0,lambda2 =100000,lmid, m=0.2,sumx,xnew_ij;
  while((lambda2-lambda1)>0.0001){
    lmid = (lambda2+lambda1)/2;
    sumx = 0;
    for(i=0;i<nely;i++){
      for(j=0;j<nelx;j++){
        ierr = MatGetValue(x,i,j,&x_ij); CHKERRQ(ierr);
        ierr = MatGetValue(dcn,i,j,&dx_ij); CHKERRQ(ierr);
        xnew_ij = x_ij*PetscPowScalar(-dx_ij/lmid,0.5);
        if(xnew_ij>(x_ij+m)){
          xnew_ij = x_ij+m;
        }//end if
        if(xnew_ij>1){
          xnew_ij = 1;
        }//end if.
        if(xnew_ij<(x_ij-m)){
          xnew_ij = x_ij-m;
        }//end if.
        if(xnew_ij<0.001){
          xnew_ij = 0.001;
        }//end if.
        sumx += xnew_ij;
        ierr = MatSetValue(xnew,i,j,xnew_ij,INSERT_VALUES); CHKERRQ(ierr);
      }//end "j" loop.
    }//end "i" loop.

    if(sumx > (volfrac*nelx*nely)){
      lambda1 = lmid;
    }//end if.
    else{
      lambda2 = lmid;
    }//end else.
  }//end while
  ierr = MatAssemblyBegin(xnew,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(xnew,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatView(xnew,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); //(DEBUG)
  //*****************END  (OPTIMALITY CRITERIA)*************************
  //Compute change.
  ierr = MatAXPY(x,-1,xnew,SAME_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = MatNorm(x,NORM_INFINITY,&change);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"Iter: %3d\tC = %10.5f\tVol = %10.5f\t change = %10.5f \n",loop,C,sumx/(nelx*nely),change);

  //Update values for next iteration.
  ierr = MatDuplicate(dcn,MAT_COPY_VALUES,&dc); CHKERRQ(ierr);
  ierr = MatDuplicate(xnew,MAT_COPY_VALUES,&x); CHKERRQ(ierr);
}//end while.


  //Cleanup (free-up memory)
  KSPDestroy(&ksp);
  MatDestroy(&K); MatDestroy(&KE); MatDestroy(&KEc); MatDestroy(&dc);
  MatDestroy(&xnew); MatDestroy(&x); MatDestroy(&dcn);
  VecDestroy(&F); VecDestroy(&U); VecDestroy(&Ue); VecDestroy(&Uc);

  return PetscFinalize();
}//end main function
