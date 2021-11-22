/* HelloWorldPETSc.c
 *
 * Simple file from petsc.wordpress.com to test if petsc is workable on this machine
 *
 * 
 */

//Includes
#include "petsc.h"

int main(int argc, char *argv[]){

	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
	PetscPrintf(PETSC_COMM_WORLD,"Hello World\n");

	//Some Variables
	PetscErrorCode ierr; //error tracker
	PetscViewer viewer; //petsc viewer for printing	
	
	PetscInt xD = 100,
		 yD = 100;
	PetscInt xD2, yD2, xDl, yDl;
	PetscMPIInt size, rank;
	
	Mat A;
	PetscInt i,j,k;
	PetscScalar st = 3, fp_out;
	PetscInt firstR,lastR,firstC,lastC;
	PetscInt ct = 0;	
	
	char fileName[80];

	//Viewer config
	//PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	//PetscViewerSetType(viewer, PETSCVIEWERASCII);
	//PetscViewerFileSetMode(viewer, FILE_MODE_CREATE);
	//PetscViewerFileSetName(viewer,"test.txt");

	//Rank and size
	MPI_Comm_size(PETSC_COMM_WORLD,&size);
	MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

	//Create matrices
	
	ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,xD,yD); CHKERRQ(ierr);
	ierr = MatSetType(A,MATDENSE); CHKERRQ(ierr);
	ierr = MatSetUp(A); CHKERRQ(ierr);
	ierr = MatGetSize(A, &xD2, &yD2); CHKERRQ(ierr); //global sizes to populate
	ierr = MatGetLocalSize(A, &xDl, &yDl); CHKERRQ(ierr); //local sizes to populate

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Global sizes %i %i\n",xD2,yD2);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Rank %i Local sizes %i %i\n",rank,xDl,yDl);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);	

	if(rank==1){ //populate matrix if master
	for(i=0;i<xD;i++){//row
	  for(j=0;j<yD;j++){//col
	    	//PetscPrintf(PETSC_COMM_WORLD,"row %i col %i",i,j);
		ierr = MatSetValue(A,i,j,0,INSERT_VALUES); CHKERRQ(ierr);
	  }//end "j" loop
	}//end "i" loop
	}
	
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	
	ierr = MatGetOwnershipRange(A,&firstR,&lastR);
	ierr = MatGetOwnershipRangeColumn(A,&firstC,&lastC);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Rank %i Local rows %i %i\n",rank,firstR,lastR);
	ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Rank %i Local cols %i %i\n",rank,firstC,lastC);
	ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);	


	//Printing
	for(k=0;k<xD2;k++){
         sprintf(fileName,"tmpDat/MatDat%06u.dat",k);
	 PetscPrintf(PETSC_COMM_WORLD,"Working on %s\n",fileName);
	
	 PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
	 PetscViewerASCIIOpen(PETSC_COMM_WORLD,fileName,&viewer);
	 PetscViewerASCIIPushSynchronized(viewer); //activate synchronized printing sesh	

	 //if(rank==0){	
	 for(i=firstR;i<lastR;i++){//row
	    for(j=0;j<yD2;j++){//col
		if(j==k){//matrix update portion (not needed for printing)
		  MatSetValue(A,i,j,st,INSERT_VALUES); //add a 3
		  if(j!=0){MatSetValue(A,i,j-1,0,INSERT_VALUES);}//remove old 3
		  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
		  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
		}
 		ierr = MatGetValues(A,1,&i,1,&j,&fp_out); CHKERRQ(ierr);
 		PetscViewerASCIISynchronizedPrintf(viewer,"%6.2f",fp_out);
 	   }//end for j - col
 	   PetscViewerASCIISynchronizedPrintf(viewer,"\n");
 	 }//end for i - row
	 //}
	
	 PetscViewerFlush(viewer); //when stuff is synchronized
	 PetscViewerASCIIPopSynchronized(viewer); //end synch sesh
	 PetscViewerDestroy(&viewer);
	}//end for k

	//end Printing
	
	char cmdStr[80];
	//if(rank==0){//only for master
	PetscInt step = (PetscInt) xD2/size;
	PetscInt begin = rank*step;
	PetscInt end = rank+1*step;
	if(rank==size){end = xD2;}

	 PetscPrintf(PETSC_COMM_WORLD,"CreatingImages\n");
         for(ct=begin;ct<end;ct++){
	  sprintf(cmdStr,"gnuplot -e \"a=%i;load 'gpScripts/genAniMain.gp'\"",ct);
          system(cmdStr);
	 }
	PetscBarrier(NULL);
	if(rank==0){
         PetscPrintf(PETSC_COMM_WORLD,"Compiling gif\n");
         system("convert -delay 5 tmpAni/heat*.jpg -loop 0 aniOut.gif");
	}
	
	PetscBarrier(NULL);
	PetscFinalize();

	return 0;
}
