MODULE OpenCMISSSetup

  USE OpenCMISS
  USE OPENCMISS_Iron

#ifndef NOMPIMOD
  USE MPI
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

CONTAINS

  SUBROUTINE CreateBasis(Basis,BasisUserNumber,BasisType,BasisInterpType,NumDim,Err)
  
    TYPE(CMFE_BasisType), INTENT(OUT) :: Basis
    INTEGER(CMISSIntg), INTENT(OUT) :: Err
    INTEGER(CMISSIntg), INTENT(IN) :: BasisUserNumber
    INTEGER(CMISSIntg), INTENT(IN) :: BasisType,BasisInterpType
    INTEGER(CMISSIntg), INTENT(IN) :: NumDim
    
    CALL cmfe_Basis_Initialise(Basis,Err)
    CALL cmfe_Basis_CreateStart(BasisUserNumber,Basis,Err)
    !set the basis to bilinear simplex
    CALL cmfe_Basis_TypeSet(Basis,BasisType,Err)
    CALL cmfe_Basis_NumberOfXiSet(Basis,NumDim,Err)
    IF(NumDim.EQ.1) THEN
      CALL cmfe_Basis_InterpolationXiSet(Basis,(/BasisInterpType/), Err)    
    ELSEIF(NumDim.EQ.2) THEN
      CALL cmfe_Basis_InterpolationXiSet(Basis,(/BasisInterpType, & 
      & BasisInterpType/), Err)
    ELSEIF(NumDim.EQ.3) THEN
      CALL cmfe_Basis_InterpolationXiSet(Basis,(/BasisInterpType, & 
      & BasisInterpType,BasisInterpType/), Err)
    ELSE
      WRITE(*,*) 'Invalid number of dimensions used'
    ENDIF

    !Finish the creation of the basis
    CALL cmfe_Basis_CreateFinish(Basis,Err)
  
    RETURN
  END SUBROUTINE CreateBasis
END MODULE OpenCMISSSetup
