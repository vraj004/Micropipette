!> \file
!> $Id: DiffusionExample.f90 1528 2010-09-21 01:32:29Z chrispbradley $
!> \author Chris Bradley
!> \brief This is an example program to solve a diffusion equation using OpenCMISS calls.
!>
!> \section LICENSE
!>
!> Version: MPL 1.1/GPL 2.0/LGPL 2.1
!>
!> The contents of this file are subject to the Mozilla Public License
!> Version 1.1 (the "License"); you may not use this file except in
!> compliance with the License. You may obtain a copy of the License at
!> http://www.mozilla.org/MPL/
!>
!> Software distributed under the License is distributed on an "AS IS"
!> basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
!> License for the specific language governing rights and limitations
!> under the License.
!>
!> The Original Code is OpenCMISS
!>
!> The Initial Developer of the Original Code is University of Auckland,
!> Auckland, New Zealand and University of Oxford, Oxford, United
!> Kingdom. Portions created by the University of Auckland and University
!> of Oxford are Copyright (C) 2007 by the University of Auckland and
!> the University of Oxford. All Rights Reserved.
!>
!> Contributor(s):
!>
!> Alternatively, the contents of this file may be used under the terms of
!> either the GNU General Public License Version 2 or later (the "GPL"), or
!> the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
!> in which case the provisions of the GPL or the LGPL are applicable instead
!> of those above. If you wish to allow use of your version of this file only
!> under the terms of either the GPL or the LGPL, and not to allow others to
!> use your version of this file under the terms of the MPL, indicate your
!> decision by deleting the provisions above and replace them with the notice
!> and other provisions required by the GPL or the LGPL. If you do not delete
!> the provisions above, a recipient may use your version of this file under
!> the terms of any one of the MPL, the GPL or the LGPL.
!>

!> \example ClassicalField/ReactionDiffusion/ReactionDiffusionNoSource1D/src/ReactionDiffusionNoSource1DExample.f90
!! Example program to solve a diffusion equation using OpenCMISS calls.
!! \htmlinclude ClassicalField/ReactionDiffusion/ReactionDiffusionNoSource1D/history.html
!<

!> Main program
PROGRAM ActinWaves

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE OpenCMISSSetup
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif


  !Program Parameters
  REAL(CMISSRP), PARAMETER ::  store_coeff = 1.0_CMISSRP
  !Boundary Markers
  INTEGER(CMISSIntg), PARAMETER ::  MEMBRANE_MARKER=1
  INTEGER(CMISSIntg), PARAMETER ::  CYTOSOL_MARKER=0

  INTEGER(CMISSIntg), PARAMETER :: CoordinateSystemUserNumber=1
  INTEGER(CMISSIntg), PARAMETER :: RegionUserNumber=2
  INTEGER(CMISSIntg), PARAMETER :: BasisUserNumber=3

  INTEGER(CMISSIntg), PARAMETER :: MeshUserNumber=5
  INTEGER(CMISSIntg), PARAMETER :: DecompositionUserNumber=6
  INTEGER(CMISSIntg), PARAMETER :: GeometricFieldUserNumber=7
  INTEGER(CMISSIntg), PARAMETER :: EquationsSetFieldUserNumber=8

  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveFieldUserNumber=9
  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveMaterialsFieldUserNumber=10
  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveEquationsSetUserNumber=11
  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveEquationsSetFieldUserNumber=23
  INTEGER(CMISSIntg), PARAMETER :: NPF_ActiveSourceFieldUserNumber=12

  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveFieldUserNumber=13
  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveMaterialsFieldUserNumber=14
  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveEquationsSetUserNumber=15
  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveEquationsSetFieldUserNumber=24
  INTEGER(CMISSIntg), PARAMETER :: NPF_InActiveSourceFieldUserNumber=16

  INTEGER(CMISSIntg), PARAMETER :: FActinFieldUserNumber=25
  INTEGER(CMISSIntg), PARAMETER :: stwoFieldUserNumber=26
  INTEGER(CMISSIntg), PARAMETER :: soneFieldUserNumber=27
  INTEGER(CMISSIntg), PARAMETER :: kzeroFieldUserNumber=28
  INTEGER(CMISSIntg), PARAMETER :: hFieldUserNumber=29

  INTEGER(CMISSIntg), PARAMETER :: ProblemUserNumber=17

  INTEGER(CMISSIntg), PARAMETER :: CellMLUserNumber=18
  INTEGER(CMISSIntg), PARAMETER :: CellMLModelsFieldUserNumber=19
  INTEGER(CMISSIntg), PARAMETER :: CellMLStateFieldUserNumber=20
  INTEGER(CMISSIntg), PARAMETER :: CellMLIntermediateFieldUserNumber=21
  INTEGER(CMISSIntg), PARAMETER :: CellMLParametersFieldUserNumber=22


  !MECHANICS VARIABLES
  INTEGER(CMISSIntg), PARAMETER :: MechStiffnessFieldUserNumber=30
  INTEGER(CMISSIntg), PARAMETER :: DeformedFieldUserNumber=31
  INTEGER(CMISSIntg), PARAMETER :: MechEquationSetUserNumber=32
  INTEGER(CMISSIntg), PARAMETER :: MechEquationsSetFieldUserNumber=33
  INTEGER(CMISSIntg), PARAMETER :: MechProblemUserNumber=34
  INTEGER(CMISSIntg), PARAMETER :: FibreFieldUserNumber=35
  INTEGER(CMISSIntg), PARAMETER :: PressureBasisUserNumber=36
    
  
  !cmfe type variables used for this program

  !mesh/domain related fields
  TYPE(cmfe_BasisType) :: Basis
  TYPE(cmfe_CoordinateSystemType) :: CoordinateSystem,WorldCoordinateSystem
  TYPE(cmfe_DecompositionType) :: Decomposition
  TYPE(cmfe_FieldType) :: GeometricField
  TYPE(cmfe_FieldsType) :: Fields
  TYPE(cmfe_RegionType) :: Region,WorldRegion
  TYPE(cmfe_MeshType) :: Mesh
  TYPE(cmfe_MeshElementsType) :: MeshElements
  TYPE(cmfe_NodesType) :: Nodes
  !simulation equations fields
  TYPE(cmfe_EquationsType) :: NPF_ActiveEquations, NPF_InActiveEquations
  TYPE(cmfe_EquationsSetType) :: NPF_ActiveEquationsSet, NPF_InActiveEquationsSet
  TYPE(cmfe_FieldType) :: NPF_ActiveEquationsSetField, NPF_InActiveEquationsSetField
  TYPE(cmfe_FieldType) :: NPF_ActiveField, NPF_InActiveField, FActinField,stwoField,soneField
  TYPE(cmfe_FieldType) :: kzeroField,NPF_ActiveMaterialsField, NPF_InActiveMaterialsField,hField
  TYPE(cmfe_FieldType) :: NPF_ActiveSourceField, NPF_InActiveSourceField


  !equations solving related fields
  TYPE(cmfe_BoundaryConditionsType) :: BoundaryConditions
  TYPE(cmfe_ProblemType) :: Problem
  TYPE(cmfe_ControlLoopType) :: ControlLoop
  TYPE(cmfe_SolverType) :: Solver, LinearSolver
  TYPE(cmfe_SolverEquationsType) :: SolverEquations

  !CellML related cmiss fields
  TYPE(cmfe_CellMLType) :: CellML
  TYPE(cmfe_CellMLEquationsType) :: CellMLEquations
  TYPE(cmfe_FieldType) :: CellMLModelsField,CellMLStateField,CellMLIntermediateField,CellMLParametersField

  !Mechanics Related Fields
  TYPE(cmfe_MeshElementsType) :: PressureElements
  TYPE(cmfe_BasisType) :: PressureBasis
  TYPE(cmfe_BoundaryConditionsType) :: MechBCs
  TYPE(cmfe_ProblemType) :: MechProblem  
  TYPE(cmfe_ControlLoopType) :: LoadLoop
  TYPE(cmfe_EquationsType) :: MechEquations
  TYPE(cmfe_EquationsSetType) :: MechEquationsSet
  TYPE(cmfe_FieldType) :: FibreField,MechStiffnessField,DeformedField,MechEquationsSetField
  TYPE(cmfe_SolverType) :: MechSolver,MechLinearSolver
  TYPE(cmfe_SolverEquationsType) :: MechSolverEquations
  
  
#ifdef WIN32
  !Quickwin type
  LOGICAL :: QUICKWIN_STATUS=.FALSE.
  TYPE(WINDOWCONFIG) :: QUICKWIN_WINDOW_CONFIG
#endif
  
  !Other program variables
  REAL(CMISSRP),ALLOCATABLE,DIMENSION(:,:) :: NodeCoords
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: ElemMap
  INTEGER(CMISSIntg),ALLOCATABLE,DIMENSION(:,:) :: NodeNums

  INTEGER(CMISSIntg) :: NUMBER_OF_ELEMENTS,CONDITION
  INTEGER(CMISSIntg) :: NUMBER_OF_DOMAINS,NODE_NUMBER
  INTEGER(CMISSIntg),DIMENSION(12) :: MechBCNODES_PULL
  INTEGER(CMISSIntg),DIMENSION(6) :: MechBCNODES
  
  INTEGER(CMISSIntg) :: MPI_IERROR
  INTEGER :: node,i,st,outputfreq,NUMBER_OF_NODES,NUMBER_OF_ATTRIBUTES,NUMBER_OF_COORDS, &
    & BOUNDARY_MARKER,nodedomain,NODES_PER_ELE,ELE_ATTRIBUTES,element,GeometricMeshComponent
  REAL(CMISSRP) :: init_NPFActive, init_NPFInActive, Dx_NPFActive,Dy_NPFActive,Dx_NPFInActive,Dy_NPFInActive
  REAL(CMISSRP) :: startT,endT,Tstep,ODE_TIME_STEP,nodex,nodey, VALUE, init_FActin,kzero,stwo,sone
  REAL(CMISSRP) :: additive_perturb_val, perturb_pos_startx,perturb_NPFActive,perturb_NPFInActive
  INTEGER(CMISSIntg) :: NumberOfComputationalNodes,ComputationalNodeNumber
  INTEGER(CMISSIntg) :: EquationsSetIndex,MechEquationsSetIndex,CellMLIndex,constantModelIndex
  INTEGER(CMISSIntg) :: Err
  LOGICAL :: EXPORT_FIELD,INPIPETTE
  CHARACTER(250) :: CELLID,NODEFILE,ELEMFILE,ActinPolymSignalModel  
#ifdef WIN32
  !Initialise QuickWin
  QUICKWIN_WINDOW_CONFIG%TITLE="General Output" !Window title
  QUICKWIN_WINDOW_CONFIG%NUMTEXTROWS=-1 !Max possible number of rows
  QUICKWIN_WINDOW_CONFIG%MODE=QWIN$SCROLLDOWN
  !Set the window parameters
  QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
  !If attempt fails set with system estimated values
  IF(.NOT.QUICKWIN_STATUS) QUICKWIN_STATUS=SETWINDOWCONFIG(QUICKWIN_WINDOW_CONFIG)
#endif
  

!_________________________________________________________________________________________________
  !Problem INPUTS. PARAMETERS FROM Holmes et.al. 2012
  !MESH FILES
  open(unit=9,file='inputs.txt',status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening inputs file',st
    STOP
  ELSE
    PRINT *,'inputs file opened correctly'
    READ(9,*) !file info
    READ(9,*) !comment on next set of variables - reading in file info
    READ(9,*) NODEFILE,ELEMFILE
    READ(9,*)
    READ(9,*) ActinPolymSignalModel
    READ(9,*)
    READ(9,*) init_NPFActive,Dx_NPFActive,Dy_NPFActive
    READ(9,*)
    READ(9,*) init_NPFInActive,Dx_NPFInActive,Dy_NPFInActive
    READ(9,*)
    READ(9,*) init_FActin
    READ(9,*)
    READ(9,*) sone,stwo,kzero
    READ(9,*)
    READ(9,*) additive_perturb_val,perturb_pos_startx
    READ(9,*)
    READ(9,*) startT,endT,Tstep,ODE_TIME_STEP
    READ(9,*) 
    READ(9,*) outputfreq
  ENDIF
  CLOSE(9)
  !Write the params out to screen for double checking

  WRITE(*,*) 'Node file:',NODEFILE
  WRITE(*,*) 'Element file:',ELEMFILE
  WRITE(*,*) 'CellML Model File:', ActinPolymSignalModel
  !cell initial conditions
  WRITE(*,*) 'Initial [NPF_Active]i = ',init_NPFActive !dependent field (cytosolic NPFActive).
  WRITE(*,*) 'NPF_Active Diff Coeff in x = ',Dx_NPFActive
  WRITE(*,*) 'NPF_Active Diff Coeff in y = ',Dy_NPFActive
  WRITE(*,*) 'Initial [NPF_InActive]i = ',init_NPFInActive !dependent field (cytosolic NPFInActive).
  WRITE(*,*) 'NPF_InActive Diff Coeff in x = ',Dx_NPFInActive
  WRITE(*,*) 'NPF_InActive Diff Coeff in y = ',Dy_NPFInActive

  WRITE(*,*) 'Tstart=',startT
  WRITE(*,*) 'Tend=',endT
  WRITE(*,*) 'Tstep=',Tstep
  WRITE(*,*) 'ODE_Tstep=',ODE_TIME_STEP
  WRITE(*,*) 'Output Frequency=',outputfreq


  
  EXPORT_FIELD=.TRUE.
!_________________________________________________________________________________________________
  !Intialise OpenCMISS
  CALL cmfe_Initialise(WorldCoordinateSystem,WorldRegion,Err)
  CALL cmfe_ErrorHandlingModeSet(CMFE_ERRORS_TRAP_ERROR,Err)
  !Get the computational nodes information
  CALL cmfe_ComputationalNumberOfNodesGet(NumberOfComputationalNodes,Err)
  CALL cmfe_ComputationalNodeNumberGet(ComputationalNodeNumber,Err)
  NUMBER_OF_DOMAINS=NumberOfComputationalNodes

  !Set all diganostic levels on for testing
  CALL MPI_BCAST(NUMBER_OF_DOMAINS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  !Start the creation of a new RC coordinate system
  CALL cmfe_CoordinateSystem_Initialise(CoordinateSystem,Err)
  CALL cmfe_CoordinateSystem_CreateStart(CoordinateSystemUserNumber,CoordinateSystem,Err)
  !Set the coordinate system to be 2D
  CALL cmfe_CoordinateSystem_DimensionSet(CoordinateSystem,2,Err)
  !Finish the creation of the coordinate system
  CALL cmfe_CoordinateSystem_CreateFinish(CoordinateSystem,Err)
  

  !Start the creation of the region
  CALL cmfe_Region_Initialise(Region,Err)
  CALL cmfe_Region_CreateStart(RegionUserNumber,WorldRegion,Region,Err)
  !Set the regions coordinate system to the 1D RC coordinate system that we have created
  CALL cmfe_Region_CoordinateSystemSet(Region,CoordinateSystem,Err)
  !Finish the creation of the region
  CALL cmfe_Region_CreateFinish(Region,Err)
  
!_________________________________________________________________________________________________
  !Start the creation of a biilinear-simplex basis

  CALL CreateBasis(Basis,BasisUserNumber,CMFE_BASIS_SIMPLEX_TYPE,CMFE_Basis_Quadratic_Simplex_Interpolation,2,Err)
  CALL CreateBasis(PressureBasis,PressureBasisUserNumber,CMFE_BASIS_SIMPLEX_TYPE,CMFE_Basis_Linear_Simplex_Interpolation,2,Err)
  

  !Time to create a mesh
  !Read in nodes.
  open(unit=10,file=NODEFILE,status='old',action='read',iostat=st)
  IF(st>0)then
    print *,'Error opening node file',st
    STOP
  ELSE
    PRINT *,'Node file opened correctly'
    READ(10,*) NUMBER_OF_NODES, NUMBER_OF_COORDS, NUMBER_OF_ATTRIBUTES, BOUNDARY_MARKER
    ALLOCATE(NodeNums(NUMBER_OF_NODES,2))
    ALLOCATE(NodeCoords(NUMBER_OF_NODES,NUMBER_OF_COORDS))
    DO i = 1,NUMBER_OF_NODES
      READ(10,*) NodeNums(i,1),NodeCoords(i,1),NodeCoords(i,2),NodeNums(i,2)
    ENDDO
  ENDIF
  CLOSE(10)
  PRINT *, 'Total Nodes',NUMBER_OF_NODES
  !Read in elements
  OPEN(unit=11,file=ELEMFILE,status='old',action='read',iostat=st)
  IF(st>0)THEN
    PRINT *,'Error opening element file',st
    STOP
  ELSE
    PRINT *,'Element file opened successfully'
    READ(11,*) NUMBER_OF_ELEMENTS,NODES_PER_ELE,ELE_ATTRIBUTES
    ALLOCATE(ElemMap(NUMBER_OF_ELEMENTS,7))
    DO i = 1,NUMBER_OF_ELEMENTS
      READ(11,*) ElemMap(i,1),ElemMap(i,2),ElemMap(i,3),ElemMap(i,4), &
        & ElemMap(i,5),ElemMap(i,6),ElemMap(i,7)
    ENDDO
  ENDIF 
  CLOSE(11)
  PRINT *,'Total Elements',NUMBER_OF_ELEMENTS
  CALL MPI_BCAST(NUMBER_OF_ELEMENTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,MPI_IERROR)

  CALL cmfe_Nodes_Initialise(Nodes,Err)
  CALL cmfe_Nodes_CreateStart(Region,NUMBER_OF_NODES,Nodes,Err)
  CALL cmfe_Nodes_CreateFinish(Nodes,Err)
  
  CALL cmfe_Mesh_Initialise(Mesh,Err)
  CALL cmfe_Mesh_CreateStart(MeshUserNumber,Region,NUMBER_OF_COORDS,Mesh,Err)
  CALL cmfe_Mesh_NumberOfElementsSet(Mesh,NUMBER_OF_ELEMENTS,Err)
  CALL cmfe_Mesh_NumberOfComponentsSet(Mesh,2,Err)
  
  CALL cmfe_MeshElements_Initialise(MeshElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,1,Basis,MeshElements,Err)
  CALL cmfe_MeshElements_CreateStart(Mesh,2,PressureBasis,PressureElements,Err)

  DO i = 1,NUMBER_OF_ELEMENTS
    element = ElemMap(i,1)
    CALL cmfe_MeshElements_NodesSet(MeshElements,element,(/ElemMap(i,2),ElemMap(i,3), &
      &   ElemMap(i,4),ElemMap(i,7),ElemMap(i,5),ElemMap(i,6)/),Err)
    CALL cmfe_MeshElements_NodesSet(PressureElements,element, &
      & (/ElemMap(i,2),ElemMap(i,3),ElemMap(i,4)/),Err)
  ENDDO
  CALL cmfe_MeshElements_CreateFinish(MeshElements,Err)
  CALL cmfe_MeshElements_CreateFinish(PressureElements,Err)
  CALL cmfe_Mesh_CreateFinish(Mesh,Err)

  !Create a decomposition
  CALL cmfe_Decomposition_Initialise(Decomposition,Err)
  CALL cmfe_Decomposition_CreateStart(DecompositionUserNumber,Mesh,Decomposition,Err)
  !Set the decomposition to be a general decomposition with the specified number of domains
  CALL cmfe_Decomposition_TypeSet(Decomposition,cmfe_DECOMPOSITION_CALCULATED_TYPE,Err)
  CALL cmfe_Decomposition_NumberOfDomainsSet(Decomposition,NumberOfComputationalNodes,Err)
  !Finish the decomposition
  CALL cmfe_Decomposition_CreateFinish(Decomposition,Err)

  !Start to create a default (geometric) field on the region
  CALL cmfe_Field_Initialise(GeometricField,Err)
  CALL cmfe_Field_CreateStart(GeometricFieldUserNumber,Region,GeometricField,Err)
  !Set the decomposition to use
  CALL cmfe_Field_MeshDecompositionSet(GeometricField,Decomposition,Err)
  !Set the domain to be used by the field components. We have 2 field components in 1 mesh component
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,2,1,Err)

  !Finish creating the field
  CALL cmfe_Field_CreateFinish(GeometricField,Err)

  !Set the geometric field values

  DO i = 1,NUMBER_OF_NODES
    node = NodeNums(i,1)
    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,node,1,NodeDomain,Err)
    IF(NodeDomain==ComputationalNodeNumber) THEN
      nodex = NodeCoords(i,1)
      nodey = NodeCoords(i,2)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,1,nodex,Err)
      CALL cmfe_Field_ParameterSetUpdateNode(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,1, &
       &   1,node,2,nodey,Err)
     ENDIF
    ENDDO
  CALL cmfe_Field_ParameterSetUpdateStart(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)
  CALL cmfe_Field_ParameterSetUpdateFinish(GeometricField,cmfe_FIELD_U_VARIABLE_TYPE,cmfe_FIELD_VALUES_SET_TYPE,Err)

  IF(EXPORT_FIELD) THEN
    CALL cmfe_Fields_Initialise(Fields,Err)
    CALL cmfe_Fields_Create(Region,Fields,Err)
    CALL cmfe_Fields_NodesExport(Fields,"Cell_Geom","FORTRAN",Err)
    CALL cmfe_Fields_ElementsExport(Fields,"Cell_Geom","FORTRAN",Err)
    CALL cmfe_Fields_Finalise(Fields,Err)
  ENDIF 

!__________________________________________________________________________________________________________
 !SET UP MECHANICS 

!__________________________________________________________________________________________________________

  !Create the mechannics equations_set
  CALL cmfe_Field_Initialise(FibreField,Err)
  CALL cmfe_Field_CreateStart(FibreFieldUserNumber,Region,FibreField,Err)
  CALL cmfe_Field_TypeSet(FibreField,CMFE_FIELD_FIBRE_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(FibreField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(FibreField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(FibreField,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL cmfe_Field_VariableLabelSet(FibreField,CMFE_FIELD_U_VARIABLE_TYPE,"Fibre",Err)
  CALL cmfe_Field_CreateFinish(FibreField,Err)

  !Create the stiffness field
  CALL cmfe_Field_Initialise(MechStiffnessField,Err)
  CALL cmfe_Field_CreateStart(MechStiffnessFieldUserNumber,Region,MechStiffnessField,Err)
  CALL cmfe_Field_TypeSet(MechStiffnessField,CMFE_FIELD_MATERIAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(MechStiffnessField,Decomposition,Err)        
  CALL cmfe_Field_GeometricFieldSet(MechStiffnessField,GeometricField,Err)
  CALL cmfe_Field_NumberOfVariablesSet(MechStiffnessField,1,Err)
  CALL cmfe_Field_NumberOfComponentsSet(MechStiffnessField,CMFE_FIELD_U_VARIABLE_TYPE,2,Err)  
  CALL cmfe_Field_ComponentMeshComponentSet(MechStiffnessField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(MechStiffnessField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL cmfe_Field_VariableLabelSet(MechStiffnessField,CMFE_FIELD_U_VARIABLE_TYPE,"Stiffness",Err)
  CALL cmfe_Field_CreateFinish(MechStiffnessField,Err)
  
  !Set Mooney-Rivlin constants c10 and c01 to 2.0 and 6.0 respectively.
  CALL cmfe_Field_ComponentValuesInitialise(MechStiffnessField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,1,1.0_CMISSRP,Err)
  CALL cmfe_Field_ComponentValuesInitialise(MechStiffnessField,CMFE_FIELD_U_VARIABLE_TYPE, &
    & CMFE_FIELD_VALUES_SET_TYPE,2,0.0_CMISSRP,Err)
  !CALL cmfe_Field_ComponentValuesInitialise(MechStiffnessField,CMFE_FIELD_V_VARIABLE_TYPE, &
  !  & CMFE_FIELD_VALUES_SET_TYPE,1,0.0_CMISSRP,Err)

  CALL cmfe_Field_Initialise(DeformedField,Err)
  CALL cmfe_Field_CreateStart(DeformedFieldUserNumber,Region,DeformedField,Err)
  CALL cmfe_Field_TypeSet(DeformedField,CMFE_FIELD_GEOMETRIC_GENERAL_TYPE,Err)
  CALL cmfe_Field_MeshDecompositionSet(DeformedField,Decomposition,Err)
  CALL cmfe_Field_GeometricFieldSet(DeformedField,GeometricField,Err)
  CALL cmfe_Field_DependentTypeSet(DeformedField,CMFE_FIELD_DEPENDENT_TYPE,Err)
  CALL cmfe_Field_NumberOfVariablesSet(DeformedField,2,Err)
  CALL cmfe_Field_VariableTypesSet(DeformedField,[CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_DELUDELN_VARIABLE_TYPE],Err)
  CALL cmfe_Field_VariableLabelSet(DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,"Deformation",Err)

  CALL cmfe_Field_NumberOfComponentsSet(DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,3,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,2,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,3,2,Err)

  CALL cmfe_Field_NumberOfComponentsSet(DeformedField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DeformedField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,1,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DeformedField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,2,1,Err)
  CALL cmfe_Field_ComponentMeshComponentSet(DeformedField,CMFE_FIELD_DELUDELN_VARIABLE_TYPE,3,2,Err)

  CALL cmfe_Field_CreateFinish(DeformedField,Err)


  CALL cmfe_Field_Initialise(MechEquationsSetField,Err)
  CALL cmfe_EquationsSet_CreateStart(MechEquationSetUserNumber,Region,FibreField,[CMFE_EQUATIONS_SET_ELASTICITY_CLASS, &
    & CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE,CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE],MechEquationsSetFieldUserNumber, &
    & MechEquationsSetField,MechEquationsSet,Err)
  CALL cmfe_EquationsSet_CreateFinish(MechEquationsSet,Err)
  CALL cmfe_EquationsSet_DependentCreateStart(MechEquationsSet,DeformedFieldUserNumber,DeformedField,Err)
  CALL cmfe_EquationsSet_DependentCreateFinish(MechEquationsSet,Err)
  CALL cmfe_EquationsSet_MaterialsCreateStart(MechEquationsSet,MechStiffnessFieldUserNumber,MechStiffnessField,Err)  
  CALL cmfe_EquationsSet_MaterialsCreateFinish(MechEquationsSet,Err)
  


  !Create the equations set equations
  CALL cmfe_Equations_Initialise(MechEquations,Err)
  CALL cmfe_EquationsSet_EquationsCreateStart(MechEquationsSet,MechEquations,Err)
  CALL cmfe_Equations_SparsityTypeSet(MechEquations,CMFE_EQUATIONS_SPARSE_MATRICES,Err)
  CALL cmfe_Equations_OutputTypeSet(MechEquations,CMFE_EQUATIONS_NO_OUTPUT,Err)
  CALL cmfe_EquationsSet_EquationsCreateFinish(MechEquationsSet,Err)

  !Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 1,DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,1,Err)
  CALL cmfe_Field_ParametersToFieldParametersComponentCopy(GeometricField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE, &
    & 2,DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,2,Err)
  CALL cmfe_Field_ComponentValuesInitialise(DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,CMFE_FIELD_VALUES_SET_TYPE,3,0.0_CMISSRP, &
    & Err)
    
  !Define the problem
  CALL cmfe_Problem_Initialise(MechProblem,Err)
  CALL cmfe_Problem_CreateStart(MechProblemUserNumber,[CMFE_PROBLEM_ELASTICITY_CLASS,CMFE_PROBLEM_FINITE_ELASTICITY_TYPE, &
    & CMFE_PROBLEM_NO_SUBTYPE],MechProblem,Err)
  CALL cmfe_Problem_CreateFinish(MechProblem,Err)

  !Create the problem control loop
  CALL cmfe_Problem_ControlLoopCreateStart(MechProblem,Err)
  CALL cmfe_ControlLoop_Initialise(LoadLoop,Err)
  CALL cmfe_Problem_ControlLoopGet(MechProblem,CMFE_CONTROL_LOOP_NODE,LoadLoop,Err)
  CALL cmfe_ControlLoop_TypeSet(LoadLoop,CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE,Err)
  CALL cmfe_ControlLoop_MaximumIterationsSet(LoadLoop,6,Err)
  CALL cmfe_Problem_ControlLoopCreateFinish(MechProblem,Err)

  !Create the problem solvers
  CALL cmfe_Solver_Initialise(MechSolver,Err)
  CALL cmfe_Solver_Initialise(MechLinearSolver,Err)
  CALL cmfe_Problem_SolversCreateStart(MechProblem,Err)
  CALL cmfe_Problem_SolverGet(MechProblem,CMFE_CONTROL_LOOP_NODE,1,MechSolver,Err)
  CALL cmfe_Solver_OutputTypeSet(MechSolver,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(MechSolver,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(MechSolver,MechLinearSolver,Err)
  CALL cmfe_Solver_LinearTypeSet(MechLinearSolver,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
  CALL cmfe_Problem_SolversCreateFinish(MechProblem,Err)

  !Create the problem solver equations
  CALL cmfe_Solver_Initialise(MechSolver,Err)
  CALL cmfe_SolverEquations_Initialise(MechSolverEquations,Err)
  CALL cmfe_Problem_SolverEquationsCreateStart(MechProblem,Err)
  CALL cmfe_Problem_SolverGet(MechProblem,CMFE_CONTROL_LOOP_NODE,1,MechSolver,Err)
  CALL cmfe_Solver_SolverEquationsGet(MechSolver,MechSolverEquations,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(MechSolverEquations,MechEquationsSet,MechEquationsSetIndex,Err)
  CALL cmfe_Problem_SolverEquationsCreateFinish(MechProblem,Err)


  !Prescribe boundary conditions (absolute nodal parameters)
  !pure stretching
!  MechBCNODES = (/1,2,3,98,195,236,883,912,982,1019,1054,1056/)
!  MechBCNODES_PULL=(/49,50,51,52,196,255,428,863,888,1057,1143,1147/)
! Micropipette version 0
  MechBCNODES_PULL = (/1,2,3,98,195,236,883,912,982,1019,1054,1056/)
  MechBCNODES=(/4,5,96,97,345,580/)
!  MechBCNODES_PULL = (/1,195,236,982,1019,1056/)
!  MechBCNODES=(/4,5,96,97,345,580/)

  CALL cmfe_BoundaryConditions_Initialise(MechBCs,Err)
  CALL cmfe_SolverEquations_BoundaryConditionsCreateStart(MechSolverEquations,MechBCs,Err)

  DO node=1,12
    NODE_NUMBER=MechBCNODES_PULL(node)
    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,-0.15_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,2, &
      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
  ENDDO    

  DO node=1,6
    NODE_NUMBER=MechBCNODES(node)
    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,1, &
      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,2, &
      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
  ENDDO    
  
!    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,1,2, &
!      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
!    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,50,2, &
!      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
!
!    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,25,1, &
!      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
!
!    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,75,1, &
!      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
  
!  DO node=1,12
!    NODE_NUMBER=MechBCNODES(node)
!    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,1, &
!      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
!    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,2, &
!      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
!    NODE_NUMBER=MechBCNODES_PULL(node)
!    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,1, &
!      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.5_CMISSRP,Err)
!    CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!      & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,2, &
!      & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
!  ENDDO    
    
!  DO node = 1,NUMBER_OF_NODES
!    INPIPETTE = .FALSE.
!    NODE_NUMBER = NodeNums(node,1)
!    
!    CALL cmfe_Decomposition_NodeDomainGet(Decomposition,NODE_NUMBER,1,NodeDomain,Err)
!    IF(NodeDomain==ComputationalNodeNumber) THEN  
!      IF(NodeNums(node,2).EQ.1) THEN
!        DO i=1,2
!          IF(NODE_NUMBER.EQ.MechBCNODES(i)) THEN
!            INPIPETTE = .TRUE.
!          ENDIF
!        ENDDO 
!        IF(INPIPETTE) THEN
!          CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField, &
!            & CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,1, &
!            & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,-0.0_CMISSRP,Err)
!        
!        ELSE
!          CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,1, &
!            & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
!
!          CALL cmfe_BoundaryConditions_AddNode(MechBCs,DeformedField,CMFE_FIELD_U_VARIABLE_TYPE,1,1,NODE_NUMBER,2, &
!            & CMFE_BOUNDARY_CONDITION_FIXED_INCREMENTED,0.0_CMISSRP,Err)
!        ENDIF
!      ENDIF
!    ENDIF    
!  ENDDO
  CALL cmfe_SolverEquations_BoundaryConditionsCreateFinish(MechSolverEquations,Err)

!__________________________________________________________________________________________________________

  !Solve problem
  CALL cmfe_Problem_Solve(MechProblem,Err)
  !Output solution
  CALL cmfe_Fields_Initialise(Fields,Err)
  CALL cmfe_Fields_Create(Region,Fields,Err)
  CALL cmfe_Fields_NodesExport(Fields,"Micropipette","FORTRAN",Err)
  CALL cmfe_Fields_ElementsExport(Fields,"Micropipette","FORTRAN",Err)
  CALL cmfe_Fields_Finalise(Fields,Err)

  CALL cmfe_Finalise(Err)
  
  WRITE(*,'(A)') "Program successfully completed."
  

  STOP
  
END PROGRAM ActinWaves
