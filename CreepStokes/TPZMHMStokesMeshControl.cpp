//
//  TPZMHMixedMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMStokesMeshControl.h"
#include "TPZMHMeshControl.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"
#include "TPZMatLaplacian.h"
#include "TPZLagrangeMultiplier.h"

#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"

TPZMHMStokesMeshControl::TPZMHMStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices) : TPZMHMixedMeshControl(gmesh,coarseindices), fBCTractionMatIds()
{
    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCoarseAveragePressMesh = new TPZCompMesh(fGMesh);
    fCoarseDistrFluxMesh = new TPZCompMesh(fGMesh);
    fsetCoarseAverageMultipliers = false;
    
    fBCTractionMatIds.clear();
    for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
    {
        fBCTractionMatIds[*it]=*it-10;
    }
    
}

TPZMHMStokesMeshControl::TPZMHMStokesMeshControl(int dimension) : TPZMHMixedMeshControl(dimension), fBCTractionMatIds(){

    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCoarseAveragePressMesh = new TPZCompMesh(fGMesh);
    fCoarseDistrFluxMesh = new TPZCompMesh(fGMesh);
    fsetCoarseAverageMultipliers = false;
    
    fBCTractionMatIds.clear();
    for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
    {
        fBCTractionMatIds[*it]=*it-10;
    }
}

TPZMHMStokesMeshControl::TPZMHMStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMixedMeshControl(gmesh), fBCTractionMatIds(){

    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCoarseAveragePressMesh = new TPZCompMesh(fGMesh);
    fCoarseDistrFluxMesh = new TPZCompMesh(fGMesh);
    fsetCoarseAverageMultipliers = false;
    
    fBCTractionMatIds.clear();
    for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
    {
        fBCTractionMatIds[*it]=*it-10;
    }
}


TPZMHMStokesMeshControl::TPZMHMStokesMeshControl(const TPZMHMStokesMeshControl &copy) : TPZMHMixedMeshControl(copy){
    
    this->operator=(copy);
}

TPZMHMStokesMeshControl &TPZMHMStokesMeshControl::operator=(const TPZMHMStokesMeshControl &cp){

    fDistrFluxMesh = cp.fDistrFluxMesh;
    fAveragePressMesh = cp.fAveragePressMesh;
    fCoarseAveragePressMesh = cp.fCoarseAveragePressMesh;
    fCoarseDistrFluxMesh = cp.fCoarseDistrFluxMesh;
    fBCTractionMatIds = cp.fBCTractionMatIds;
    fsetCoarseAverageMultipliers = cp.fsetCoarseAverageMultipliers;
    return *this;
}


TPZMHMStokesMeshControl::~TPZMHMStokesMeshControl()
{

}

void TPZMHMStokesMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    
    //InsertPeriferalMaterialObjects();  //Skeleton com dirichlet v=0
    CreateHDivMHMMesh();
    
    InsertBCSkeleton();
    InsertInternalSkeleton();

#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    
    InsertPeriferalPressureMaterialObjects();
    CreatePressureAndTractionMHMMesh();
    
    InsertDistributedFluxMaterialObjects();
    CreateDistributedFluxMHMMesh();
   
    InsertPeriferalAveragePressMaterialObjects();
    CreateAveragePressMHMMesh();

    if (fsetCoarseAverageMultipliers) {
        CreateLagrangeMultiplierMesh();
    }

  //  CreateCoarseAveragePressMHMMesh();
  //  CreateCoarseDistributedFluxMHMMesh();
    
    CreateMultiPhysicsMHMMesh();
    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();
    
    fCMesh->ComputeNodElCon();
#ifdef PZDEBUG
    {
        std::ofstream out("Friendly.txt");
        PrintFriendly(out);
    }
    CheckMeshConsistency();
#endif
    
    //if (usersubstructure) {
    //HideTheElements();
    //GroupAndCondense(fCMesh.operator->());
    //}
    fNumeq = fCMesh->NEquations();
    
}

void TPZMHMStokesMeshControl::CreatePressureAndTractionMHMMesh(){
    
    CreatePressureMHMMesh();
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    
    int64_t nskeletonconnects = fPressureFineMesh->NConnects();
    int porder = fpOrderInternal;
    TPZCompMesh * cmeshTraction = fPressureFineMesh.operator->();
    gmesh->ResetReference();
    cmeshTraction->SetName("PressureAndTractionMesh");
    cmeshTraction->SetDimModel(gmesh->Dimension()-1);
    cmeshTraction->SetAllCreateFunctionsDiscontinuous();
    cmeshTraction->SetDefaultOrder(porder-1);
    int meshdim = cmeshTraction->Dimension();
    
    std::set<int> matids;
    TPZMaterial *mat = cmeshTraction->FindMaterial(fTractionMatId);
    if (mat && mat->Dimension() == meshdim) {
        matids.insert(fTractionMatId);
    }
    
    for (auto it:fMaterialBCIds) {
        int dsmatid = fBCTractionMatIds[it];
        TPZMaterial *mat = cmeshTraction->FindMaterial(fBCTractionMatIds[it]);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(fBCTractionMatIds[it]);
        }
    }
    
    cmeshTraction->AutoBuild(matids);
    fPressureFineMesh->ExpandSolution();
    
//    if(0)
//    {
        std::ofstream out("PressureAndTractionFineMesh.txt");
        fPressureFineMesh->Print(out);
//    }
    

#ifdef PZDEBUG
    {
        int64_t nel = fGMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel && gel->MaterialId() == 1) {
                if (gel->Dimension() != fGMesh->Dimension()) {
                    DebugStop();
                }
            }
        }
    }
#endif
    
    int64_t nc = cmeshTraction->NConnects();
    if(nskeletonconnects != 0){
  //      DebugStop();
    }
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshTraction->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    int64_t nel = cmeshTraction->NElements();
//    for (int64_t el=0; el<nel; el++)
//    {
//        TPZCompEl *cel = cmeshTraction->Element(el);
//#ifdef PZDEBUG
//        if (! cel) {
//            DebugStop();
//        }
//#endif
//        TPZGeoEl *gel = cel->Reference();
//        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
//        {
//            continue;
//        }
//#ifdef PZDEBUG
//        if (fGeoToMHMDomain[gel->Index()] == -1) {
//            DebugStop();
//        }
//#endif
//
//        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
//    }

}

void TPZMHMStokesMeshControl::InsertInternalSkeleton(){
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (gel->Dimension()!=fGMesh->Dimension()-1) {
            continue;
        }
        if (gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        if (gel->MaterialId()==fInternalWrapMatId) {
            TPZGeoElBC(gelside, fTractionMatId);
        }
        if (gel->MaterialId()==fSkeletonMatId) {
            TPZGeoElBC(gelside, fTractionMatId);
        }
        
    }

}

void TPZMHMStokesMeshControl::InsertBCSkeleton(){
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        
        for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
        {
            if (gel->MaterialId()==*it) {
                TPZGeoElBC(gelside, fBCTractionMatIds[*it]);
            }
        }
    }
    
}


void TPZMHMStokesMeshControl::InsertPeriferalPressureMaterialObjects()
{
    TPZCompMesh * cmeshPressure = fPressureFineMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshPressure->MaterialVec().find(matid) == cmeshPressure->MaterialVec().end())
        {
            TPZNullMaterial *matl2 = new TPZNullMaterial((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshPressure->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
    // Material for interior traction:
    
    TPZVecL2 *matTraction = new TPZVecL2(fTractionMatId);
    matTraction->SetDimension(fGMesh->Dimension()-1);
    cmeshPressure->InsertMaterialObject(matTraction);

    for (auto it:fMaterialBCIds)
    {
        if (fBCTractionMatIds.size()!=fMaterialBCIds.size()) {
            DebugStop();
        }
        
        int matid= fBCTractionMatIds[it];
        if (cmeshPressure->MaterialVec().find(matid) == cmeshPressure->MaterialVec().end())
        {
            TPZVecL2 *matBCTraction = new TPZVecL2(matid);
            matBCTraction->SetDimension(fGMesh->Dimension()-1);
            cmeshPressure->InsertMaterialObject(matBCTraction);
        }
        
    }
    
}

void TPZMHMStokesMeshControl::InsertPeriferalAveragePressMaterialObjects(){
    
    TPZCompMesh *cmeshAverPressure = fAveragePressMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshAverPressure->MaterialVec().find(matid) == cmeshAverPressure->MaterialVec().end())
        {
            TPZVecL2 *matl2 = new TPZVecL2((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshAverPressure->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
    
}

void TPZMHMStokesMeshControl::CreateAveragePressMHMMesh(){
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshAverPressute = fAveragePressMesh.operator->();
    gmesh->ResetReference();
    cmeshAverPressute->SetName("AveragePressureMesh");
    cmeshAverPressute->SetDimModel(gmesh->Dimension());
    cmeshAverPressute->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshAverPressute->SetDefaultOrder(porder);
    
    int meshdim = cmeshAverPressute->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshAverPressute->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshAverPressute->AutoBuild(matids);
    fAveragePressMesh->ExpandSolution();
    
    if(1)
    {
        std::ofstream out("AveragePressureMesh.txt");
        fAveragePressMesh->Print(out);
    }
    
    
    int64_t nel = fAveragePressMesh->NElements();
    for(int64_t i=0; i<nel; i++){
        TPZCompEl *cel = cmeshAverPressute->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) DebugStop();
        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
        {
            DebugStop();
        }
        celdisc->SetTotalOrderShape();
        celdisc->SetFalseUseQsiEta();
    }
    
    int64_t nc = cmeshAverPressute->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshAverPressute->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshAverPressute->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            continue;
        }
#ifdef PZDEBUG
        if (fGeoToMHMDomain[gel->Index()] == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
    }
    
    
    return;
    
}

void TPZMHMStokesMeshControl::InsertPeriferalCoarseAveragePressMaterialObjects(){
    
    TPZCompMesh *cmeshCoarseAverPressure = fCoarseAveragePressMesh.operator->();
    
    TPZVecL2 *matl2 = new TPZVecL2(fSkeletonMatId);
    matl2->SetNStateVariables(fNState);
    matl2->SetDimension(fGMesh->Dimension());
    cmeshCoarseAverPressure->InsertMaterialObject(matl2);
    
}

void TPZMHMStokesMeshControl::CreateCoarseAveragePressMHMMesh(){
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshAverPressute = fCMeshConstantPressure.operator->();
    gmesh->ResetReference();
    cmeshAverPressute->SetName("CoarseAveragePressureMesh");
    cmeshAverPressute->SetDimModel(gmesh->Dimension());
    cmeshAverPressute->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshAverPressute->SetDefaultOrder(porder);
    
    int meshdim = cmeshAverPressute->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshAverPressute->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshAverPressute->AutoBuild(matids);
    fAveragePressMesh->ExpandSolution();
    
    if(1)
    {
        std::ofstream out("CoarseAveragePressureMesh.txt");
        fAveragePressMesh->Print(out);
    }
    
    
//    int64_t nel = fAveragePressMesh->NElements();
//    for(int64_t i=0; i<nel; i++){
//        TPZCompEl *cel = cmeshAverPressute->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        if(!celdisc) DebugStop();
//        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
//        {
//            DebugStop();
//        }
//        celdisc->SetTotalOrderShape();
//        celdisc->SetFalseUseQsiEta();
//    }
    
    int64_t nc = cmeshAverPressute->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshAverPressute->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    
//    for (int64_t el=0; el<nel; el++)
//    {
//        TPZCompEl *cel = cmeshAverPressute->Element(el);
//#ifdef PZDEBUG
//        if (! cel) {
//            DebugStop();
//        }
//#endif
//        TPZGeoEl *gel = cel->Reference();
//        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
//        {
//            continue;
//        }
//#ifdef PZDEBUG
//        if (fGeoToMHMDomain[gel->Index()] == -1) {
//            DebugStop();
//        }
//#endif
//
//        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
//    }
    
    
    return;
    
}



void TPZMHMStokesMeshControl::InsertDistributedFluxMaterialObjects(){
    
    TPZCompMesh *cmeshDistributedFlux = fDistrFluxMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshDistributedFlux->MaterialVec().find(matid) == cmeshDistributedFlux->MaterialVec().end())
        {
            TPZVecL2 *matl2 = new TPZVecL2((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshDistributedFlux->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
}


void TPZMHMStokesMeshControl::CreateDistributedFluxMHMMesh(){

    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshDistributedFlux = fDistrFluxMesh.operator->();
    gmesh->ResetReference();
    cmeshDistributedFlux->SetName("DistributedFluxMesh");
    cmeshDistributedFlux->SetDimModel(gmesh->Dimension());
    cmeshDistributedFlux->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshDistributedFlux->SetDefaultOrder(porder);
    
    int meshdim = cmeshDistributedFlux->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshDistributedFlux->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshDistributedFlux->AutoBuild(matids);
    fDistrFluxMesh->ExpandSolution();
    
    if(1)
    {
        std::ofstream out("DistributedFluxMesh.txt");
        fDistrFluxMesh->Print(out);
    }
    
    
    int64_t nel = fDistrFluxMesh->NElements();
    for(int64_t i=0; i<nel; i++){
        TPZCompEl *cel = cmeshDistributedFlux->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) DebugStop();
        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
        {
            DebugStop();
        }
        celdisc->SetTotalOrderShape();
        celdisc->SetFalseUseQsiEta();
    }
    
    int64_t nc = cmeshDistributedFlux->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshDistributedFlux->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshDistributedFlux->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            continue;
        }
#ifdef PZDEBUG
        if (fGeoToMHMDomain[gel->Index()] == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
    }
    
    
    return;
    
}

void TPZMHMStokesMeshControl::InsertCoarseDistributedFluxMaterialObjects(){
    
    TPZCompMesh *cmeshDistributedFlux = fDistrFluxMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshDistributedFlux->MaterialVec().find(matid) == cmeshDistributedFlux->MaterialVec().end())
        {
            TPZVecL2 *matl2 = new TPZVecL2((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshDistributedFlux->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
}

void TPZMHMStokesMeshControl::CreateCoarseDistributedFluxMHMMesh(){
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshDistributedFlux = fCMeshLagrange.operator->();
    gmesh->ResetReference();
    cmeshDistributedFlux->SetName("CoarseDistributedFluxMesh");
    cmeshDistributedFlux->SetDimModel(gmesh->Dimension());
    cmeshDistributedFlux->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshDistributedFlux->SetDefaultOrder(porder);
    
    int meshdim = cmeshDistributedFlux->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshDistributedFlux->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshDistributedFlux->AutoBuild(matids);
    fDistrFluxMesh->ExpandSolution();
    
    if(1)
    {
        std::ofstream out("CoarseDistributedFluxMesh.txt");
        fDistrFluxMesh->Print(out);
    }
    
    
//    int64_t nel = fDistrFluxMesh->NElements();
//    for(int64_t i=0; i<nel; i++){
//        TPZCompEl *cel = cmeshDistributedFlux->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        if(!celdisc) DebugStop();
//        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
//        {
//            DebugStop();
//        }
//        celdisc->SetTotalOrderShape();
//        celdisc->SetFalseUseQsiEta();
//    }
    
    int64_t nc = cmeshDistributedFlux->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshDistributedFlux->ConnectVec()[ic].SetLagrangeMultiplier(1);
    }
    gmesh->ResetReference();
    
//    for (int64_t el=0; el<nel; el++)
//    {
//        TPZCompEl *cel = cmeshDistributedFlux->Element(el);
//#ifdef PZDEBUG
//        if (! cel) {
//            DebugStop();
//        }
//#endif
//        TPZGeoEl *gel = cel->Reference();
//        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
//        {
//            continue;
//        }
//#ifdef PZDEBUG
//        if (fGeoToMHMDomain[gel->Index()] == -1) {
//            DebugStop();
//        }
//#endif
//
//        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
//    }
    
    return;
    
}


void TPZMHMStokesMeshControl::CreateMultiPhysicsMHMMesh()
{
    TPZManVector<TPZCompMesh *,6 > cmeshes(4);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    cmeshes[2] = fDistrFluxMesh.operator->();
    cmeshes[3] = fAveragePressMesh.operator->();
    if (fsetCoarseAverageMultipliers) {
        cmeshes.Resize(6);
        cmeshes[4] = fCMeshLagrange.operator->();
        cmeshes[5] = fCMeshConstantPressure.operator->();
    }
    
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    
    // Malha computacional
    TPZCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    
    BuildMultiPhysicsMesh();
    TPZManVector<TPZCompMesh * ,3> meshvector;
    
    if(0)
    {
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
        std::ofstream out3("HDivMesh.txt");
        fFluxMesh->Print(out3);
        std::ofstream out4("PressureMesh.txt");
        fPressureFineMesh->Print(out4);
    }
    
    meshvector = cmeshes;
    
    JoinSubdomains(meshvector, MixedFluxPressureCmesh);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);
    
    if(1)
    {
        std::ofstream file("multiphysics.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(MixedFluxPressureCmesh, file,true);
        std::ofstream out("multiphysics_before_condense.txt");
        MixedFluxPressureCmesh->Print(out);
    }

    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
    //CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-1);
    //CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-2);
    
    CreateMultiPhysicsInterfaceElements();
    CreateMultiPhysicsBCInterfaceElements();

    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
    GroupAndCondense(MixedFluxPressureCmesh);

#ifdef PZDEBUG
    if(1)
    {
        std::ofstream out("multiphysics.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    
    return;
    
}

void TPZMHMStokesMeshControl::CreateMultiPhysicsInterfaceElements(){
    
    TPZCompMesh *cmesh = fCMesh.operator->();
    TPZVec<int> m_interfaceVector_ids(2,0);
    m_interfaceVector_ids[0] = fLagrangeMatIdLeft;
    m_interfaceVector_ids[1] = fLagrangeMatIdRight;
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        int meshdim = fGMesh->Dimension();
        int matid = gel->MaterialId();
        
        if (matid != fTractionMatId) {
            continue;
        }
        
        if (gel->HasSubElement() == 1) {
            continue;
        }
        
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZCompElSide celside = gelside.Reference();
        
        if (!celside ) {
            DebugStop();
        }
        
        TPZStack<TPZGeoElSide> neighbourset;
       // gelside.AllNeighbours(neighbourset);
        
        //Find elements with the same mesh dimension :
        TPZGeoElSide neigh = gelside.Neighbour();
        while(neigh != gelside)
        {
            if (neigh.Element()->Dimension()!=fGMesh->Dimension()) {
                neigh = neigh.Neighbour();
                continue;
            }
            
            neighbourset.Push(neigh);
            neigh = neigh.Neighbour();
        }
        
        //gelside.LowerLevelCompElementList2(1);
        
        int nneighs = neighbourset.size();
        if (nneighs!=2) {
            DebugStop();
        }
        
        TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
        LeftElIndices[0]=0;
        RightElIndices[0]=1;
        
        for(int stack_i=0; stack_i <nneighs; stack_i++){
            TPZGeoElSide neigh = neighbourset[stack_i];
            
            if (neigh.Element()->Dimension()!=meshdim) {
                continue;
            }
            
            TPZCompElSide celneigh = neigh.Reference();

            
            if(neigh.Element()->HasSubElement()) {

                TPZStack<TPZGeoElSide> subel;
                neigh.GetAllSiblings(subel);
                
                for (int i_sub =0; i_sub<subel.size(); i_sub++) {
                    
                    TPZCompElSide cel_sub_neigh = subel[i_sub].Reference();
                    
                    TPZGeoElBC gbc_sub(subel[i_sub],m_interfaceVector_ids[stack_i]);
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc_sub.CreatedElement(),index,cel_sub_neigh,celside);
                    elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                    
#ifdef PZDEBUG
                    std::cout << "****Created an interface element between volumetric element " << subel[i_sub].Element()->Index() <<
                    " side " << subel[i_sub].Side() <<
                    " and Skeleton element " << gelside.Element()->Index() << std::endl;
#endif

                }
                
                
            }else{
                
                TPZGeoElBC gbc(gelside,m_interfaceVector_ids[stack_i]);
                int64_t index;
                
                TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celneigh,celside);
                elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);

#ifdef PZDEBUG
                std::cout << "Created an interface element between volumetric element " << neigh.Element()->Index() <<
                " side " << neigh.Side() <<
                " and interior 1D element " << gelside.Element()->Index() << std::endl;

#endif
                
            }
            
        }
        
    }
    
    
}

void TPZMHMStokesMeshControl::CreateMultiPhysicsBCInterfaceElements(){
    
    int matBCinterface = fLagrangeMatIdLeft;
    for (auto it : fMaterialBCIds) {
        
        int matfrom = fBCTractionMatIds[it];
        TPZCompMesh *cmesh = fCMesh.operator->();
        
        int64_t nel = fGMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            int meshdim = fGMesh->Dimension();
            int matid = gel->MaterialId();
            
            if (matid != matfrom) {
                continue;
            }
            
            int nsides = gel->NSides();
            TPZGeoElSide gelside(gel,nsides-1);
            TPZCompElSide celside = gelside.Reference();
            
            TPZStack<TPZGeoElSide> neighbourset;
            gelside.AllNeighbours(neighbourset);
            
            int nneighs = neighbourset.size();
//            if(nneighs!=2){
//                //    DebugStop();
//            }
            
            TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
            LeftElIndices[0]=0;
            RightElIndices[0]=1;
            
            for(int stack_i=0; stack_i <nneighs; stack_i++){
                TPZGeoElSide neigh = neighbourset[stack_i];
                if (neigh.Element()->Dimension()!=meshdim) {
                    continue;
                }
                
                TPZCompElSide celneigh = neigh.Reference();
                if (!celside || !celneigh) {
                    //    DebugStop();
                }
                int64_t neigh_index = neigh.Element()->Index();
                if (neigh.Element()->Dimension()!=meshdim){
                    continue;
                }
                
                if (neigh.Element()->HasSubElement()) {
                    
                    TPZStack<TPZGeoElSide > subelside;
                    neigh.GetAllSiblings(subelside);
                    
                    for (int i_sub = 0; i_sub<subelside.size(); i_sub++) {
                        
                        TPZCompElSide cel_sub_neigh = subelside[i_sub].Reference();
                        
                        TPZGeoElBC gbc_sub(subelside[i_sub],matBCinterface);
                        int64_t index;
                        
                        TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc_sub.CreatedElement(),index,cel_sub_neigh,celside);
                        elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);

#ifdef PZDEBUG
                        std::cout << "****Created an BC interface element between volumetric element " << subelside[i_sub].Element()->Index() <<
                        " side " << subelside[i_sub].Side() <<
                        " and boundary 1D element " << gelside.Element()->Index() << std::endl;
#endif
                        
                    }
                    
                }else{
                    
                    TPZGeoElBC gbc(gelside,matBCinterface);
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celneigh,celside);
                    elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                    
#ifdef PZDEBUG
                    std::cout << "Created an BC interface element between volumetric element " << neigh.Element()->Index() <<
                    " side " << neigh.Side() <<
                    " and boundary 1D element " << gelside.Element()->Index() << std::endl;
#endif
                }
                
            }
        }
        
    }
    
    
}

void TPZMHMStokesMeshControl::GroupAndCondense(TPZCompMesh *cmesh_m){
    
    //Criando agrupamento de elementos
    
    int64_t ncompel = cmesh_m->ElementVec().NElements();
    int dim = cmesh_m->Reference()->Dimension();
    
    std::vector<int64_t> GroupIndex;
    TPZStack<TPZElementGroup *> elgroups;
    int count = 0;
    int64_t index =0;
    
    for(int64_t el = 0; el < ncompel; el++){
        
        TPZCompEl *cel = cmesh_m->Element(el);
        if (cel->Dimension()!=dim) {
            continue;
        }
        //GroupIndex[el] = cel->Index();
        count++;
        GroupIndex.resize(count);
        GroupIndex[count-1]=cel->Index();
        TPZElementGroup *GroupEl = new TPZElementGroup(*cmesh_m,index);
        elgroups.Push(GroupEl);
        elgroups[count-1]->AddElement(cel);
    }
    
    
    //Inserindo as respectivas interfaces e condições de contorno
    
    for(int64_t el = 0; el < ncompel; el++){
        TPZCompEl *cel = cmesh_m->Element(el);
        
        TPZMultiphysicsInterfaceElement *interel = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (interel) {
            TPZCompEl *Leftel = interel->LeftElement();
            
            if (Leftel->Dimension()!=dim) {
                continue;
            }
            int leftindex = Leftel->Index();
            
            for(int64_t iel = 0; iel < GroupIndex.size(); iel++){
                if (leftindex==GroupIndex[iel]) {
                    elgroups[iel]->AddElement(cel);
                }
            }
        }
        
        if (!cel) {
            continue;
        }
        
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension()==dim-1) {
            
            TPZBndCond *elBC = dynamic_cast<TPZBndCond *>(cel->Material());
            if (!elBC) {
                continue;
            }
            
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gel, gel->NSides() - 1);
            
            gelside.EqualLevelCompElementList(celstack, 0, 0);
            
            for (auto &celstackindex : celstack) {
                if (celstackindex.Reference().Element()->Dimension() == dim) {
                    int bcindex = celstackindex.Element()->Index();
                    
                    for(int64_t iel = 0; iel < GroupIndex.size(); iel++){
                        if (bcindex==GroupIndex[iel]) {
                            elgroups[iel]->AddElement(cel);
                        }
                    }
                }
            }
        }
    }
    
    cmesh_m->ComputeNodElCon();
    // create condensed elements
    // increase the NumElConnected of one pressure connects in order to prevent condensation

    int nenvel = elgroups.NElements();
    for (int64_t ienv=0; ienv<nenvel; ienv++) {
        TPZElementGroup *elgr = elgroups[ienv];

        int nc = elgroups[ienv]->GetElGroup()[0]->NConnects();
        elgroups[ienv]->GetElGroup()[0]->Connect(nc-1).IncrementElConnected();
        if (fsetCoarseAverageMultipliers) {
                elgroups[ienv]->GetElGroup()[0]->Connect(nc-2).IncrementElConnected();
                elgroups[ienv]->GetElGroup()[0]->Connect(nc-3).IncrementElConnected();
        }
        new TPZCondensedCompEl(elgr);
    }
    
    
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->ExpandSolution();
    
}
