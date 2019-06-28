//
//  TPZMHMixedMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMStokesMeshControl.h"
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

#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"

TPZMHMStokesMeshControl::TPZMHMStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices) : TPZMHMixedMeshControl(gmesh,coarseindices), fBCTractionMatIds()
{
    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    
    fBCTractionMatIds.clear();
    for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
    {
        fBCTractionMatIds[*it]=*it-10;
    }
    
}

TPZMHMStokesMeshControl::TPZMHMStokesMeshControl(int dimension) : TPZMHMixedMeshControl(dimension), fBCTractionMatIds(){

    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);

    fBCTractionMatIds.clear();
    for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
    {
        fBCTractionMatIds[*it]=*it-10;
    }
}

TPZMHMStokesMeshControl::TPZMHMStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMixedMeshControl(gmesh), fBCTractionMatIds(){

    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);

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
    fBCTractionMatIds = cp.fBCTractionMatIds;
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
    
    InsertPeriferalMaterialObjects();
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
    
    InsertPeriferalAveragePressMaterialObjects();
    CreateAveragePressMHMMesh();
    
    InsertDistributedFluxMaterialObjects();
    CreateDistributedFluxMHMMesh();
    
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
    
    if (usersubstructure) {
        HideTheElements();
    }
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
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshTraction->Element(el);
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
    int porder = fpOrderInternal;
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
    
    if(0)
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
    int porder = fpOrderInternal;
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
    
    if(0)
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


void TPZMHMStokesMeshControl::CreateMultiPhysicsMHMMesh()
{
    TPZManVector<TPZCompMesh *,4 > cmeshes(4);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    cmeshes[2] = fAveragePressMesh.operator->();
    cmeshes[3] = fDistrFluxMesh.operator->();

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
    
#ifdef PZDEBUG
    if(1)
    {
        std::ofstream file("cmeshmphys.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(MixedFluxPressureCmesh, file,true);
        std::ofstream out("cmeshmphys.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif
    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
    CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-1);
    //CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-2);
    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    
    if(1)
    {
        std::ofstream out("multiphysics.txt");
        MixedFluxPressureCmesh->Print(out);
    }
    
    return;
    
}
