/*
 *  TPZCouplingDSMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 11/10/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZCouplingDSMaterial.h"
#include "TPZBndCond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"


TPZCouplingDSMaterial::TPZCouplingDSMaterial() : TBase(){
    fk=1;
}

////////////////////////////////////////////////////////////////////

TPZCouplingDSMaterial::TPZCouplingDSMaterial(int matid, int dimension, STATE viscosity,STATE permeability, STATE theta) : 
TBase(matid),fViscosity(viscosity),fTheta(theta),fDimension(dimension)
{
    fk=permeability;
}

////////////////////////////////////////////////////////////////////

TPZCouplingDSMaterial::TPZCouplingDSMaterial(const TPZCouplingDSMaterial &mat) : 
TBase(mat), fViscosity(mat.fViscosity), fTheta(mat.fTheta),fDimension(mat.fDimension)
{
       fk= mat.fk;
}

////////////////////////////////////////////////////////////////////

TPZCouplingDSMaterial::~TPZCouplingDSMaterial(){
    
    
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavec_left, 
std::map<int, TPZMaterialDataT<STATE>> &datavec_right)
{
    data.fNeedsNormal = true;
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

////////////////////////////////////////////////////////////////////

int TPZCouplingDSMaterial::VariableIndex(const std::string &name) {
    
    if (!strcmp("P", name.c_str()))  return 0;
    if (!strcmp("V", name.c_str()))  return 1;
    if (!strcmp("f", name.c_str()))         return 2;
    if (!strcmp("V_exact", name.c_str()))   return 3;
    if (!strcmp("P_exact", name.c_str()))   return 4;
   
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

////////////////////////////////////////////////////////////////////

int TPZCouplingDSMaterial::NSolutionVariables(int var) {
    
    switch(var) {
        
        case 0:
            return 1; // Pressure, Scalar
        case 1:
            return this->Dimension(); // Velocity, Vector
        case 2:
            return this->Dimension(); // f, Vector
        case 3:
            return this->Dimension(); // V_exact, Vector
        case 4:
            return this->Dimension(); // P_exact, Vector

        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) {
    

    int vindex = this->VIndex();
    int pindex = this->PIndex();
    
    TPZManVector<STATE,3> v_h = datavec[vindex].sol[0];
    REAL p_h = datavec[pindex].sol[0][0];
    TPZFNMatrix<9,STATE> gradu(2,1);
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        
        case 0: //Pressure
        {
            Solout[0] = p_h;
        }
            break;
        
        case 1: //Velocity
        {
            Solout[0] = v_h[0]; // Vx
            Solout[1] = v_h[1]; // Vy
        }
            break;
        case 2: //f
        {
            TPZVec<STATE> f;
            if(this->HasForcingFunction()){
                this->fForcingFunction(datavec[vindex].x, f);
            }
            Solout[0] = f[0]; // fx
            Solout[1] = f[1]; // fy
        }
            break;
        
        case 3: //v_exact
        {
            TPZVec<STATE> sol;
            if(this->HasExactSol()){
                this->fExactSol(datavec[vindex].x, sol, gradu);
            }
            Solout[0] = sol[0]; // vx
            Solout[1] = sol[1]; // vy
        }
            break;
        
        case 4: //p_exact
        {
            TPZVec<STATE> sol;
            if(this->HasExactSol()){
                this->fExactSol(datavec[pindex].x, sol, gradu);
            }
            Solout[0] = sol[2]; // px
        }
            break;
            
    
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::Write(TPZStream &buf, int withclassid) {
    
    TPZMaterial::Write(buf, withclassid);
 
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::Read(TPZStream &buf, void *context) {
    
    TPZMaterial::Read(buf, context);

}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi){
    
   
    TPZFMatrix<REAL> &dphiV = dataV.dphix;
    
    const int dim = this->Dimension();
    
    GradPhi.clear();
    GradPhi.resize(dim);
    
    //for each shape
    for(int shape = 0; shape < dphiV.Rows(); shape++){
        
        TPZFMatrix<REAL> GPhi(dim,dim,0.);
        
        for(int i = 0; i < dim; i++){
            
            for(int j = 0; j < dim; j++){
                
                GPhi(i,j) = dphiV(j,shape);// itapopo H1 ??
            
            }//j
        }//i
        
        GradPhi[shape] = GPhi;
        
    }//shape
    
}

// Contricucao dos elementos internos

void TPZCouplingDSMaterial::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
 
    return;
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavec[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[vindex]);
    }

    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    
    
    TPZManVector<REAL> n = datavec[0].normal;
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    TPZVec<double> f;
    TPZFMatrix<STATE> phiVi(fDimension,1,0.0),phiVj(fDimension,1,0.0);
    
    for(int i = 0; i < nshapeV; i++ )
    {
        int iphi = datavec[vindex].fVecShapeIndex[i].second;
        int ivec = datavec[vindex].fVecShapeIndex[i].first;
        TPZFNMatrix<4> phiVti(1,1,0.);
        for (int e=0; e<fDimension; e++) {
            phiVi(e,0) = phiV(iphi,0)*datavec[vindex].fDeformedDirections(e,ivec);
            
        }

        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            
            TPZFNMatrix<4> phiVtj(1,1,0.);
            for (int e=0; e<fDimension; e++) {
                phiVj(e,0) = phiV(jphi,0)*datavec[vindex].fDeformedDirections(e,jvec);
            }

            STATE val = Inner(phiVi,phiVj);
            ek(i,j) += weight * fViscosity * pow(fk,-1./2.)* val ;
            
        }
        
        
    }


}


void TPZCouplingDSMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
       
    DebugStop();
    
}




////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::ContributeInterface(const TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavecleft, 
std::map<int, TPZMaterialDataT<STATE>> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){


#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space for datavecleft Darcy
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright Stokes
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();

    
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    // V - right
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    // P - right
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    //TPZManVector<REAL,3> &n = data.normal;
    TPZManVector<REAL,3> t(2,0.);
    t[0]=-data.normal[1];
    t[1]=data.normal[0];
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);

    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    TPZVec<double> f, v1, v2;
    TPZFMatrix<STATE> phiV1i(fDimension,1,0.0),phiV1j(fDimension,1,0.0);
    

    
    TPZFMatrix<STATE> phiV2i(fDimension,1,0.0),phiV2j(fDimension,1,0.0);
    
    for(int i2 = 0; i2 < nshapeV2; i2++ )
    {
        
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        TPZFNMatrix<4> phiV2ti(1,1,0.),phiV2ni(1,1,0.);
        for (int e=0; e<fDimension; e++) {
            phiV2i(e,0) = phiV2(iphi2,0)*datavecright[vindex].fDeformedDirections(e,ivec2);
            
        }

        
        phiV2ti(0,0) = phiV2i(0,0)*t[0]+phiV2i(1,0)*t[1];
        phiV2ni(0,0) = phiV2i(0,0)*data.normal[0]+phiV2i(1,0)*data.normal[1];
        // matrix A - velocity * test-funtion velocity

        
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            
            STATE fact = weight * phiP1j(0,0)* phiV2ni(0,0);
            
            
            ek(nshapeV1+nshapeP1+i2,nshapeV1+j1) += weight * fact*0.;
            
            
        }

        
        
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            
            TPZFNMatrix<4> phiV2tj(1,1,0.);
            for (int e=0; e<fDimension; e++) {
                phiV2j(e,0) = phiV2(jphi2,0)*datavecright[vindex].fDeformedDirections(e,jvec2);
            }
            phiV2tj(0,0) = phiV2j(0,0)*t[0]+phiV2j(1,0)*t[1];
            
            
            STATE val = phiV2ti(0,0) * phiV2tj(0,0) ;
            
            
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += weight * fViscosity * val*0.;


        }

    }
    
}



void TPZCouplingDSMaterial::ContributeBCInterface(const TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavec, REAL weight, 
TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc){
   
    DebugStop();
    
}



////////////////////////////////////////////////////////////////////
template <typename TVar2>
TVar2 TPZCouplingDSMaterial::Inner(TPZFMatrix<TVar2> &S, TPZFMatrix<TVar2> &T){
    
    //inner product of two tensors

    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }

    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    TVar2 Val = 0;
    
    for(int i = 0; i < S.Cols(); i++){
        for(int j = 0; j < S.Cols(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}


////////////////////////////////////////////////////////////////////
template <typename TVar2>
TVar2 TPZCouplingDSMaterial::InnerVec(TPZFMatrix<TVar2> &S, TPZFMatrix<TVar2> &T){
    
    //inner product of two vectors
    
    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    TVar2 Val = 0;
    
    for(int j = 0; j < S.Cols(); j++){
        for(int i = 0; i < S.Rows(); i++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}



////////////////////////////////////////////////////////////////////

STATE TPZCouplingDSMaterial::Tr( TPZFMatrix<REAL> &GradU ){
 
#ifdef DEBUG
    if( GradU.Rows() != GradU.Cols() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0.;
    
    for(int i = 0; i < GradU.Rows(); i++){
        Val += GradU(i,i);
    }
    
    return Val;
}


/// transform a H1 data structure to a vector data structure
void TPZCouplingDSMaterial::FillVecShapeIndex(TPZMaterialData &data)
{
    data.fDeformedDirections.Resize(fDimension,fDimension);
    data.fDeformedDirections.Identity();
    data.fVecShapeIndex.Resize(fDimension*data.phi.Rows());
    for (int d=0; d<fDimension; d++) {
        for (int i=0; i<data.phi.Rows(); i++) {
            data.fVecShapeIndex[i*fDimension+d].first = d;
            data.fVecShapeIndex[i*fDimension+d].second = i;
        }
    }
}

void TPZCouplingDSMaterial::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors)
{
    DebugStop();
}
