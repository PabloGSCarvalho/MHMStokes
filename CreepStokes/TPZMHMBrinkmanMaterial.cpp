/*
 *  TPZMHMBrinkmanMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMHMBrinkmanMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"

using namespace std;

void TPZMHMBrinkmanMaterial::FillDataRequirementsInterface(TPZMaterialData &data, TPZVec<TPZMaterialData > &datavec_left, TPZVec<TPZMaterialData > &datavec_right)
{
    int nref_left = datavec_left.size();
    for(int iref = 0; iref<nref_left; iref++){
        datavec_left[iref].SetAllRequirements(false);
        datavec_left[iref].fNeedsNormal = true;
        datavec_right[iref].fNeedsNormal = true;
    }
}

void TPZMHMBrinkmanMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    //DebugStop();
    //return;
    // Verificar que
    // os termos mistos devem estar sem viscosidade!
    
    
    //2 = 1 Vel space + 1 Press space for datavecleft
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    TPZManVector<REAL, 3>  &normalLeft = datavecleft[vindex].normal;
    
    data.fNeedsNormal = true;
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    TPZFNMatrix<3, STATE> normalM(fDimension,1,0.);
    for (int e=0; e<fDimension; e++) {
        normalM(e,0)=normal[e];
    }
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    TPZManVector<REAL, 3> tangent(fDimension,0.);
    TPZFNMatrix<3, STATE> tangentV(fDimension,1,0.);
    for(int i=0; i<fDimension; i++) tangent[i] = data.axes(0,i);
    for(int i=0; i<fDimension; i++) tangentV(i,0) = data.axes(0,i);
    
    int nshapeV , nshapeLambda;
    
    nshapeV = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeLambda = datavecright[pindex].phi.Rows();
    
    for(int i1 = 0; i1 < nshapeV; i1++)
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
        TPZFNMatrix<9, STATE> phiVi(fDimension,1);
        for (int e=0; e<fDimension; e++) {
            phiVi(e,0)=datavecleft[vindex].fNormalVec(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
        }
        
        TPZFNMatrix<3, STATE> phiV1tti(fDimension,1,0.),normalM(fDimension,1,0.);
        
        // K12 e K21 - (test V left) * (trial Lambda right)
        for(int j1 = 0; j1 < nshapeLambda; j1++)
        {
            // Var. Sigma, Sn :
            TPZFNMatrix<9, STATE> lambda_j(fDimension,1,0.);
            // Tangencial comp. vector (t x t)Sn :
            for (int e=0; e<fDimension; e++) {
                lambda_j(e,0)= datavecright[pindex].phi(j1,0)*tangent[e];
            }
            
            STATE fact = weight * InnerVec(phiVi,lambda_j);
            ek(i1,j1) +=fact;
            ek(j1,i1) +=fact;
        }
        

    }
    
    std::ofstream plotfileM("ekInterfaceH.txt");
    ek.Print("KintH = ",plotfileM,EMathematicaInput);
    
    
}



void TPZMHMBrinkmanMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    //return;
    
    STATE rhsnorm = Norm(ef);
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }
    
    
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    TPZFNMatrix<3, STATE> direction(fDimension,1,0.);
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    
    // Getting the linear combination or finite element approximations
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    
    int nshapeV, nshapeP;
    nshapeV = phiV.Rows();//datavec[vindex].fVecShapeIndex.NElements();
    nshapeP = datavec[pindex].phi.Rows();
    
    if((nshapeV != 0) && (bc.Type() < 10))
    {
        for(int i=0; i<fDimension; i++)
        {
            direction(0,0) = -datavec[vindex].axes(0,1);
            direction(1,0) = datavec[vindex].axes(0,0);
        }
    }
    else if(nshapeV != 0)
    {
        for(int i=0; i<fDimension; i++) direction(i,0) = datavec[vindex].axes(0,i);
    }
    //Adaptação para Hdiv
    int ekr= ek.Rows();
    
    //Vefifica se HDiv
    if(ekr!=nshapeP+nshapeV){
        DebugStop();
    }
    
    
    
    TPZFNMatrix<9, STATE> phiVi(fDimension,1,0.),phiVni(1,1,0.), phiVj(fDimension,1,0.),phiVnj(1,1,0.), phiPi(fDimension,1),phiPj(fDimension,1);
    
    TPZFNMatrix<3,STATE> v_2=bc.Val2();
    TPZFNMatrix<3,STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    switch (bc.Type()) {
        case 0: //Dirichlet for continuous formulation
        case 10:
        {
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                TPZFMatrix<STATE> gradu;
                bc.ForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D = vbc[2];
            }
            
            if(fSpace==1){
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    
                    REAL vh_n = v_h[0];
                    REAL v_n = direction[0] * v_2[0] + direction[1] * v_2[1];
                    
                    ef(i,0) += -weight * gBigNumber * (vh_n - v_n) * phiV(i,0);
                    
                    for(int j = 0; j < nshapeV; j++){
                        
                        ek(i,j) += weight * gBigNumber * phiV(j,0) * phiV(i,0);
                        
                    }
                    
                }
                
                
                
            }else{
                DebugStop();
            }
            
            
            
        }
            break;
            
        case 1: //Neumann for continuous formulation
        case 11:
        {
            
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                bc.ForcingFunction()->Execute(datavec[pindex].x,vbc);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D  = vbc[2]*0.;
                
                
            }
            
            
            for(int i = 0; i < nshapeV; i++ )
            {
                
                //Adaptação para Hdiv
                
                STATE val = 0;
                for(int i=0; i<fDimension; i++) val += direction(i,0)*v_2(i,0);
                
                ef(i,0) += weight * val * datavec[vindex].phi(i,0);
                
            }
            
            
        }
            
            
            
            break;
            
        case 5: //Ponto pressao
        {
            p_D = bc.Val2()(0,0);
            
            
            for(int i = 0; i < nshapeP; i++ )
            {
                
                
                ef(i) += 1.0 * p_D * phiP(i,0);
                
                for(int j = 0; j < nshapeP; j++){
                    
                    ek(i,j) += 1.0 * (phiP(i,0) * phiP(j,0));
                    
                }
                
            }
            
        }
            break;
            
            
            //        case 5: // put a value on the diagonal of the pressure to mark a reference pressure
            //            ek(0,0) += 1.;
            //            break;
            
        default:
        {
            std::cout << "Boundary not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
    
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }
    
}


void TPZMHMBrinkmanMaterial::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    return;
    DebugStop();
    
}
