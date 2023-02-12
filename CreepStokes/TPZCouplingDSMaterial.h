/*
 *  TPZCouplingDSMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 11/10/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */


#include "TPZMatWithMem.h"
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "TPZBndCond.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZMatBase.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatInterfaceCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "pzfunction.h"

#ifndef TPZCOUPLINGDSMATERIAL
#define TPZCOUPLINGDSMATERIAL

class TPZCouplingDSMaterial :
    public TPZMatBase<STATE,
                      TPZMatCombinedSpacesT<STATE>,
                      TPZMatErrorCombinedSpaces<STATE>,
                      TPZMatInterfaceCombinedSpaces<STATE>>
{
    using TBase = TPZMatBase<STATE,
                             TPZMatCombinedSpacesT<STATE>,
                             TPZMatErrorCombinedSpaces<STATE>,
                             TPZMatInterfaceCombinedSpaces<STATE>>;    

private:

    STATE fViscosity;
    
    // Medium permeability. Coeficient which multiplies the gradient operator
    STATE fk;
    
    // Termo contrario a beta na sua formulacao (para ser conforme a literatura)
    STATE fTheta;

    int fDimension;
public:
    
    TPZCouplingDSMaterial();
    
    [[maybe_unused]] TPZCouplingDSMaterial(int matid, int dimension, STATE viscosity,STATE permeability, STATE theta);
    
    TPZCouplingDSMaterial(const TPZCouplingDSMaterial &mat);
    
    ~TPZCouplingDSMaterial();
    
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    void FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data,
                                       std::map<int, TPZMaterialDataT<STATE>> &datavec_left,
                                       std::map<int, TPZMaterialDataT<STATE>> &datavec_right) override;

    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    void SetPermeability(REAL perm) {
        fk = perm;
    }
    
    std::string Name() {
        return "TPZStokesMaterial";
    }
    
    int Dimension() const {return 2;}
    
    /** returns the number of state variables associated with the material */
    virtual int NStateVariables() const override {return 4;} // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
    
    void Print(std::ostream &out = std::cout);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout) override {};

    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout,
                           TPZCompEl *left, TPZCompEl *right) override {};
    
    int VIndex(){ return 0; }

    int PIndex(){ return 1; }
    
    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    template <typename TVar2>
    TVar2 Inner(TPZFMatrix<TVar2> &S, TPZFMatrix<TVar2> &T);
    
    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    template <typename TVar2>
    TVar2 InnerVec(TPZFMatrix<TVar2> &S, TPZFMatrix<TVar2> &T);
    
    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<REAL> &GradU );
    
    /** transpose of the tensor GradU = Div(U)*/
    STATE Transpose(TPZFMatrix<REAL> &GradU );
    
    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi);
    
    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);
       
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    void ContributeBCInterface(const TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavecleft, REAL weight, 
    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    void ContributeInterface(const TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavecleft, 
    std::map<int, TPZMaterialDataT<STATE>> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    void Write(TPZStream &buf, int withclassid);
    
    void Read(TPZStream &buf, void *context);
    
    int NEvalErrors() {return 6;}
    
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
        
};

#endif
