/*
 *  TPZDarcyPMaterial.h
 *  PZ
 *
 *  Created by Pablo G S Carvalho on 08/09/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZDARCYPMATERIAL
#define TPZDARCYPMATERIAL

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

class TPZDarcyPMaterial :
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
    
    int fDimension;
    
    /// Aproximation Space for velocity
    int fSpace;
    
    /// viscosidade
    STATE fViscosity;
    
    /** @brief Medium permeability. Coeficient which multiplies the gradient operator*/
    REAL fk;
    
    /// termo contrario a beta na sua formulacao (para ser conforme a literatura)
    STATE fTheta;
    
    
public:
    
    TPZDarcyPMaterial();
    
    [[maybe_unused]] TPZDarcyPMaterial(int matid, int dimension, int space, STATE viscosity, STATE permeability, STATE theta);
    
    TPZDarcyPMaterial(const TPZDarcyPMaterial &mat);
    
    ~TPZDarcyPMaterial();
    
    //void FillDataRequirements(TPZMaterialDataT<STATE> &data) const override;
    
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    void FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data,
                                       std::map<int, TPZMaterialDataT<STATE>> &datavec_left,
                                       std::map<int, TPZMaterialDataT<STATE>> &datavec_right) override;

    void FillBoundaryConditionDataRequirements(int type,TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
    void SetPermeability(REAL perm) {
        fk = perm;
    }
    
    /** returns the name of the material */
    std::string Name() {
        return "TPZDarcyPMaterial";
    }
    
    int Dimension() const {return 2;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() const override {return 4;} // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
    
    void Print(std::ostream &out = std::cout);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout) override;

    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout,
                           TPZCompEl *left, TPZCompEl *right) override;

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
    void FillGradPhi(TPZMaterialDataT<STATE> &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi);
    
    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialDataT<STATE> &data);
       
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
