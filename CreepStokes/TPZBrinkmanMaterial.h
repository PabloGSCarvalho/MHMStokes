/*
 *  TPZBrinkmanMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

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

#ifndef TPZBrinkmanMATERIAL
#define TPZBrinkmanMATERIAL


class TPZBrinkmanMaterial :
    public TPZMatBase<STATE,
                      TPZMatCombinedSpacesT<STATE>,
                      TPZMatErrorCombinedSpaces<STATE>,
                      TPZMatInterfaceCombinedSpaces<STATE>>
{
    using TBase = TPZMatBase<STATE,
                             TPZMatCombinedSpacesT<STATE>,
                             TPZMatErrorCombinedSpaces<STATE>,
                             TPZMatInterfaceCombinedSpaces<STATE>>;  

protected:
    
    /// dimension of the material
    int fDimension;
    
    /// Aproximation Space for velocity
    int fSpace;
    
    /// viscosidade
    STATE fViscosity;
    
    /** @brief Medium permeability. Coeficient which multiplies the gradient operator*/
    STATE fk;
    
    /// termo contrario a beta na sua formulacao (para ser conforme a literatura)
    STATE fTheta;
    
    STATE fSigma;
    
//    TPZTransform<STATE> f_T;
//    
//    TPZTransform<STATE> f_InvT;
    
public:
    
    bool NeedsNormalVecFad = true;
    /**
     * Empty Constructor
     */
    TPZBrinkmanMaterial();
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    [[maybe_unused]] TPZBrinkmanMaterial(int matid, int dimension, int space, STATE viscosity, STATE theta, STATE Sigma);
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZBrinkmanMaterial(const TPZBrinkmanMaterial &mat);
    
    /**
     * Destructor
     */
    ~TPZBrinkmanMaterial();
    
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
        return "TPZBrinkmanMaterial";
    }
    
    STATE GetViscosity() const{
        return fViscosity;
    }

    void SetViscosity(STATE visco){
        fViscosity = visco;
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {
        return fDimension;
    }
    
    /** returns the number of state variables associated with the material */
    virtual int NStateVariables() const override {
        
        if(fSpace==1){
            return 1;
        }else{
            return 2;
        }
        
    } // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);

    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    // Returns the solution associated with the var index based on the finite element approximation around one interface element
    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout) override {};

    void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                           const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                           const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                           int var, TPZVec<STATE> &Solout,
                           TPZCompEl *left, TPZCompEl *right) override {};

    /** index of velocity */
    int VIndex(){ return 0; }
    
    /** index of pressure */
    int PIndex(){ return 1; }
    
    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    template <class TVar2>
    TVar2 Inner(TPZFMatrix<TVar2> &S, TPZFMatrix<TVar2> &T);
    
    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    STATE InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T);
    
    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<STATE> &GradU );
    
    /** transpose of the tensor GradU = Div(U)*/
    STATE Transpose(TPZFMatrix<STATE> &GradU );
    
    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi);
    
    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);

    // Contribute Methods being used - Multiphysics
    void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    // It computes a contribution to the stiffness matrix and load vector at one BC integration point.
    void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    // It computes a contribution to the stiffness matrix and load vector at one BC interface integration point
    void ContributeBCInterface(const TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavecleft, REAL weight, 
    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    // It computes a contribution to the stiffness matrix and load vector at one internal interface integration point
    void ContributeInterface(const TPZMaterialDataT<STATE> &data, std::map<int, TPZMaterialDataT<STATE>> &datavecleft, 
    std::map<int, TPZMaterialDataT<STATE>> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

    // Save the element data to a stream
    void Write(TPZStream &buf, int withclassid);
    
    // Read the element data from a stream
    void Read(TPZStream &buf, void *context);
    
    int NEvalErrors() {return 6;}
    
    //It computes errors contribution in differents spaces.
    void Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) override;
        
};

#endif
