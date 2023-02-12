/*
 *  TPZStokesMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZSTOKESMATERIAL
#define TPZSTOKESMATERIAL

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

class TPZStokesMaterial :
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
    
public:
    
    
    /**
     * Empty Constructor
     */
    TPZStokesMaterial();
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    [[maybe_unused]] TPZStokesMaterial(int matid, int dimension, int space, STATE viscosity, STATE theta, STATE Sigma);
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZStokesMaterial(const TPZStokesMaterial &mat);
    
    /**
     * Destructor
     */
    ~TPZStokesMaterial();
    
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
        return "TPZStokesMaterial";
    }

    /** returns the integrable dimension of the material */
    int Dimension() const {return 2;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() const override {return 4;} // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout) override;

    // Returns the solution associated with the var index based on the finite element approximation around one interface element
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
