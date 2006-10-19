// -*- C++ -*-
// $Id: unprec_ovlap_contfrac5d_nonhermop_array_w.h,v 3.1 2006-10-19 16:01:32 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) linear operator
 */

#ifndef __unprec_ovlap_contfrac5d_nonhermop_array_w_h__
#define __unprec_ovlap_contfrac5d_nonhermop_array_w_h__

#include "linearop.h"
#include "unprec_wilstype_fermact_w.h" 


namespace Chroma 
{ 
  //! Unpreconditioned Extended-Overlap (N&N) linear operator
  /*!
   * \ingroup linop
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+m_q)/(1-m_q)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   *
   * This operator implements  hep-lat/0005004
   */

  class UnprecOvlapContFrac5DNonHermOpArray : public UnprecLinearOperatorArray<LatticeFermion, 
					      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    UnprecOvlapContFrac5DNonHermOpArray(const UnprecWilsonTypeFermAct<T,P,Q>& S_aux,
					Handle< FermState<T,P,Q> > state,
					const Real& _m_q,
					int _N5,
					const Real& _scale_fac,
					const multi1d<Real>& _alpha,
					const multi1d<Real>& _beta,
					int _NEig,
					const multi1d<Real>& _EigValFunc,
					const multi1d<LatticeFermion>& _EigVec) :
      M(S_aux.linOp(state)), fbc(state->getFermBC()), 
      m_q(_m_q), N5(_N5), scale_fac(_scale_fac), alpha(_alpha),
      beta(_beta), NEig(_NEig), EigVec(_EigVec), EigValFunc(_EigValFunc)
      {}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Destructor is automatic
    ~UnprecOvlapContFrac5DNonHermOpArray() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

    //! Only defined on the entire lattice
    const OrderedSubset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (multi1d<LatticeFermion>& chi, 
		     const multi1d<LatticeFermion>& psi, 
		     enum PlusMinus isign) const;

  private:
    Handle< LinearOperator<LatticeFermion> > M;  // Wilson(esque) Op
    Handle< FermBC<T,P,Q> >  fbc;
    const Real m_q;
    const int  N5;    // Size of the 5th dimension
    const Real scale_fac;
    const multi1d<Real> alpha;
    const multi1d<Real> beta;
    const int NEig;
    const multi1d<LatticeFermion> EigVec;
    const multi1d<Real> EigValFunc;
  };

} // End Namespace Chroma


#endif
