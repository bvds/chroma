// $Id: nef_quarkprop4_w.h,v 1.4 2004-12-29 22:13:41 edwards Exp $
/*! \file
 * \brief Full quark propagator solver for domain wall fermions
 *
 * Given a complete propagator as a source, this does all the inversions needed
 */

#ifndef __nef_quarkprop4_w_h__
#define __nef_quarkprop4_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"

namespace Chroma
{
  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to Domain Wall fermions (Array) fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param t_src    time slice of source ( Read )
   * \param j_decay  direction of decay ( Read )
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  void
  nef_quarkProp4(LatticePropagator& q_sol, 
		 XMLWriter& xml_out,
		 const LatticePropagator& q_src,
		 int t_src, int j_decay,
		 const UnprecDWFermActBaseArray<LatticeFermion>& S_f,
		 Handle<const ConnectState> state,
		 const InvertParam_t& invParam,
		 int& ncg_had);

  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to Domain Wall fermions (Array) fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param t_src    time slice of source ( Read )
   * \param j_decay  direction of decay ( Read )
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  void
  nef_quarkProp4(LatticePropagator& q_sol, 
		 XMLWriter& xml_out,
		 const LatticePropagator& q_src,
		 int t_src, int j_decay,
		 const EvenOddPrecDWFermActBaseArray<LatticeFermion>& S_f,
		 Handle<const ConnectState> state,
		 const InvertParam_t& invParam,
		 int& ncg_had);

}

#endif
