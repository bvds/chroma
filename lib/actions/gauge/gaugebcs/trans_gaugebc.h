// -*- C++ -*-
/*! \file
 *  \brief Simple gauge boundary conditions
 */

#ifndef __trans_gaugebc_h__
#define __trans_gaugebc_h__

#include "gaugebc.h"

namespace Chroma 
{ 
  
  /*! @ingroup gaugebcs */
  namespace TransGaugeBCEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }

  /*! @ingroup gaugebcs */
  struct TransGaugeBCParams 
  { 
    TransGaugeBCParams();
    TransGaugeBCParams(XMLReader& xml, const std::string& path);
    multi1d<Real> phases;
    multi1d<int> trans_dirs;
  };

  /*! @ingroup gaugebcs */
  void read(XMLReader& xml, const std::string& path, TransGaugeBCParams& p); 
  
  /*! @ingroup gaugebcs */
  void write(XMLWriter& xml, const std::string& path, const TransGaugeBCParams& p);


  //! Concrete class for gauge actions with fixed links on transverse boundaries
  /*! @ingroup gaugebcs
   *
   *  For two transverse directions, the 
   *  links are set to some fixed constant.
   *  
   */
  class TransGaugeBC : public GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:

    //! Destructor is automatic
    ~TransGaugeBC() {}

    TransGaugeBC(const TransGaugeBCParams& p); 

    //! Modify U fields in place
    // Used in lib/update/heatbath/mciter.cc
    void modify( multi1d<LatticeColorMatrix>& u) const;

    //! Zero the U fields in place on the masked links
    // Used in lib/actions/gauge/gaugeacts/
    // for various force terms.
    void zero( multi1d<LatticeColorMatrix>& u) const;

    //! Says if there are non-trivial BC links
    // Used only in clover and wilson actions, see lib/actions/ferm/linop/ 
    bool nontrivialP() const {return true;}

  private:
    multi1d<LatticeBoolean>      mask;
    multi1d<LatticeColorMatrix>  fixedLinks;

  };
}

#endif
