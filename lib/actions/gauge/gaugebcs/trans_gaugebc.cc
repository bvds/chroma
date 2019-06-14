/*! \file
 *  \brief Simple gauge boundary conditions
 */

#include "actions/gauge/gaugebcs/trans_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"


namespace Chroma {

  namespace TransGaugeBCEnv 
  { 
    //! Calllback function to register with the factory
    GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeBC(XMLReader& xml, 
										       const std::string& path)
    {
      QDPIO::cout << "Factory Callback: Creating TransGaugeBC " << std::endl;
      return new TransGaugeBC(TransGaugeBCParams(xml, path));
    }

    const std::string name = "TRANS_GAUGEBC";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeBCFactory::Instance().registerObject(name, createGaugeBC);
	registered = true;
      }
      return success;
    }
  }

  TransGaugeBCParams::TransGaugeBCParams(XMLReader& xml, 
					   const std::string& path) {
      
    XMLReader paramtop(xml, path);
    
    try {
      // List of diagonal phases for a link
      read(paramtop, "phases", phases);
      // The two directions to apply the boundary conditions.
      read(paramtop, "transverse_dirs", trans_dirs);
    }
    catch( const std::string& e ) { 
      QDPIO::cerr << "Error reading XML: " << e << std::endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Creating TransGaugeBCParams with phases: ";
    for(int i=0; i<Nc; i++)
      QDPIO::cout << phases[i] << " ";
    QDPIO::cout << " in directions: " << trans_dirs[0]
		<< " " << trans_dirs[1] << std::endl;
  }


  void read(XMLReader& xml, const std::string& path, TransGaugeBCParams& p) {
    TransGaugeBCParams tmp(xml, path);
    p=tmp;
  }
  
  void write(XMLWriter& xml, const std::string& path, const TransGaugeBCParams& p) { 
    push(xml, path);
    write(xml, "phases", p.phases);
    write(xml, "transverse_dirs", p.trans_dirs);
    pop(xml);
  }

  TransGaugeBC::TransGaugeBC(const TransGaugeBCParams& p)
  {
    mask.resize(Nd);
    fixedLinks.resize(Nd);
    // Set the mask and default link values.
    for(int mu=0; mu<Nd; ++mu)
      {
	mask[mu] = false;
	// Set Longitudinal links on the boundary to 1.
	fixedLinks[mu] = 1;
	for(int i=0; i<2; ++i)
	  {
	    int nu = p.trans_dirs[i];
	    // All lattice sites, labeled by the nu-direction coordinate.
	    LatticeInteger lnu = Layout::latticeCoordinate(nu);
	    /* Mask links on the surface of boundary
	     * Also, mask links that are outside the boundary.
	     */
	    mask[mu] |= lnu==Layout::lattSize()[nu]-1;
	    if(mu != nu)
	      mask[mu] |= lnu==0;
	  }
      }

    /*
     * Set the links winding around the transverse boundary.
     * Loops winding once around the transverse direction add up
     * to the phase set in the xml file.
     * Note the given phases should obey |lambda_i - lambda_j| < 2 pi
     * and sum_i lambda_i = 0.
     */
    Real rescale = 1.0;
    for(int i=0; i<2; ++i)
      rescale /= 2.0*(Layout::lattSize()[p.trans_dirs[i]]-1);
    for(int i=0; i<2; ++i)
      {
	int mu = p.trans_dirs[i];
	// Mask for links outside the boundary
	LatticeInteger llong = Layout::latticeCoordinate(mu);
	LatticeBoolean notEdge = llong==Layout::lattSize()[mu]-1;
	// Want phases to add up as we loop around the transverse boundary.
	// Mask for links where we take the complex conjugate.
	LatticeInteger lperp = Layout::latticeCoordinate(p.trans_dirs[1-i]);
	LatticeBoolean oddEdge = (lperp==0) != (i==0);
	LatticeReal phase;
	LatticeReal tmpPhase;
	for(int j = 0; j < Nc;  ++j)
	  {
	    phase = p.phases[j]*rescale;
	    // Take the complex conjugate
	    tmpPhase = -phase;
	    copymask(phase, oddEdge, tmpPhase);
	    // Set links ouside the boundary to 1.
	    tmpPhase = 0;
	    copymask(phase, notEdge, tmpPhase);
	    pokeColor(fixedLinks[mu], cmplx(cos(phase), sin(phase)), j, j);
	  }
      }
  }

  //! Modify U fields in place
  void TransGaugeBC::modify( multi1d<LatticeColorMatrix>& u) const
  {
    START_CODE();

    for(int mu=0; mu < u.size(); ++mu)
      copymask(u[mu], mask[mu], fixedLinks[mu]);
 
    END_CODE();
  }
  
  //! Zero the U fields in place on the masked links
  void TransGaugeBC::zero( multi1d<LatticeColorMatrix>& u) const
  {
    START_CODE();

    LatticeColorMatrix z = QDP::zero;

    for(int mu=0; mu < u.size(); ++mu)
      copymask(u[mu], mask[mu], z);

    END_CODE();
  }

} // End namespace Chroma 
