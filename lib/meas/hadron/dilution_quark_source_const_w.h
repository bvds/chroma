// -*- C++ -*-
// $Id: dilution_quark_source_const_w.h,v 1.1 2007-12-14 06:53:42 edwards Exp $
/*! \file
 * \brief Dilution operator constructed from MAKE_SOURCE and PROPAGATOR calls
 */

#ifndef __dilution_quark_source_const_h__
#define __dilution_quark_source_const_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace DilutionQuarkSourceConstEnv 
  {
    extern const std::string name;
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      //! Solution files for each quark
      struct QuarkFiles_t
      {
	//! Time dilution components  
	struct TimeDilutions_t
	{
	  struct SpinDilutions_t
	  {
	    multi1d<std::string> dilution_files;  /*!< dilution files for this spin and time dilution*/
	  };

	  multi1d<SpinDilutions_t> spin_files;  /*!< Spin dilution files for this time dilution */
	};
	
	multi1d<TimeDilutions_t> time_files;
      };

      QuarkFiles_t  quark;       /*!< All the solutions that are needed for a single quark */
    }; // struct Params


    //! Dilutions constructed for propagator solutions over MAKE_SOURCE dilutions
    class Dilute : public DilutionOperator<LatticeFermion>
    {
    public:
      //! Virtual destructor to help with cleanup;
      ~Dilute() {}

      //! Get an iterator for this dilution
      /*! This uses the t0 for the time slice where the creation operator is created */
      const_iterator begin(int t0) const;

      //! Get an iterator for this dilution
      const_iterator end() const;

      //! The decay direction
      int getDecayDir() const;

      //! The seed identifies this quark
      const Seed& getSeed() const;

      //! Is this element of the dilution zero?
      bool zero(const_iterator iter) const;

      //! Return the original source vector
      /*! NOTE: this might be slow */
      LatticeFermion source() const;
    
      //! Return the original source vector
      LatticeFermion dilutedSource(const_iterator iter) const;
    
      //! Return the original source vector
      LatticeFermion dilutedSolution(const_iterator iter) const;
    };
    
  } // namespace 
}

#endif
