// -*- C++ -*-
// $Id: inline_plaquette.h,v 3.3 2007-04-09 18:41:57 bjoo Exp $
/*! \file
 *  \brief Inline plaquette
 */

#ifndef __inline_plaquette_h__
#define __inline_plaquette_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlinePlaquetteEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  /*! \ingroup inlineglue */
  struct InlinePlaquetteParams 
  {
    InlinePlaquetteParams();
    InlinePlaquetteParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      GroupXML_t    cgs;      /*!< Gauge State */
    } param;

    struct NamedObject_t
    {
      std::string   gauge_id;
    } named_obj;

    std::string xml_file;
  };


  /*! \ingroup inlineglue */
  class InlinePlaquette : public AbsInlineMeasurement 
  {
  public:
    ~InlinePlaquette() {}
    InlinePlaquette(const InlinePlaquetteParams& p) : params(p) {}
    InlinePlaquette(const InlinePlaquette& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    void operator()(unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    void func(const unsigned long update_no, 
	      XMLWriter& xml_out);

  private:


    InlinePlaquetteParams params;
  };

};

#endif
