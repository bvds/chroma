// $Id: t_prec_contfrac.cc,v 1.1 2004-10-01 17:48:31 bjoo Exp $

#include <iostream>
#include <cstdio>
#include "chroma.h"

using namespace QDP;
using namespace std;
using namespace Chroma;

struct App_input_t {
  Cfg_t        cfg;
  std::string  stateInfo;
  multi1d<int> nrow;
  multi1d<int> boundary;
  UnprecOvlapContFrac5DFermActParams p_unprec;
  EvenOddPrecOvlapContFrac5DFermActParams p_prec;
};

// Reader for input parameters
void read(XMLReader& xml, const string& path, App_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read the input
  try
  {
    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    XMLReader xml_state_info(inputtop, "StateInfo");
    std::ostringstream os;
    xml_state_info.print(os);
    input.stateInfo = os.str();
    
    read(inputtop, "boundary", input.boundary);
    read(inputtop, "nrow", input.nrow);
    read(inputtop, "UnprecFermAct", input.p_unprec);
    read(inputtop, "PrecFermAct", input.p_prec);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }
}


int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  App_input_t input;
  XMLReader xml_in("DATA");

  try {
    read(xml_in, "/ContFracTest", input);
  }
   catch( const string& e) { 
    QDPIO::cerr << "Caught Exception : " << e << endl;
    QDP_abort(1);
  }


  // Setup the lattice
  Layout::setLattSize(input.nrow);
  Layout::create();

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;
  gaugeStartup(gauge_file_xml, gauge_xml, u, input.cfg);

  XMLFileWriter xml_out("XMLDAT");
  push(xml_out,"ContFracTest");


  // Measure the plaquette on the gauge
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "plaquette");
  write(xml_out, "w_plaq", w_plaq);
  write(xml_out, "s_plaq", s_plaq);
  write(xml_out, "t_plaq", t_plaq);
  write(xml_out, "link", link);
  pop(xml_out);

  // Create a FermBC
  Handle< FermBC< multi1d< LatticeFermion> > >  fbc(new SimpleFermBC< multi1d< LatticeFermion> >(input.boundary));
 
  // Initialize fermion actions
  UnprecOvlapContFrac5DFermActArray S_unprec(fbc, input.p_unprec);
  EvenOddPrecOvlapContFrac5DFermActArray S_prec(fbc, input.p_prec);


  // Create an overlap state
  std::istringstream state_info_is(input.stateInfo);
  XMLReader state_info_xml(state_info_is);
  string state_info_path="/StateInfo";

  Handle< const ConnectState > state(S_prec.createState(u, state_info_xml, state_info_path));

  // Make an unprec linOp
  Handle< const LinearOperator< multi1d<LatticeFermion> > > M_u( S_unprec.linOp(state) );
  
  // Make the prec linOp
  Handle< const EvenOddPrecLinearOperator< multi1d<LatticeFermion> > > M_e( S_prec.linOp(state));

  QDPIO::cout << "Unprec LinOp size = " << M_u->size() << endl;
  QDPIO::cout << "Prec   LinOp size = " << M_e->size() << endl;

  int N5 = M_u->size();
  multi1d<LatticeFermion> s(N5);
  multi1d<LatticeFermion> Mu_s(N5);
  multi1d<LatticeFermion> Me_s(N5);
  multi1d<LatticeFermion> r(N5);

  for(int i=0; i<N5; i++) { 
    gaussian(s[i]);
    Mu_s[i]=zero;
    Me_s[i]=zero;
    r[i] = zero;
  }

  (*M_u)(Mu_s, s, PLUS);
  (*M_e).unprecLinOp(Me_s, s, PLUS);

  for(int i=0; i < N5; i++) { 
    r[i] = Mu_s[i] - Me_s[i];
    QDPIO::cout << "i = " << i << " || r [" << i << "] || = " 
		<< sqrt(norm2(r[i])) << endl;
  }

  Mu_s = zero;
  Me_s = zero;
 (*M_u)(Mu_s, s, MINUS);
  (*M_e).unprecLinOp(Me_s, s, MINUS);

  for(int i=0; i < N5; i++) { 
    r[i] = Mu_s[i] - Me_s[i];
    QDPIO::cout << "i = " << i << " || r [" << i << "] || = " 
		<< sqrt(norm2(r[i])) << endl;
  }

  // Now some solver tests
  multi1d<LatticeFermion> source(N5);
  for(int i=0; i < N5; i++) { 
    gaussian(source[i]);
  }
  Double source_norm = sqrt(norm2(source));
  for(int i=0; i < N5; i++) { 
    source[i] /= source_norm;
  }

  multi1d<LatticeFermion> unprec_source(N5);
  multi1d<LatticeFermion> prec_source(N5);
  for(int i=0; i < N5; i++) { 
    unprec_source[i] = source[i];
    prec_source[i] = source[i];
  }

  multi1d<LatticeFermion> unprec_soln(N5);
  multi1d<LatticeFermion> prec_soln(N5);
  unprec_soln = zero;
  prec_soln = zero;

  Real RsdCG = Real(1.0e-6);
  int  MaxCG = 500;
  int  n_count_prec=0;
  int  n_count_unprec=0;

  // Solve the unpreconditioned system:
  {
    multi1d<LatticeFermion> tmp(N5);
    (*M_u)(tmp, unprec_source, MINUS);
    InvCG2(*M_u, tmp, unprec_soln, RsdCG, MaxCG, n_count_unprec);
  }

  // Solve the preconditioned system.
  //
  // Set up source on odd checkerboards
  //  S_e = source_e 
  //  S_o = source_o - QoeQee^{-1} source_e
  {
    multi1d<LatticeFermion> tmp(N5);
    multi1d<LatticeFermion> tmp2(N5);
    tmp=zero;
    tmp2=zero;
    M_e->evenEvenInvLinOp(tmp, prec_source, PLUS);
    M_e->oddEvenLinOp(tmp2, tmp, PLUS);

    multi1d<LatticeFermion> tmp_source(N5);
    for(int i=0; i < N5; i++) { 
      // Take the even piece of the normal source
      tmp_source[i][rb[0]] = prec_source[i];

      // Take the odd piece  - QoeQee^{-1} even piece
      tmp_source[i][rb[1]] = prec_source[i] - tmp2[i];
    } 

    // CGNE on the odd source only
    multi1d<LatticeFermion> cgne_source(N5);
    multi1d<LatticeFermion> tmp_soln(N5);
    cgne_source=zero;
    tmp_soln=zero;

    (*M_e)(cgne_source, tmp_source, MINUS);
    InvCG2(*M_e, cgne_source, tmp_soln, RsdCG, MaxCG, n_count_prec);

    // Now reconstruct the solutions. 
    tmp = zero;
    tmp2 = zero;
    // tmp1 = D_eo tmp_soln_o
    M_e->evenOddLinOp(tmp, tmp_soln, PLUS);
    
    // tmp2 = S_e - tmp_1
    for(int i=0; i < N5; i++) { 
      tmp2[i][rb[0]] = tmp_source[i] - tmp[i];
    }
    M_e->evenEvenInvLinOp(prec_soln, tmp2, PLUS);
    
    // Copy off parts
    for(int i=0; i < N5; i++) { 
     prec_soln[i][rb[1]] = tmp_soln[i];
    }
  }

  // Check the solutions
  (*M_u)(r, unprec_soln, PLUS);
  Double r_norm =Double(0);
  for(int i=0; i < N5; i++) {
    r[i] -= source[i];
    r_norm += norm2(r[i]);
  }
  QDPIO::cout << "|| source - M_u unprec_soln || = " << sqrt(r_norm) << endl;

  (*M_u)(r, prec_soln, PLUS);
  r_norm=0;
  for(int i=0; i < N5; i++) {
    r[i] -= source[i];
    r_norm += norm2(r[i]);
  }
  QDPIO::cout << "|| source - M_u prec_soln || = " << sqrt(r_norm) << endl;

 (*M_e).unprecLinOp(r, unprec_soln, PLUS);
  r_norm=0;
  for(int i=0; i < N5; i++) {
    r[i] -= source[i];
    r_norm += norm2(r[i]);
  }
  QDPIO::cout << "|| source - M_e unprec_soln || = " << sqrt(r_norm) << endl;

 (*M_e).unprecLinOp(r, prec_soln, PLUS);
  r_norm=0;
  for(int i=0; i < N5; i++) {
    r[i] -= source[i];
    r_norm += norm2(r[i]);
  }
  QDPIO::cout << "|| source - M_e prec_soln || = " << sqrt(r_norm) << endl;

  // Test EvenEvenInv LinOp
  {
    multi1d<LatticeFermion> f_even(N5);
    multi1d<LatticeFermion> t1(N5);
    multi1d<LatticeFermion> t2(N5);
    
    f_even = zero;
    t1=zero;
    t2=zero;
    for(int i=0; i < N5; i++) { 
      f_even[i][rb[1]] = zero ; // Zilch out the odd ones
      f_even[i][rb[0]] = source[i];
    }
    
    M_e->evenEvenLinOp(t1, f_even, PLUS);
    M_e->evenEvenInvLinOp(t2, t1, PLUS);

    r_norm=0;
    for(int i=0; i < N5; i++) {
      r[i][rb[0]] = t2[i];
      r[i][rb[0]] -= f_even[i];

      r_norm += norm2(r[i], rb[0]);
  }
  QDPIO::cout << "|| source_even - Qee^{-1}Qee sorce_even || = " 
	      << sqrt(r_norm) << endl;

  }
  pop(xml_out);
  QDP_finalize();
    
  exit(0);
}
