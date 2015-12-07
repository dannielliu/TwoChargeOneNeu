#ifndef Physics_Analysis_EE3pi_H
#define Physics_Analysis_EE3pi_H 

#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
//#include "VertexFit/ReadBeamParFromDb.h"
//#include "EventTag/EventTagSvc.h"
class EE3pi : public Algorithm {

public:
  EE3pi(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

private:



  //   ReadBeamParFromDb m_reader;
  bool m_runmc;
  bool m_mcmatch;
  double m_ecms;
  double m_espread;
  double m_beamangle;
  double mpi0;
  double mrho;
  double meta;
  double metap;
  double m_vr0cut;
  double m_vz0cut;
  double m_coscut;
  double m_p0cut;

  bool m_useBarrel;
  bool m_useEndcap;
  double m_energyBThreshold;
  double m_energyEThreshold;
  double m_gammaPhiCut;
  double m_gammaThetaCut;
  double m_gammaTrkCut;
  double m_gammathCut;
  double m_gammatlCut;
  int  m_gammanumlCut;
  int  m_gammanumhCut;



  NTuple::Tuple*  m_tuple;      // truth
  NTuple::Array<double>  m_truthppi0;
  NTuple::Array<double>  m_truthppip;
  NTuple::Array<double>  m_truthppim;
  NTuple::Item<double>  m_mcminv;

  NTuple::Item<long>  m_ncharg;
  NTuple::Item<long>  m_nneu;
  NTuple::Item<long>  m_ngch; 
  NTuple::Item<long>  m_run; 
  NTuple::Item<long>  m_rec;
  NTuple::Item<long>  m_evttag; 
  NTuple::Item<long>   m_idxmc;
  NTuple::Array<int>  m_pdgid;
  NTuple::Array<int>  m_motheridx;

  NTuple::Array<double>  m_recppip;
  NTuple::Array<double>  m_recppim;

  //  IEventTagSvc* m_EventTagSvc;


};

#endif 
