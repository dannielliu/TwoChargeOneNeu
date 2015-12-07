#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"



#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"




#include "McTruth/McParticle.h"
#include "McTruth/DecayMode.h"
#include "McTruth/MdcMcHit.h"
#include "McTruth/TofMcHit.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/TofMcHit.h"
#include "McTruth/EmcMcHit.h"
#include "McTruth/MucMcHit.h"
#include "McTruth/McEvent.h"


#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "EE3piAlg/EE3pi.h"
#include "VertexFit/Helix.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/IVertexDbSvc.h"
#include "ParticleID/ParticleID.h"
#include "VertexFit/SecondVertexFit.h"

#include <vector>
//const double twopi = 6.2831853;
//const double pi = 3.1415927;
//const double mpi0 = 0.1349766;
const double me = 0.000511;
const double mu = 0.105658;
const double mpi = 0.13957;
const double mka = 0.493677;
const double xmass[5] = {0.000511, 0.105658, 0.139570,0.493677, 0.938272};
//const double velc = 29.9792458;  tof_path unit in cm.
const double velc = 299.792458;   // tof path unit in mm
typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
typedef std::vector<Hep3Vector> Vp3;
static int counter[10]={0,0,0,0,0,0,0,0,0,0};
/////////////////////////////////////////////////////////////////////////////

EE3pi::EE3pi(const std::string& name, ISvcLocator* pSvcLocator) :
		Algorithm(name, pSvcLocator) {

				//Declare the properties 
				declareProperty("runmc",m_runmc = false);
				declareProperty("mcmatch",m_mcmatch = true);  
				declareProperty("ecms",m_ecms = 3.686);
				declareProperty("espread",m_espread=0.0013);
				declareProperty("beamangle",m_beamangle = 0.022); 

				declareProperty("mpi0",mpi0 = 0.1349766); 
				declareProperty("meta",meta = 0.547853); 
				declareProperty("metap",metap = 0.95766);
				declareProperty("mrho",mrho = 0.77549);
                //charged tracks
			////////declareProperty("Vr0cut", m_vr0cut=1.0);
			////////declareProperty("Vz0cut", m_vz0cut=10.0);
				declareProperty("Vr0cut", m_vr0cut=10.0);
				declareProperty("Vz0cut", m_vz0cut=30.0);
				declareProperty("Coscut", m_coscut=0.93);
                                declareProperty("MomentumCut", m_p0cut=2.0);               
                //gamma  
				declareProperty("UseBarrelEmc",    m_useBarrel        = true);
				declareProperty("UseEndcapEmc",    m_useEndcap        = true);
				declareProperty("BarrelEnergyThreshold", m_energyBThreshold=0.025);
				declareProperty("EndcapEnergyThreshold", m_energyEThreshold=0.025);
				declareProperty("GammaPhiCut", m_gammaPhiCut=20.0);
				declareProperty("GammaThetaCut", m_gammaThetaCut=20.0);
				declareProperty("GammaTrkCut", m_gammaTrkCut=0.0);
				declareProperty("GammaTimeWinHCut", m_gammathCut=14.0);
				declareProperty("GammaTimeWinLCut", m_gammatlCut=0.0);
				declareProperty("GammaNumLCut", m_gammanumlCut=2);
				declareProperty("GammaNumHCut", m_gammanumhCut=8);

	     //final 6c

		}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode EE3pi::initialize(){
	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in initialize()" << endmsg;

	StatusCode status;
	//   static const bool CREATEIFNOTTHERE(true);
	//   StatusCode EventTagSvcStatus = service("EventTagSvc", m_EventTagSvc, CREATEIFNOTTHERE);
	//   if (!EventTagSvcStatus.isSuccess() || 0 ==m_EventTagSvc ) {
	//     log << MSG::ERROR << " Could not initialize Decay code service" << endreq;
	//     return EventTagSvcStatus;
	//   }


	NTuplePtr nt4(ntupleSvc(), "GETAPETAP/EE3pi");
	if ( nt4 ) m_tuple = nt4;
	else {
		m_tuple = ntupleSvc()->book ("GETAPETAP/EE3pi", CLID_ColumnWiseTuple, "ks N-Tuple example");
		if ( m_tuple )    {
			status = m_tuple->addItem ("run",  m_run); 
			status = m_tuple->addItem ("rec",  m_rec); 
			status = m_tuple->addItem ("evttag",  m_evttag); 
			status = m_tuple->addItem("indexmc",          m_idxmc, 0, 100);
			status = m_tuple->addIndexedItem("pdgid",     m_idxmc, m_pdgid);
			status = m_tuple->addIndexedItem("motheridx", m_idxmc, m_motheridx);

			status = m_tuple->addItem ("ngch" ,   m_ngch);
			status = m_tuple->addItem ("ncharg",   m_ncharg);
			status = m_tuple->addItem ("nneu",    m_nneu);
			
			status = m_tuple->addItem ("mcpi0", 4, m_truthppi0);
			status = m_tuple->addItem ("mcpip", 4, m_truthppip);
			status = m_tuple->addItem ("mcpim", 4, m_truthppim);
			status = m_tuple->addItem ("mcminv", m_mcminv);
			
			status = m_tuple->addItem ("recpip", 4, m_recppip);
			status = m_tuple->addItem ("recpim", 4, m_recppim);
			
		}
		else    { 
			log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple) << endmsg;
			return StatusCode::FAILURE;
		}
	}
	//
	//--------end of book--------
	//

	log << MSG::INFO << "successfully return from initialize()" <<endmsg;
	return StatusCode::SUCCESS;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
StatusCode EE3pi::execute() {
	const double beamEnergy = m_ecms/2.;
	const HepLorentzVector p_cms(m_ecms*sin(m_beamangle*0.5),0.0,0.0,m_ecms);
	//const HepLorentzVector p_cms(0.040547,0,0,3.68632);
	const Hep3Vector u_cms = -p_cms.boostVector(); 



	StatusCode sc=StatusCode::SUCCESS;

	MsgStream log(msgSvc(), name());
	log << MSG::INFO << "in execute()" << endreq;

	counter[0]++;


	SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
	if (!eventHeader) 
	{
		log << MSG::FATAL << "Could not find Event Header" << endreq;
		return StatusCode::SUCCESS;
	}

	m_run = eventHeader->runNumber();
	m_rec = eventHeader->eventNumber();
	m_evttag=eventHeader->eventTag();
	//if( m_run >= 9755 && m_run<= 9766 ) return StatusCode::SUCCESS;
	//if( m_run >= 34298 && m_run<= 34313 ) return StatusCode::SUCCESS;
	//if( m_run >= 34461 && m_run<= 34477 ) return StatusCode::SUCCESS;
	//if( m_run >= 35086 && m_run<= 35098 ) return StatusCode::SUCCESS;


	Hep3Vector vgammamc_ext, vpi0_gammamc_high_ext, vpi0_gammamc_low_ext;
	if (eventHeader->runNumber()<0)
	{
		//MC information
		//
		Vp4 mcppip,mcppim;

		mcppip.clear();
		mcppim.clear();
		
		SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(), "/Event/MC/McParticleCol");

		int m_numParticle = 0;
		Vint mcidxi;
		std::vector<int> mcidx, pdgid;
		mcidx.clear();
		pdgid.clear();
		if (!mcParticleCol)
		{
			std::cout << "Could not retrieve McParticelCol" << std::endl;
			return StatusCode::FAILURE;
		}

		bool psipDecay = false;
		//int rootIndex = -1;
		int isrtag=0;
		int KKtag = 0;
		HepLorentzVector mcpip, mcpim, mcpi0;		

		Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
		for (; iter_mc != mcParticleCol->end(); iter_mc++)
		{
			HepLorentzVector  mctrue_track = (*iter_mc)->initialFourMomentum();
			HepLorentzVector mctrack_iniposition = (*iter_mc)->initialPosition();
			HepLorentzVector mctrack_finposition = (*iter_mc)->finalPosition();

	////////	if ((*iter_mc)->primaryParticle()) continue;
	////////	if (!(*iter_mc)->decayFromGenerator()) continue;
	////////	//if ( ((*iter_mc)->mother()).trackIndex()<3 ) continue;
	////////	int imom=443;
	////////	if(m_ecms<3.2)imom=443;
	////////	else if (m_ecms>3.6)imom=100443;
	////////	if ((*iter_mc)->particleProperty()==imom)  
	////////		//if ((*iter_mc)->particleProperty()==100443) 
	////////	{
	////////		psipDecay = true;
	////////		rootIndex = (*iter_mc)->trackIndex();
	////////		//    m_evttag=(m_EventTagSvc->getChainCode(*iter_mc))<<8;

	////////	}
	////////	if (!psipDecay) continue;
			mcidxi.push_back((*iter_mc)->trackIndex());
			mcidx.push_back(((*iter_mc)->mother()).trackIndex());
			pdgid.push_back((*iter_mc)->particleProperty());
		////////long mcidxi=(*iter_mc)->trackIndex() - rootIndex;	
		////////long mcidx = ((*iter_mc)->mother()).trackIndex() - rootIndex;
		////////long pdgid = (*iter_mc)->particleProperty();
		////////m_motheridx[m_numParticle] = (long)mcidx;
		////////m_pdgid[m_numParticle] = (long)pdgid;
		////////m_numParticle += 1; 

		        //if ((*iter_mc)->particleProperty() == 321 || (*iter_mc)->particleProperty() == 211){
		        if ((*iter_mc)->particleProperty() == 321){
                         	m_truthppip[0] = mctrue_track.px();
		        	m_truthppip[1] = mctrue_track.py();
		        	m_truthppip[2] = mctrue_track.pz();
		        	m_truthppip[3] = mctrue_track.e();
                                mcpip = mctrue_track;
				//cout<<"Mother_ID"<<((*iter_mc)->mother()).particleProperty()<<endl;
		        } 
		        //if ((*iter_mc)->particleProperty() == -321||(*iter_mc)->particleProperty() == -211){
		        if ((*iter_mc)->particleProperty() == -321){
		        	m_truthppim[0] = mctrue_track.px();
		        	m_truthppim[1] = mctrue_track.py();
		        	m_truthppim[2] = mctrue_track.pz();
		        	m_truthppim[3] = mctrue_track.e();
                                mcpim = mctrue_track;
		        }
		        
			if ((*iter_mc)->particleProperty() == 111){
		        	m_truthppi0[0] = mctrue_track.px();
		        	m_truthppi0[1] = mctrue_track.py();
		        	m_truthppi0[2] = mctrue_track.pz();
		        	m_truthppi0[3] = mctrue_track.e();
                                mcpi0 = mctrue_track;
		        }

 
		}
		//cout<<"mcIGam.size() ="<<mcIGam.size()<<endl;
		//cout<<"mcIEta.size() ="<<mcIEta.size()<<endl;

		m_numParticle = mcidx.size();
		m_idxmc = m_numParticle;
		for (int i=0;i<mcidx.size();i++){
			m_motheridx[i] = (short)mcidx[i];
			m_pdgid[i] = (short)pdgid[i];
			//m_numParticle += 1; 
		}
		m_mcminv = (mcpip+mcpim+mcpi0).m();



		//		m_tuple0->write();               
	}

	//cout<<"vpi0_gammamc_high_ext ="<<vpi0_gammamc_high_ext<<endl;
	//cout<<"vpi0_gammamc_low_ext ="<<vpi0_gammamc_low_ext<<endl;

	SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
	if (!evtRecEvent) 
	{
		log << MSG::FATAL << "Could not find EvtRecEvent" << endreq;
		return StatusCode::SUCCESS;
	}
	log << MSG::INFO <<"ncharg, nneu, tottks = " 
		<< evtRecEvent->totalCharged() << " , "
		<< evtRecEvent->totalNeutral() << " , "
		<< evtRecEvent->totalTracks() <<endreq;
	//  if(evtRecEvent->totalNeutral()>30)return sc;

	m_ncharg  = evtRecEvent->totalCharged();

	m_nneu = evtRecEvent->totalNeutral();



        HepPoint3D vx(0., 0., 0.);
        HepSymMatrix Evx(3, 0);

        IVertexDbSvc*  vtxsvc;
        Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
        if(vtxsvc->isVertexValid()){
        	double* dbv = vtxsvc->PrimaryVertex();
        	double*  vv = vtxsvc->SigmaPrimaryVertex();
        	//	if (m_reader.isRunNumberValid( m_run)) {
        	//   HepVector dbv = m_reader.PrimaryVertex( m_run);
        	//    HepVector vv = m_reader.SigmaPrimaryVertex( m_run);
        	vx.setX(dbv[0]);
        	vx.setY(dbv[1]);
        	vx.setZ(dbv[2]);
        	Evx[0][0]=vv[0]*vv[0];
        	Evx[0][1]=vv[0]*vv[1];
        	Evx[1][1]=vv[1]*vv[1];
        	Evx[1][2]=vv[1]*vv[2];
        	Evx[2][2]=vv[2]*vv[2];
        }

        VertexParameter vx_db;
        vx_db.setVx(vx);
        vx_db.setEvx(Evx);




	SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
	if (!evtRecTrkCol) 
	{
		log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
		return StatusCode::SUCCESS;
	}

	Vint iGood;
	iGood.clear();

	int nCharge = 0;
	for(int i = 0; i < evtRecEvent->totalCharged(); i++){
		EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
		if(!(*itTrk)->isMdcTrackValid()) continue;
		//if(!(*itTrk)->isMdcKalTrackValid()) continue;
		RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
		double pch=mdcTrk->p();
		double x0=mdcTrk->x();
		double y0=mdcTrk->y();
		double z0=mdcTrk->z();
		double costheta = mdcTrk->theta();
		double phi_mdc = mdcTrk->phi();
		double phi0=mdcTrk->helix(1);
		double xv=vx.x();
		double yv=vx.y();
		double zv=vx.z();
		double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);
////////	double m_vx0 = x0;
////////	double m_vy0 = y0;
////////	double m_vz0 = z0;
////////	double m_vr0 = Rxy;
		HepVector a = mdcTrk->helix();
		HepSymMatrix Ea = mdcTrk->err();
		HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
		HepPoint3D IP(xv,yv,zv); 
		VFHelix helixip(point0,a,Ea); 
		helixip.pivot(IP);
		HepVector vecipa = helixip.a();
		double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
		double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
		double  Rvphi0=vecipa[1];
////////	m_rvxy0=Rvxy0;
////////	m_rvz0=Rvz0;
		//		m_rvphi0=Rvphi0;
		//m_costheta = costheta;
		//		m_phi_mdc = phi_mdc;
////////	m_pt = pch;
////////	m_rvx0 = Rvxy0 * cos(Rvphi0);
////////	m_rvy0 = Rvxy0 * sin(Rvphi0);
////////	m_phi0 = Rvphi0;	
		//		m_tuple1->write();
		if(fabs(Rvz0) >= m_vz0cut) continue;
		if(fabs(Rvxy0) >= m_vr0cut) continue;
		//		if(pch >= m_p0cut) continue;
		double cost = cos(mdcTrk->theta());
		if(fabs(cost) >= m_coscut ) continue;

		iGood.push_back(i);//   iGood.push_back((*itTrk)->trackId());
		nCharge += mdcTrk->charge();

	}


	//
	// Finish Good Charged Track Selection
	//
	int nGood = iGood.size();
	m_ngch=nGood;
	log << MSG::DEBUG << "ngood, totcharge = " << nGood << " , " << nCharge << endreq;

	if( nGood != 2 || nCharge !=0){
		return StatusCode::SUCCESS;
	}
	counter[1]++;
	
	ParticleID *pid = ParticleID::instance();
	int nk=0;

	for (int i=0; i<2; i++){
		EvtRecTrackIterator itTrk1=evtRecTrkCol->begin() + iGood.at(i);
		RecMdcTrack *piTrk = (*itTrk1)->mdcTrack();

		pid->init();
		pid->setMethod(pid->methodProbability());
		pid->setChiMinCut(4);
		pid->setRecTrack(*itTrk1);
		pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2()); // use PID sub-system
		pid->identify(pid->onlyPionKaonProton());
		pid->calculate();
		if(!(pid->IsPidInfoValid())) continue;
		if((pid->probKaon() < pid->probProton()) || (pid->probKaon() < pid->probPion())) continue;
		nk++;

		if (piTrk->charge()>0) {
			HepLorentzVector p4pip = piTrk->p4(xmass[0]);
			m_recppip[0] = p4pip.px();
			m_recppip[1] = p4pip.py();
			m_recppip[2] = p4pip.pz();
			m_recppip[3] = p4pip.e();
        	}
		else {
			HepLorentzVector p4pim = piTrk->p4(xmass[0]);
			m_recppim[0] = p4pim.px();
			m_recppim[1] = p4pim.py();
			m_recppim[2] = p4pim.pz();
			m_recppim[3] = p4pim.e();
		}
	}
	if (nk!=2) return StatusCode::SUCCESS;
	counter[2]++;

	m_tuple->write();

	return StatusCode::SUCCESS;
	}


	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
	StatusCode EE3pi::finalize() {
		cout<<"Alg is  : EE3pi   "<<endl;
		cout<<"total number:         "<<counter[0]<<endl;
		cout<<"good track  :         "<<counter[1]<<endl;
		cout<<"pid 3 gt:             "<<counter[2]<<endl;
		cout<<"vertex valid:         "<<counter[3]<<endl;
		MsgStream log(msgSvc(), name());
		log << MSG::INFO << "in finalize()" << endmsg;
		return StatusCode::SUCCESS;
	}

