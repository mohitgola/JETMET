// -*- C++ -*-
//
//
// Package:    jet/MiniAnalyzer
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc jet/MiniAnalyzer/plugins/MiniAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mohit Gola
//         Created:  Wed, 20 Jan 2016 03:01:06 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "TTree.h"
#include "TFile.h"
//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MiniAnalyzer(const edm::ParameterSet&);
      ~MiniAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

     edm::EDGetTokenT<pat::JetCollection> jetToken_;
     edm::EDGetTokenT<pat::METCollection> metToken_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

TFile * file;
TTree* tree;
std::vector<double> NHFa;
std::vector<double> NEMFa;
std::vector<double> CHFa;
std::vector<double> MUFa;
std::vector<double> CEMFa;
std::vector<double> NumCountsa;
std::vector<double> NumNeutralParticlesa;
std::vector<double> CHMa;
std::vector<double> jet_eta;
std::vector<double> jet_pt;
std::vector<double> jet_phi;
std::vector<double> met_eta;
std::vector<double> met_pt;
std::vector<double> met_phi;
std::vector<double> met_SumEt;

std::vector<double> NHFb;
std::vector<double> NEMFb;
std::vector<double> CHFb;
std::vector<double> MUFb;
std::vector<double> CEMFb;
std::vector<double> NumCountsb;
std::vector<double> NumNeutralParticlesb;
std::vector<double> CHMb;











};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
 jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
 metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))

{
   //now do what ever initialization is needed

}


MiniAnalyzer::~MiniAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);

for(const pat::Jet &j : *jets)

{
float  NHF1 = j.neutralHadronEnergyFraction();
float  NHF2 = j.neutralHadronEnergyFraction();

float NEMF1 = j.neutralEmEnergyFraction();
float NEMF2 = j.neutralEmEnergyFraction();

float CHF1 = j.chargedHadronEnergyFraction();
float CHF2 = j.chargedHadronEnergyFraction();

float MUF1 = j.muonEnergyFraction();
float MUF2 = j.muonEnergyFraction();

float CEMF1 = j.chargedEmEnergyFraction();
float CEMF2 = j.chargedEmEnergyFraction();

float NumCounts1 = j.chargedMultiplicity()+j.neutralMultiplicity();
float NumCounts2 = j.chargedMultiplicity()+j.neutralMultiplicity();

float NumNeutralParticles1 =j.neutralMultiplicity();
float NumNeutralParticles2 =j.neutralMultiplicity();

float CHM1 = j.chargedMultiplicity();
float CHM2 = j.chargedMultiplicity();  
jet_eta.push_back(j.eta());
jet_pt.push_back(j.pt());
jet_phi.push_back(j.phi()); 


 if(fabs(j.eta())<=3.0 )
	{

		if ((NHF1<0.99 && NEMF1<0.99 && NumCounts1>1) && ((fabs(j.eta())<=2.4 && CHF1>0 && CHM1>0 && CEMF1<0.99) || fabs(j.eta())>2.4) && fabs(j.eta())<=3.0 )

		{	NHFa.push_back(NHF1);
                	NEMFa.push_back(NEMF1);
                	CHFa.push_back(CHF1);
                	MUFa.push_back(MUF1);
	                CEMFa.push_back(CEMF1);
        	        NumCountsa.push_back(NumCounts1);
        		NumNeutralParticlesa.push_back(NumNeutralParticles1);
	                CHMa.push_back(CHM1);
		}


	}



 if(fabs(j.eta())>3.0 )
        {

	if(NEMF2<0.90 && NumNeutralParticles2>10 && fabs(j.eta())>3.0 ) 
                {       NHFb.push_back(NHF2);
                        NEMFb.push_back(NEMF2);
                        CHFb.push_back(CHF2);
                        MUFb.push_back(MUF2);
                        CEMFb.push_back(CEMF2);
                        NumCountsb.push_back(NumCounts2);
                        NumNeutralParticlesb.push_back(NumNeutralParticles2);
                        CHMb.push_back(CHM2);
                }


        }

}

	
		edm::Handle<pat::METCollection> mets;
                iEvent.getByToken(metToken_, mets);

for(const pat::MET &m : *mets)

	{
		met_eta.push_back(m.eta());
		met_pt.push_back(m.pt());
                met_phi.push_back(m.phi());
                met_SumEt.push_back(m.sumEt());

	}








	tree->Fill();


jet_eta.clear();
jet_pt.clear();
jet_phi.clear();
met_eta.clear();
met_pt.clear();
met_phi.clear();
met_SumEt.clear();
NHFa.clear();
NEMFa.clear();
CHFa.clear(); 
MUFa.clear();
CEMFa.clear();
NumCountsa.clear();
NumNeutralParticlesa.clear();
CHMa.clear();

NHFb.clear();
NEMFb.clear();
CHFb.clear();
MUFb.clear();
CEMFb.clear();
NumCountsb.clear();
NumNeutralParticlesb.clear();
CHMb.clear();





#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob()
{

file = new TFile("test.root","RECREATE");
tree = new TTree("T","Electron and muon");
tree->Branch("jet_eta",&jet_eta);
tree->Branch("jet_pt",&jet_pt);
tree->Branch("jet_phi",&jet_phi);
tree->Branch("met_eta",&met_eta);
tree->Branch("met_pt",&met_pt);
tree->Branch("met_phi",&met_phi);
tree->Branch("met_SumEt",&met_SumEt);

tree->Branch("NHFa",&NHFa);
tree->Branch("NEMFa",&NEMFa);
tree->Branch("CHFa",&CHFa);
tree->Branch("CEMFa",&CEMFa);
tree->Branch("NumCountsa",&NumCountsa);
tree->Branch("NumNeutralParticlesa",&NumNeutralParticlesa);
tree->Branch("MUFa",&MUFa);
tree->Branch("CHMa",&CHMa);
tree->Branch("NHFb",&NHFb);
tree->Branch("NEMFb",&NEMFb);
tree->Branch("CHFb",&CHFb);
tree->Branch("CEMFb",&CEMFb);
tree->Branch("NumCountsb",&NumCountsb);
tree->Branch("NumNeutralParticlesb",&NumNeutralParticlesb);
tree->Branch("MUFb",&MUFb);
tree->Branch("CHMb",&CHMb);

}




// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzer::endJob() 
{
file->Write();
file->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
MiniAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MiniAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MiniAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MiniAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
