// Class:      HLTTauAnalyzer
//
// Original Author:  Sandeep Bhowmik
//         Created:  Tue, 12 Mar 2019 18:38:39 GMT
//
#include "FWCore/Framework/interface/one/EDAnalyzer.h" 
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"     
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include <FWCore/Framework/interface/Frameworkfwd.h>
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

#include "DataFormats/VertexReco/interface/Vertex.h"  
#include "DataFormats/L1Trigger/interface/Tau.h" 
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <DataFormats/PatCandidates/interface/Tau.h> 

#include "DataFormats/TauReco/interface/PFTau.h"                // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"             // reco::PFTauRef, reco::PFTauCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"   // reco::PFTauDiscriminator

#include "DataFormats/L1TParticleFlow/interface/HPSPFTau.h"     // l1t::L1HPSPFTau
#include "DataFormats/L1TParticleFlow/interface/HPSPFTauFwd.h"  // l1t::HPSPFTauCollection
#include "DataFormats/L1TParticleFlow/interface/PFTau.h"       // l1t::PFTauCollection for NNTau

class L1andHLTTauAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit L1andHLTTauAnalyzer(const edm::ParameterSet&);
  ~L1andHLTTauAnalyzer();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  void Initialize();
  
  TTree *tree_;
  std::string treeName_;
  
  ULong64_t       indexevents_;
  Int_t           runNumber_;
  Int_t           lumi_;
  float MC_weight_;
  double stitching_weight_;
  double genVertex_;
  std::vector<float> genTauPt_;
  std::vector<float> genTauEta_;
  std::vector<float> genTauPhi_;
  std::vector<int> genTauCharge_;
  std::vector<Bool_t> isGenMatched_;
  std::vector<int> recoTauDecayMode_;
  std::vector<float> recoTauPt_;
  std::vector<float> recoTauEta_;
  std::vector<float> recoTauPhi_;
  std::vector<int> recoTauCharge_;
  std::vector<Bool_t> isRecoMatched_;
  std::vector<int> recoGMTauDecayMode_;
  std::vector<float> recoGMTauPt_;
  std::vector<float> recoGMTauEta_;
  std::vector<float> recoGMTauPhi_;
  std::vector<int> recoGMTauCharge_;
  std::vector<Bool_t> isRecoGMMatched_;
  std::vector<int> hltTauType_;
  std::vector<float> hltTauPt_;
  std::vector<float> hltTauEta_;
  std::vector<float> hltTauPhi_;
  std::vector<int> hltTauCharge_;
  std::vector<float> hltTauIso_;
  std::vector<Bool_t> hltTauTightIso_;
  std::vector<Bool_t> hltTauMediumIso_;
  std::vector<Bool_t> hltTauLooseIso_;
  std::vector<Bool_t> hltTauVLooseIso_;
  std::vector<Bool_t> hltTauTightRelIso_;
  std::vector<Bool_t> hltTauMediumRelIso_;
  std::vector<Bool_t> hltTauLooseRelIso_;
  std::vector<Bool_t> hltTauVLooseRelIso_;
  std::vector<float> hltTauZ_;
  std::vector<float> hltTauLeadTrackPt_;
  std::vector<int> l1hpsTauType_;
  std::vector<float> l1hpsTauPt_;
  std::vector<float> l1hpsTauEta_;
  std::vector<float> l1hpsTauPhi_;
  std::vector<int> l1hpsTauCharge_;
  std::vector<float> l1hpsTauIso_;
  std::vector<Bool_t> l1hpsTauTightIso_;
  std::vector<Bool_t> l1hpsTauMediumIso_;
  std::vector<Bool_t> l1hpsTauLooseIso_;
  std::vector<Bool_t> l1hpsTauVLooseIso_;
  std::vector<Bool_t> l1hpsTauTightRelIso_;
  std::vector<Bool_t> l1hpsTauMediumRelIso_;
  std::vector<Bool_t> l1hpsTauLooseRelIso_;
  std::vector<Bool_t> l1hpsTauVLooseRelIso_;
  std::vector<float> l1hpsTauZ_;
  std::vector<float> l1hpsTauLeadTrackPt_;
  std::vector<int> l1nnTauType_;
  std::vector<float> l1nnTauPt_;
  std::vector<float> l1nnTauEta_;
  std::vector<float> l1nnTauPhi_;
  std::vector<int> l1nnTauCharge_;
  std::vector<float> l1nnTauIso_;
  std::vector<Bool_t> l1nnTauTightIso_;
  std::vector<Bool_t> l1nnTauMediumIso_;
  std::vector<Bool_t> l1nnTauLooseIso_;
  std::vector<Bool_t> l1nnTauVLooseIso_;
  std::vector<Bool_t> l1nnTauTightRelIso_;
  std::vector<Bool_t> l1nnTauMediumRelIso_;
  std::vector<Bool_t> l1nnTauLooseRelIso_;
  std::vector<Bool_t> l1nnTauVLooseRelIso_;
  std::vector<float> l1nnTauZ_;
  std::vector<float> l1nnTauLeadTrackPt_;


  bool createHistRoorFile_;
  std::string histRootFileName_;
  TFile* histRootFile_;
  TH1F* hist_genTauPt_;
  TH1F* hist_genTauEta_;
  TH1F* hist_genTauPhi_;
  TH1F* hist_isGenMatched_;
  TH1F* hist_recoTauPt_;
  TH1F* hist_recoTauEta_;
  TH1F* hist_recoTauPhi_;
  TH1F* hist_isRecoMatched_;
  TH1F* hist_recoGMTauPt_;
  TH1F* hist_recoGMTauEta_;
  TH1F* hist_recoGMTauPhi_;
  TH1F* hist_isRecoGMMatched_;
  TH1F* hist_hltTauPt_;
  TH1F* hist_hltTauEta_;
  TH1F* hist_hltTauPhi_;
  TH1F* hist_hltTauReso_vs_Gen_;
  TH1F* hist_hltTauReso_vs_Reco_;
  TH1F* hist_hltTauReso_vs_RecoGM_;

  // ----------member data ---------------------------

  bool debug_;
  bool isGenTau_;
  bool isRecoTau_;
  bool isHLTTau_;
  bool isL1HPSTau_;
  bool isL1NNTau_;
  double min_pt_;
  double max_eta_;
  edm::EDGetTokenT<double>                         evtWeightToken_;
  edm::EDGetTokenT<GenEventInfoProduct>            genTagToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>>      genTauToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>>      recoVertexToken_;
  edm::EDGetTokenT<std::vector<pat::Tau>>          recoTauToken_;
  edm::EDGetTokenT<pat::TauRefVector>              recoGMTauToken_;
  edm::EDGetTokenT<pat::TauCollection>             hltTauToken_;
  std::string                                      hltTauSumChargedIsoToken_;
  edm::EDGetTokenT<l1t::HPSPFTauCollection>      l1hpsTauToken_;
  edm::EDGetTokenT<l1t::PFTauCollection>           l1nnTauToken_;
};


L1andHLTTauAnalyzer::L1andHLTTauAnalyzer(const edm::ParameterSet& iConfig)
  : debug_          (iConfig.getUntrackedParameter<bool>("debug", false))
  , isGenTau_       (iConfig.getUntrackedParameter<bool>("isGenTau", false))
  , isRecoTau_      (iConfig.getUntrackedParameter<bool>("isRecoTau", false))
  , isHLTTau_       (iConfig.getUntrackedParameter<bool>("isHLTTau", false))
  , isL1HPSTau_     (iConfig.getUntrackedParameter<bool>("isL1HPSTau", false))
  , isL1NNTau_      (iConfig.getUntrackedParameter<bool>("isL1NNTau", false))
  , min_pt_         (iConfig.getUntrackedParameter<double>("min_pt", 0))
  , max_eta_        (iConfig.getUntrackedParameter<double>("max_eta", 3.0))
  , evtWeightToken_ (consumes<double>                       (iConfig.getParameter<edm::InputTag>("src_evtWeight")))
  , genTagToken_    (consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("genTagToken")))
  , genTauToken_    (consumes<std::vector<reco::GenJet>>    (iConfig.getParameter<edm::InputTag>("genTauToken")))
  , recoVertexToken_(consumes<std::vector<reco::Vertex>>    (iConfig.getParameter<edm::InputTag>("recoVertexToken")))
  , recoTauToken_   (consumes<std::vector<pat::Tau>>        (iConfig.getParameter<edm::InputTag>("recoTauToken")))
  , recoGMTauToken_ (consumes<pat::TauRefVector>            (iConfig.getParameter<edm::InputTag>("recoGMTauToken")))
  , hltTauToken_    (consumes<pat::TauCollection>           (iConfig.getParameter<edm::InputTag>("hltTauToken")))
  , hltTauSumChargedIsoToken_ (iConfig.getParameter<std::string>("hltTauSumChargedIsoToken"))
  , l1hpsTauToken_  (consumes<l1t::HPSPFTauCollection>    (iConfig.getParameter<edm::InputTag>("l1hpsTauToken")))
  , l1nnTauToken_   (consumes<l1t::PFTauCollection>         (iConfig.getParameter<edm::InputTag>("l1nnTauToken")))
{
   //now do what ever initialization is needed
  treeName_             = iConfig.getParameter<std::string>("treeName");
  edm::Service<TFileService> fs;
  tree_                 = fs -> make<TTree>(treeName_.c_str(), treeName_.c_str());
  createHistRoorFile_   = iConfig.getUntrackedParameter<bool>("createHistRoorFile", false);
  histRootFileName_     = iConfig.getParameter<std::string>("histRootFileName");
  histRootFile_         = new TFile(histRootFileName_.c_str(), "RECREATE");
  return;
}


L1andHLTTauAnalyzer::~L1andHLTTauAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1andHLTTauAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Initialize();

  if(debug_){
    std::cout<<" Starting L1 nad HLTTau Analyzer ............     "<< std::endl;  
  }
    using namespace edm;

   indexevents_ = iEvent.id().event();
   runNumber_ = iEvent.id().run();
   lumi_ = iEvent.luminosityBlock();

   //edm::Handle<double>             evtWeight;
   //iEvent.getByToken(evtWeightToken_, evtWeight);
   //stitching_weight_ = *evtWeight;
   //cout<<"evtWeight "<<*evtWeight<<endl;

   if(isGenTau_){
     edm::Handle<GenEventInfoProduct>        genEvt;
     try {iEvent.getByToken(genTagToken_,    genEvt);}
     catch (...) {;}
     if(genEvt.isValid()) MC_weight_ = genEvt->weight();

     edm::Handle<std::vector<reco::GenJet>>  genTauHandle;
     iEvent.getByToken(genTauToken_,         genTauHandle);

     for(auto genTau : *genTauHandle){
       if (fabs(genTau.eta())>max_eta_)
	 continue;
       //if (fabs(genTau.pt())<min_pt_)
       //continue;
       genTauPt_.push_back(genTau.pt());
       genTauEta_.push_back(genTau.eta());
       genTauPhi_.push_back(genTau.phi());
       genTauCharge_.push_back(genTau.charge());

       hist_genTauPt_->Fill(genTau.pt());
       hist_genTauEta_->Fill(genTau.eta());
       hist_genTauPhi_->Fill(genTau.phi());

       if(debug_){
	 std::cout<<" GenTau pt "<<genTau.pt()<<" eta "<< genTau.eta()<<" phi "<< genTau.phi()<<" charge "<< genTau.charge()<<std::endl;
       }
     }
   }

   if(isL1HPSTau_){
     edm::Handle<l1t::HPSPFTauCollection>  l1hpsTauHandle;
     iEvent.getByToken(l1hpsTauToken_,        l1hpsTauHandle);
     for(auto l1hpsTau : *l1hpsTauHandle){
       l1hpsTauPt_.push_back(l1hpsTau.pt());
       l1hpsTauEta_.push_back(l1hpsTau.eta());
       l1hpsTauPhi_.push_back(l1hpsTau.phi());
       l1hpsTauCharge_.push_back(l1hpsTau.charge());
       l1hpsTauType_.push_back(l1hpsTau.tauType());
       l1hpsTauIso_.push_back(l1hpsTau.sumChargedIso());
       l1hpsTauTightIso_.push_back(l1hpsTau.passTightIso());
       l1hpsTauMediumIso_.push_back(l1hpsTau.passMediumIso());
       l1hpsTauLooseIso_.push_back(l1hpsTau.passLooseIso());
       l1hpsTauVLooseIso_.push_back(l1hpsTau.passVLooseIso());
       if(l1hpsTau.pt()!=0){
	 if(l1hpsTau.sumChargedIso()/l1hpsTau.pt() < 0.40){
	   l1hpsTauVLooseRelIso_.push_back(true);
	 }else{
	   l1hpsTauVLooseRelIso_.push_back(false);
	 }
	 if(l1hpsTau.sumChargedIso()/l1hpsTau.pt() < 0.20){
	   l1hpsTauLooseRelIso_.push_back(true);
	 }else{
	   l1hpsTauLooseRelIso_.push_back(false);
	 }
	 if(l1hpsTau.sumChargedIso()/l1hpsTau.pt() < 0.10){
	   l1hpsTauMediumRelIso_.push_back(true);
	 }else{
	   l1hpsTauMediumRelIso_.push_back(false);
	 }
	 if(l1hpsTau.sumChargedIso()/l1hpsTau.pt() < 0.05){
	   l1hpsTauTightRelIso_.push_back(true);
	 }else{
	   l1hpsTauTightRelIso_.push_back(false);
	 }
       }
       
       double l1hpsTauZ = 1000;
       if ( l1hpsTau.leadChargedPFCand().isNonnull() && l1hpsTau.leadChargedPFCand()->pfTrack().isNonnull()){
	 l1hpsTauZ = l1hpsTau.leadChargedPFCand()->pfTrack()->vertex().z();
       }
       l1hpsTauZ_.push_back(l1hpsTauZ);
       if(debug_){
	 std::cout<<" L1 HPS Tau pt "<<l1hpsTau.pt()<<" eta "<< l1hpsTau.eta()<<" phi "<< l1hpsTau.phi()<<" charge "<< l1hpsTau.charge()<<" Type "<< l1hpsTau.tauType()<<" sumChargedIso "<< l1hpsTau.sumChargedIso()<<std::endl;
	 std::cout<<" L1 HPS Tau Z " << l1hpsTauZ  << std::endl;
       }
     }
   }//if(isL1HPSTau_){

   if(isL1NNTau_){
     edm::Handle<l1t::PFTauCollection>       l1nnTauHandle;
     iEvent.getByToken(l1nnTauToken_,        l1nnTauHandle);
     for(auto l1nnTau : *l1nnTauHandle){
       l1nnTauPt_.push_back(l1nnTau.pt());
       l1nnTauEta_.push_back(l1nnTau.eta());
       l1nnTauPhi_.push_back(l1nnTau.phi());
       l1nnTauCharge_.push_back(l1nnTau.charge());
       l1nnTauType_.push_back(1);
       l1nnTauIso_.push_back(l1nnTau.chargedIso());
       l1nnTauTightIso_.push_back(l1nnTau.passTightNN());
       l1nnTauMediumIso_.push_back(l1nnTau.passTightPF());
       l1nnTauLooseIso_.push_back(l1nnTau.passLooseNN());
       l1nnTauVLooseIso_.push_back(l1nnTau.passLoosePF());
       if(l1nnTau.pt()!=0){
	 if(l1nnTau.chargedIso()/l1nnTau.pt() < 0.40){
	   l1nnTauVLooseRelIso_.push_back(true);
	 }else{
	   l1nnTauVLooseRelIso_.push_back(false);
	 }
	 if(l1nnTau.chargedIso()/l1nnTau.pt() < 0.20){
	   l1nnTauLooseRelIso_.push_back(true);
	 }else{
	   l1nnTauLooseRelIso_.push_back(false);
	 }
	 if(l1nnTau.chargedIso()/l1nnTau.pt() < 0.10){
	   l1nnTauMediumRelIso_.push_back(true);
	 }else{
	   l1nnTauMediumRelIso_.push_back(false);
	 }
	 if(l1nnTau.chargedIso()/l1nnTau.pt() < 0.05){
	   l1nnTauTightRelIso_.push_back(true);
	 }else{
	   l1nnTauTightRelIso_.push_back(false);
	 }
       }
       double l1nnTauZ = 1000;
       l1nnTauZ_.push_back(l1nnTauZ);
       if(debug_){
	 std::cout<<" L1NNTau pt "<<l1nnTau.pt()<<" eta "<< l1nnTau.eta()<<" phi "<< l1nnTau.phi()<<" charge "<< l1nnTau.charge()<<" Type "<< 1 <<" chargedIso "<< l1nnTau.chargedIso()<<std::endl;
       }
     }
   }//if(isL1NNTau_){
   
   
   if(isHLTTau_){
     edm::Handle<pat::TauCollection>            hltTauHandle;
     iEvent.getByToken(hltTauToken_,            hltTauHandle);
     
     size_t numHLTTaus = hltTauHandle->size();
     for (size_t idxHLTTau=0; idxHLTTau<numHLTTaus; ++idxHLTTau){
       pat::TauRef hltTau(hltTauHandle, idxHLTTau);
       hltTauPt_.push_back(hltTau->pt());
       hltTauEta_.push_back(hltTau->eta());
       hltTauPhi_.push_back(hltTau->phi());
       hltTauCharge_.push_back(hltTau->charge());
       double sumChargedIso = hltTau->tauID(hltTauSumChargedIsoToken_);
       hltTauIso_.push_back(sumChargedIso);
       //std::cout<<" sumChargedIso = "<< sumChargedIso<<std::endl;
       /*
       if(hltTau->pt()!=0){
	 if(sumChargedIso > 0.60){
	   hltTauVLooseRelIso_.push_back(true);
	 }else{
	   hltTauVLooseRelIso_.push_back(false);
	 }
	 if(sumChargedIso > 0.70){
	   hltTauLooseRelIso_.push_back(true);
	 }else{
	   hltTauLooseRelIso_.push_back(false);
	 }
	 if(sumChargedIso > 0.80){
	   hltTauMediumRelIso_.push_back(true);
	 }else{
	   hltTauMediumRelIso_.push_back(false);
	 }
	 if(sumChargedIso > 0.90){
	   hltTauTightRelIso_.push_back(true);
	 }else{
	   hltTauTightRelIso_.push_back(false);
	 }
       }
       */
       //for ChargedIso
       if(hltTau->pt()!=0){
	 if(sumChargedIso/hltTau->pt() < 0.40){
	   hltTauVLooseRelIso_.push_back(true);
	 }else{
	   hltTauVLooseRelIso_.push_back(false);
	 }
	 if(sumChargedIso/hltTau->pt() < 0.20){
	   hltTauLooseRelIso_.push_back(true);
	 }else{
	   hltTauLooseRelIso_.push_back(false);
	 }
	 if(sumChargedIso/hltTau->pt() < 0.10){
	   hltTauMediumRelIso_.push_back(true);
	 }else{
	   hltTauMediumRelIso_.push_back(false);
	 }
	 if(sumChargedIso/hltTau->pt() < 0.05){
	   hltTauTightRelIso_.push_back(true);
	 }else{
	   hltTauTightRelIso_.push_back(false);
	 }
       }
       double hltTauZ = 1000;
       double hltTauLeadTrackPt = 0;
       if (hltTau->leadChargedHadrCand().isNonnull() && hltTau->leadChargedHadrCand()->bestTrack()){
	 hltTauZ = hltTau->leadChargedHadrCand()->bestTrack()->vertex().z();
	 hltTauLeadTrackPt = hltTau->leadChargedHadrCand()->bestTrack()->pt();
       }
       hltTauZ_.push_back(hltTauZ);
       hltTauLeadTrackPt_.push_back(hltTauLeadTrackPt);
       
       hist_hltTauPt_->Fill(hltTau->pt());
       hist_hltTauEta_->Fill(hltTau->eta());
       hist_hltTauPhi_->Fill(hltTau->phi());
       if(debug_){
	 std::cout<<" hltTau pt "<<hltTau->pt()<<" eta "<< hltTau->eta()<<" phi "<< hltTau->phi()<<" charge "<< hltTau->charge()<<std::endl;
	 std::cout<<" hltTau Z "<< hltTauZ <<std::endl;
       }
     }
   }//if(isHLTTau_){
   
   if(isRecoTau_){
     edm::Handle<std::vector<pat::Tau>>      recoTauHandle;
     iEvent.getByToken(recoTauToken_,        recoTauHandle);

     edm::Handle<pat::TauRefVector>          recoGMTauHandle;
     iEvent.getByToken(recoGMTauToken_,      recoGMTauHandle);

     for(auto recoTau : *recoTauHandle){
       if (fabs(recoTau.eta())>max_eta_)
	 continue;
       recoTauPt_.push_back(recoTau.pt());
       recoTauEta_.push_back(recoTau.eta());
       recoTauPhi_.push_back(recoTau.phi());
       recoTauCharge_.push_back(recoTau.charge());
       recoTauDecayMode_.push_back(recoTau.decayMode());

       hist_recoTauPt_->Fill(recoTau.pt());
       hist_recoTauEta_->Fill(recoTau.eta());
       hist_recoTauPhi_->Fill(recoTau.phi());
       
       if(debug_){
	 std::cout<<" RecoTau pt "<<recoTau.pt()<<" eta "<< recoTau.eta()<<" phi "<< recoTau.phi()<<" charge "<< recoTau.charge()<<" DecayMode "<< recoTau.decayMode()<<std::endl;
       }
     } //for(auto recoTau : *recoTauHandle){
     
     for(auto recoGMTau : *recoGMTauHandle){
       if (fabs(recoGMTau->eta())>max_eta_)
	 continue;
       recoGMTauPt_.push_back(recoGMTau->pt());
       recoGMTauEta_.push_back(recoGMTau->eta());
       recoGMTauPhi_.push_back(recoGMTau->phi());
       recoGMTauCharge_.push_back(recoGMTau->charge());
       recoGMTauDecayMode_.push_back(recoGMTau->decayMode());
       
       hist_recoGMTauPt_->Fill(recoGMTau->pt());
       hist_recoGMTauEta_->Fill(recoGMTau->eta());
       hist_recoGMTauPhi_->Fill(recoGMTau->phi());
       
       if(debug_){
	 std::cout<<" RecoGMTau pt "<<recoGMTau->pt()<<" eta "<< recoGMTau->eta()<<" phi "<< recoGMTau->phi()<<" charge "<< recoGMTau->charge()<<" DecayMode "<< recoGMTau->decayMode()<<std::endl;
       }
     } //for(auto recoGMTau : *recoGMTauHandle){
   } //if(isRecoTau_){
   






   tree_ -> Fill();
}

void L1andHLTTauAnalyzer::Initialize() {
  indexevents_ = 0;
  runNumber_ = 0;
  lumi_ = 0;
  MC_weight_ = 1;
  stitching_weight_ = 1;
  genVertex_ = 0;
  genTauPt_ .clear();
  genTauEta_ .clear();
  genTauPhi_ .clear();
  genTauCharge_ .clear();
  isGenMatched_ .clear();
  recoTauPt_ .clear();
  recoTauEta_ .clear();
  recoTauPhi_ .clear();
  recoTauCharge_ .clear();
  isRecoMatched_ .clear();
  recoTauDecayMode_ .clear();
  recoGMTauPt_ .clear();
  recoGMTauEta_ .clear();
  recoGMTauPhi_ .clear();
  recoGMTauCharge_ .clear();
  isRecoGMMatched_ .clear();
  recoGMTauDecayMode_ .clear();
  hltTauPt_ .clear();
  hltTauEta_ .clear();
  hltTauPhi_ .clear();
  hltTauCharge_ .clear();
  hltTauType_ .clear();
  hltTauIso_ .clear();
  hltTauTightIso_ .clear();
  hltTauMediumIso_ .clear();
  hltTauLooseIso_ .clear();
  hltTauVLooseIso_ .clear();
  hltTauTightRelIso_ .clear();
  hltTauMediumRelIso_ .clear();
  hltTauLooseRelIso_ .clear();
  hltTauVLooseRelIso_ .clear();
  hltTauZ_ .clear();
  hltTauLeadTrackPt_.clear();
  hltTauType_ .clear();

  l1hpsTauPt_ .clear();
  l1hpsTauEta_ .clear();
  l1hpsTauPhi_ .clear();
  l1hpsTauCharge_ .clear();
  l1hpsTauType_ .clear();
  l1hpsTauIso_ .clear();
  l1hpsTauTightIso_ .clear();
  l1hpsTauMediumIso_ .clear();
  l1hpsTauLooseIso_ .clear();
  l1hpsTauVLooseIso_ .clear();
  l1hpsTauTightRelIso_ .clear();
  l1hpsTauMediumRelIso_ .clear();
  l1hpsTauLooseRelIso_ .clear();
  l1hpsTauVLooseRelIso_ .clear();
  l1hpsTauZ_ .clear();
  l1hpsTauLeadTrackPt_.clear();
  l1hpsTauType_ .clear();
  l1nnTauPt_ .clear();
  l1nnTauEta_ .clear();
  l1nnTauPhi_ .clear();
  l1nnTauCharge_ .clear();
  l1nnTauType_ .clear();
  l1nnTauIso_ .clear();
  l1nnTauTightIso_ .clear();
  l1nnTauMediumIso_ .clear();
  l1nnTauLooseIso_ .clear();
  l1nnTauVLooseIso_ .clear();
  l1nnTauTightRelIso_ .clear();
  l1nnTauMediumRelIso_ .clear();
  l1nnTauLooseRelIso_ .clear();
  l1nnTauVLooseRelIso_ .clear();
  l1nnTauZ_ .clear();
  l1nnTauLeadTrackPt_.clear();
  l1nnTauType_ .clear();


}


// ------------ method called once each job just before starting event loop  ------------
void
L1andHLTTauAnalyzer::beginJob()
{
  tree_ -> Branch("EventNumber",&indexevents_,"EventNumber/l");
  tree_ -> Branch("RunNumber",&runNumber_,"RunNumber/I");
  tree_ -> Branch("lumi",&lumi_,"lumi/I");
  tree_ -> Branch("MC_weight",&MC_weight_,"MC_weight/F");
  tree_ -> Branch("stitching_weight", &stitching_weight_, "stitching_weight/D");
  tree_ -> Branch("genTauPt",  &genTauPt_);
  tree_ -> Branch("genTauEta", &genTauEta_);
  tree_ -> Branch("genTauPhi", &genTauPhi_);
  tree_ -> Branch("genTauCharge", &genTauCharge_);
  tree_ -> Branch("isGenMatched", &isGenMatched_);
  tree_ -> Branch("recoTauPt",  &recoTauPt_);
  tree_ -> Branch("recoTauEta", &recoTauEta_);
  tree_ -> Branch("recoTauPhi", &recoTauPhi_);
  tree_ -> Branch("recoTauCharge", &recoTauCharge_);
  tree_ -> Branch("isRecoMatched", &isRecoMatched_);
  tree_ -> Branch("recoTauDecayMode", &recoTauDecayMode_);
  tree_ -> Branch("recoGMTauPt",  &recoGMTauPt_);
  tree_ -> Branch("recoGMTauEta", &recoGMTauEta_);
  tree_ -> Branch("recoGMTauPhi", &recoGMTauPhi_);
  tree_ -> Branch("recoGMTauCharge", &recoGMTauCharge_);
  tree_ -> Branch("isRecoGMMatched", &isRecoGMMatched_);
  tree_ -> Branch("recoGMTauDecayMode", &recoGMTauDecayMode_);
  tree_ -> Branch("hltTauPt",  &hltTauPt_);
  tree_ -> Branch("hltTauEta", &hltTauEta_);
  tree_ -> Branch("hltTauPhi", &hltTauPhi_);
  tree_ -> Branch("hltTauCharge", &hltTauCharge_);
  tree_ -> Branch("hltTauType", &hltTauType_);
  tree_ -> Branch("hltTauIso", &hltTauIso_);
  tree_ -> Branch("hltTauTightIso", &hltTauTightIso_);
  tree_ -> Branch("hltTauMediumIso", &hltTauMediumIso_);
  tree_ -> Branch("hltTauLooseIso", &hltTauLooseIso_);
  tree_ -> Branch("hltTauVLooseIso", &hltTauVLooseIso_);
  tree_ -> Branch("hltTauTightRelIso", &hltTauTightRelIso_);
  tree_ -> Branch("hltTauMediumRelIso", &hltTauMediumRelIso_);
  tree_ -> Branch("hltTauLooseRelIso", &hltTauLooseRelIso_);
  tree_ -> Branch("hltTauVLooseRelIso", &hltTauVLooseRelIso_);
  tree_ -> Branch("hltTauZ", &hltTauZ_);
  tree_ -> Branch("hltTauLeadTrackPt",  &hltTauLeadTrackPt_);

  tree_ -> Branch("l1hpsTauType", &l1hpsTauType_);
  tree_ -> Branch("l1hpsTauPt",  &l1hpsTauPt_);
  tree_ -> Branch("l1hpsTauEta", &l1hpsTauEta_);
  tree_ -> Branch("l1hpsTauPhi", &l1hpsTauPhi_);
  tree_ -> Branch("l1hpsTauCharge", &l1hpsTauCharge_);
  tree_ -> Branch("l1hpsTauType", &l1hpsTauType_);
  tree_ -> Branch("l1hpsTauIso", &l1hpsTauIso_);
  tree_ -> Branch("l1hpsTauTightIso", &l1hpsTauTightIso_);
  tree_ -> Branch("l1hpsTauMediumIso", &l1hpsTauMediumIso_);
  tree_ -> Branch("l1hpsTauLooseIso", &l1hpsTauLooseIso_);
  tree_ -> Branch("l1hpsTauVLooseIso", &l1hpsTauVLooseIso_);
  tree_ -> Branch("l1hpsTauTightRelIso", &l1hpsTauTightRelIso_);
  tree_ -> Branch("l1hpsTauMediumRelIso", &l1hpsTauMediumRelIso_);
  tree_ -> Branch("l1hpsTauLooseRelIso", &l1hpsTauLooseRelIso_);
  tree_ -> Branch("l1hpsTauVLooseRelIso", &l1hpsTauVLooseRelIso_);
  tree_ -> Branch("l1hpsTauZ", &l1hpsTauZ_);
  tree_ -> Branch("l1hpsTauLeadTrackPt",  &l1hpsTauLeadTrackPt_);
  tree_ -> Branch("l1hpsTauType", &l1hpsTauType_);
  tree_ -> Branch("l1nnTauPt",  &l1nnTauPt_);
  tree_ -> Branch("l1nnTauEta", &l1nnTauEta_);
  tree_ -> Branch("l1nnTauPhi", &l1nnTauPhi_);
  tree_ -> Branch("l1nnTauCharge", &l1nnTauCharge_);
  tree_ -> Branch("l1nnTauType", &l1nnTauType_);
  tree_ -> Branch("l1nnTauIso", &l1nnTauIso_);
  tree_ -> Branch("l1nnTauTightIso", &l1nnTauTightIso_);
  tree_ -> Branch("l1nnTauMediumIso", &l1nnTauMediumIso_);
  tree_ -> Branch("l1nnTauLooseIso", &l1nnTauLooseIso_);
  tree_ -> Branch("l1nnTauVLooseIso", &l1nnTauVLooseIso_);
  tree_ -> Branch("l1nnTauTightRelIso", &l1nnTauTightRelIso_);
  tree_ -> Branch("l1nnTauMediumRelIso", &l1nnTauMediumRelIso_);
  tree_ -> Branch("l1nnTauLooseRelIso", &l1nnTauLooseRelIso_);
  tree_ -> Branch("l1nnTauVLooseRelIso", &l1nnTauVLooseRelIso_);
  tree_ -> Branch("l1nnTauZ", &l1nnTauZ_);
  tree_ -> Branch("l1nnTauLeadTrackPt",  &l1nnTauLeadTrackPt_);
  tree_ -> Branch("l1nnTauType", &l1nnTauType_);

  hist_genTauPt_ = new TH1F("genTauPt","genTauPt", 100, 0., 1000.);
  hist_genTauEta_ = new TH1F("genTauEta","genTauEta",50, -3., 3.);
  hist_genTauPhi_ = new TH1F("genTauPhi","genTauPhi",50, -3., 3.);
  hist_isGenMatched_ = new TH1F("isGenMatched","isGenMatched", 3, -1., 2.);
  hist_recoTauPt_ = new TH1F("recoTauPt","recoTauPt", 100, 0., 1000.);
  hist_recoTauEta_ = new TH1F("recoTauEta","recoTauEta",50, -3., 3.);
  hist_recoTauPhi_ = new TH1F("recoTauPhi","recoTauPhi",50, -3., 3.);
  hist_isRecoMatched_ = new TH1F("isRecoMatched","isRecoMatched", 3, -1., 2.);
  hist_recoGMTauPt_ = new TH1F("recoGMTauPt","recoGMTauPt", 100, 0., 1000.);
  hist_recoGMTauEta_ = new TH1F("recoGMTauEta","recoGMTauEta",50, -3., 3.);
  hist_recoGMTauPhi_ = new TH1F("recoGMTauPhi","recoGMTauPhi",50, -3., 3.);
  hist_isRecoGMMatched_ = new TH1F("isRecoGMMatched","isRecoGMMatched", 3, -1., 2.);
  hist_hltTauPt_ = new TH1F("hltTauPt","hltTauPt", 100, 0., 1000.);
  hist_hltTauEta_ = new TH1F("hltTauEta","hltTauEta",50, -3., 3.);
  hist_hltTauPhi_ = new TH1F("hltTauPhi","hltTauPhi",50, -3., 3.);
  hist_hltTauReso_vs_Gen_ = new TH1F("hltTauReso_vs_Gen","hltTauReso_vs_Gen", 60, 0., 3.);
  hist_hltTauReso_vs_Reco_ = new TH1F("hltTauReso_vs_Reco","hltTauReso_vs_Reco", 60, 0., 3.);
  hist_hltTauReso_vs_RecoGM_ = new TH1F("hltTauReso_vs_RecoGM","hltTauReso_vs_RecoGM", 60, 0., 3.);


  return;
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1andHLTTauAnalyzer::endJob()
{
  if(createHistRoorFile_){
    histRootFile_->cd();
    
    hist_genTauPt_->Write();
    hist_genTauEta_->Write();
    hist_genTauPhi_->Write();
    hist_isGenMatched_->Write();
    hist_recoTauPt_->Write();
    hist_recoTauEta_->Write();
    hist_recoTauPhi_->Write();
    hist_isRecoMatched_->Write();
    hist_recoGMTauPt_->Write();
    hist_recoGMTauEta_->Write();
    hist_recoGMTauPhi_->Write();
    hist_isRecoGMMatched_->Write();
    hist_hltTauPt_->Write();  
    hist_hltTauEta_->Write();
    hist_hltTauPhi_->Write();
    hist_hltTauReso_vs_Gen_->Write();
    hist_hltTauReso_vs_Reco_->Write();
    hist_hltTauReso_vs_RecoGM_->Write();
    
  //  histRootFile_->Write();
    histRootFile_->Close();
  }
  

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1andHLTTauAnalyzer);
