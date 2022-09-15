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
#include <TH2F.h>
#include <TProfile2D.h>

#include "DataFormats/VertexReco/interface/Vertex.h"  
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <DataFormats/PatCandidates/interface/Tau.h> 

#include "DataFormats/TauReco/interface/PFTau.h"                // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"             // reco::PFTauRef, reco::PFTauCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"   // reco::PFTauDiscriminator

#include "DataFormats/L1TParticleFlow/interface/HPSPFTau.h"     // l1t::L1HPSPFTau
#include "DataFormats/L1TParticleFlow/interface/HPSPFTauFwd.h"  // l1t::HPSPFTauCollection
#include "DataFormats/L1TParticleFlow/interface/PFTau.h"       // l1t::PFTauCollection for NNTau

namespace {
  int
  constrainValue(int value,
                 int lowerBound,
                 int upperBound)
  {
    assert(lowerBound <= upperBound);
    value = std::max(value, lowerBound);
    value = std::min(value, upperBound);
    return value;
  }

  double
  square(double x)
  {
    return x*x;
  }

  void
  fillWOverflow(TH1F* histogram, 
                double x)
  {
    if(!histogram) return;
    if ( !histogram->GetSumw2N() ) histogram->Sumw2();
    const TAxis * const xAxis = histogram->GetXaxis();
    const int bin = constrainValue(xAxis->FindBin(x), 1, xAxis->GetNbins());
    const double binContent = histogram->GetBinContent(bin);
    const double binError   = histogram->GetBinError(bin);
    histogram->SetBinContent(bin, binContent+1);
    const double newerror = TMath::Sqrt( square(binError) + square(1) );
    histogram->SetBinError(bin, newerror);
  }
  
}

class HLTTauIsolationStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  enum ParticleType {
    X = 0,     // undefined
    h,         // charged hadron
    e,         // electron
    mu,        // muon
    gamma,     // photon
    h0,        // neutral hadron
    h_HF,      // HF tower identified as a hadron
    egamma_HF  // HF tower identified as an EM particle
  };
  enum GenMatch {
    kNoGenMatch = 0,
    kGenMatch
  };

  
  explicit HLTTauIsolationStudy(const edm::ParameterSet&);
  ~HLTTauIsolationStudy();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  TH1F* h_analyzed;

  // ----------member data ---------------------------
  bool debug_;
  edm::EDGetTokenT<std::vector<reco::Vertex>>      recoVertexToken_;
  edm::EDGetTokenT<pat::TauCollection>      hltTauToken_;
  std::string hltTauSumChargedIsoToken_;
  edm::EDGetTokenT<std::vector<reco::GenJet>>      genTauToken_;
  TH1F* h_iso;
  TH1F* h_deno;
  TH1F* h_neo;
};

HLTTauIsolationStudy::HLTTauIsolationStudy(const edm::ParameterSet& iConfig)
  : debug_          (iConfig.getUntrackedParameter<bool>("debug", false))
  ,  hltTauToken_    (consumes<pat::TauCollection>           (iConfig.getParameter<edm::InputTag>("hltTauToken")))
  , hltTauSumChargedIsoToken_ (iConfig.getParameter<std::string>("hltTauSumChargedIsoToken"))
  , genTauToken_ (consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genTauToken")))
{
   //now do what ever initialization is needed
  return;
}


HLTTauIsolationStudy::~HLTTauIsolationStudy()
{}

// member functions
//

// ------------ method called for each event  ------------
void
HLTTauIsolationStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(debug_){
    std::cout<<" Starting L1 nad HLTTau Analyzer ............     "<< std::endl;  
  }
  using namespace edm;
  
  edm::Handle<pat::TauCollection>            hltTauHandle;
  iEvent.getByToken(hltTauToken_,            hltTauHandle);

  //edm::Handle<std::vector<reco::GenJet>>            gentauHandle;
  //iEvent.getByToken(genTauToken_,            gentauHandle);

  size_t numHLTTaus = hltTauHandle->size();
  if ( numHLTTaus ) h_analyzed->Fill(1);

  std::vector<reco::GenJet *> matchedjets;

  for (size_t idxHLTTau=0; idxHLTTau<numHLTTaus; ++idxHLTTau){
    pat::TauRef hltTau(hltTauHandle, idxHLTTau);
    assert(hltTau->pt()>20);
    assert(std::fabs(hltTau->eta())<2.4);
    assert(hltTau->tauID("decayModeFindingNewDMs") == 0 || hltTau->tauID("decayModeFindingNewDMs") ==1 || hltTau->tauID("decayModeFindingNewDMs") ==2 || hltTau->tauID("decayModeFindingNewDMs") ==10 || hltTau->tauID("decayModeFindingNewDMs") ==11);
    double reliso = hltTau->tauID(hltTauSumChargedIsoToken_);// / hltTau->pt();
    fillWOverflow(h_iso, reliso);
    /*if (hltTau->tauID("chargedIsoPtSumdR03")/hltTau->pt() <0.20)
    {
      double mindr(99);
      reco::GenJet * matchedtau = 0;
      for(auto &gentau : *gentauHandle)
      {
        if ( !(gentau.pt() > 20 && fabs(gentau.eta()) < 2.3) ) continue;
        double deltar = deltaR(hltTau->eta(), hltTau->phi(), gentau.eta(), gentau.phi());
        if ( deltar <0.3 && deltar < mindr )
        {
          mindr = deltar;
          matchedtau = const_cast<reco::GenJet*>(&gentau);
        }
      }
      if (mindr != 99 && std::count(matchedjets.begin(), matchedjets.end(), matchedtau) == 0)
      {
        fillWOverflow(h_neo, matchedtau->pt());
        matchedjets.push_back(matchedtau);
      }
      }*/
  }
  /*for(const auto gentau : *gentauHandle)
  {
    if (gentau.pt() > 20 && fabs(gentau.eta()) < 2.3)
      fillWOverflow(h_deno, gentau.pt());
      }*/
}

void
HLTTauIsolationStudy::beginJob()
{
  edm::Service<TFileService> fs;
  h_analyzed = fs->make<TH1F>("analyzed", "analyzed", 1, 0.5, 1.5);
  h_analyzed->Sumw2();
  h_iso = fs->make<TH1F>("chargedIsoPtSumdR03", "chargedIsoPtSumdR03/pt", 400, 0, 100.);//2.);
  h_iso->Sumw2();
  h_deno = fs->make<TH1F>("deno", "deno", 20, 20, 100.);
  h_deno->Sumw2();
  h_neo = fs->make<TH1F>("neo", "neo", 20, 20, 100.);
  h_neo->Sumw2();
}
// ------------ method called once each job just after ending the event loop  ------------
void
HLTTauIsolationStudy::endJob()
{}
//define this as a plug-in
DEFINE_FWK_MODULE(HLTTauIsolationStudy);
