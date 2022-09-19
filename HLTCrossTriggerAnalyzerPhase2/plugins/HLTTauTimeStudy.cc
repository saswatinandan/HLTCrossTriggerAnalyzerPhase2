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
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include <TTree.h>
#include <TH2F.h>
#include <TProfile2D.h>

#include "DataFormats/VertexReco/interface/Vertex.h"  
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <DataFormats/PatCandidates/interface/Tau.h> 

#include "DataFormats/TauReco/interface/PFTau.h"                // reco::PFTau
#include "DataFormats/TauReco/interface/PFTauFwd.h"             // reco::PFTauRef, reco::PFTauCollection
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"   // reco::PFTauDiscriminator

#include "DataFormats/L1TParticleFlow/interface/HPSPFTau.h"     // l1t::L1HPSPFTau
#include "DataFormats/L1TParticleFlow/interface/HPSPFTauFwd.h"  // l1t::HPSPFTauCollection
#include "DataFormats/L1TParticleFlow/interface/PFTau.h"       // l1t::PFTauCollection for NNTau
#include "DataFormats/JetReco/interface/PFJet.h"

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
  float
  printpt(const pat::Tau& tau, std::string type)
  {
    float taupt(0);
    reco::CandidatePtrVector consts;
    if (type == "isolation" )
      consts = tau.isolationCands();
    else if (type == "signal" )
      consts = tau.signalCands();
    for (auto cand: consts) taupt+= cand->pt();
    return taupt;
  }
  bool
  check_tauproperty(const pat::Tau& recoPatTau)
  {
    return (recoPatTau.pt()>20 && std::fabs(recoPatTau.eta())<2.4 && (recoPatTau.tauID("decayModeFindingNewDMs") == 0 || recoPatTau.tauID("decayModeFindingNewDMs") ==1 || recoPatTau.tauID("decayModeFindingNewDMs") ==2 || recoPatTau.tauID("decayModeFindingNewDMs") ==10 || recoPatTau.tauID("decayModeFindingNewDMs") ==11));
  }
  /*  bool
  check_tauproperty(const reco::PFTau& recoPatTau)
  {
    return (recoPatTau.pt()>20 && std::fabs(recoPatTau.eta())<2.4 && (recoPatTau.tauID("decayModeFindingNewDMs") == 0 || recoPatTau.tauID("decayModeFindingNewDMs") ==1 || recoPatTau.tauID("decayModeFindingNewDMs") ==2 || recoPatTau.tauID("decayModeFindingNewDMs") ==10 || recoPatTau.tauID("decayModeFindingNewDMs") ==11));
    }*/
  template<typename recpartType>
  bool
  do_genmatch(const reco::GenParticle& genpart, const recpartType& recopart)
  {
    double genpt = genpart.pt();

    if ( deltaR(recopart.eta(), recopart.phi(), genpart.eta(), genpart.phi()) < 0.05 )
    {
      if ( fabs(recopart.pt() - genpt) / genpt < 0.50 );
      return true;
    }
    return false;
  }
}
class HLTTauTimeStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
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

  
  explicit HLTTauTimeStudy(const edm::ParameterSet&);
  ~HLTTauTimeStudy();
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  std::map<int, std::map<int, std::map<std::string, TH1F*> >> genMatched_histogram;
  TProfile2D* hprofile;
  TH1F* h_analyzed;
  TH1F* h_vertex;
  TH2F* h_tau_rec_vertex;
  TH1F* h_genmatched_timeinfo;
  TH1F* h_nogenmatched_timeinfo;
  TH1F* h_eta_wtime;
  TH1F* h_pt_wtime;
  TH1F* h_eta_wotime;
  TH1F* h_pt_wotime;
  TH1F* h_genmatchedavgtime;
  TH1F* h_vertexmatchedavgtime;
  TH2F* h_avgtime;
  TH1F* h_deno;
  TH1F* h_neo;
  TH1F* h_taupt;
  TH1F* h_jetpt;
  TH1F* h_jetneutralsumpt;
  TH1F* h_jetchargedpt;
  TH1F* h_pdgid;
  TH1F* h_pionpt;
  TH1F* h_kaonpt;
  TH1F* h_gamapt;
  TH1F* h_elept;
  TH1F* h_muonpt;
  TH1F* h_jettaupt;
  TH1F* h_jetchargedsumpt;
  TH1F* h_tausignal;
  TH1F* h_tauisolation;
  TH1F * h_diff;
  TH1F * h_diff_inclu;
  TH1F * h_diff_high;
  TH1F* h_tauptnocut;
  TH1F* h_taueta;
  TH1F* h_tauiso;
  TH1F* h_taudecay;
  TH1F* h_match;
  TH1F* h_jeteta;
  TH1F* h_iso;
  TH1F* h_reliso;
  std::map<int, TFileDirectory> gdir;
  ULong64_t       indexevents_;
  Int_t           runNumber_;
  Int_t           lumi_;
  double genVertex_;
  std::map<TDirectory *, std::vector<TH1 *>> gHistograms_;

  // ----------member data ---------------------------
  bool debug_;
  bool usegenJets_;
  edm::EDGetTokenT<std::vector<reco::Vertex>>      recoVertexToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection>      pfCandToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<vector<pat::Tau>>      recoPatTauToken_;
  edm::EDGetTokenT<std::vector<reco::PFJet>> recoJetToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> genvertexToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puSummaryToken_;
  edm::EDGetTokenT<std::vector<reco::PFTau>> recoPFTauToken_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> pfTauDiscrimination_byDecayModeToken_;
};
  
HLTTauTimeStudy::HLTTauTimeStudy(const edm::ParameterSet& iConfig)
  : debug_          (iConfig.getUntrackedParameter<bool>("debug", false))
  , usegenJets_ (iConfig.getParameter<bool>("usegenJets"))
  , recoVertexToken_(consumes<std::vector<reco::Vertex>>    (iConfig.getParameter<edm::InputTag>("recoVertexToken")))
  ,  pfCandToken_ (consumes<reco::PFCandidateCollection>    (iConfig.getParameter<edm::InputTag>("recoParticleToken")))
  , genParticleToken_ (consumes<reco::GenParticleCollection>    (iConfig.getParameter<edm::InputTag>("genParticleToken")))
  , recoPatTauToken_ (consumes<vector<pat::Tau>>    (iConfig.getParameter<edm::InputTag>("recoPatTauToken")))
  , recoJetToken_ (consumes<vector<reco::PFJet>>    (iConfig.getParameter<edm::InputTag>("genJetToken")))
  , genvertexToken_ (consumes<vector<SimVertex>>    (iConfig.getParameter<edm::InputTag>("genvertexToken")))
  , puSummaryToken_ (consumes<vector<PileupSummaryInfo>>    (iConfig.getParameter<edm::InputTag>("puSummaryToken")))
  , recoPFTauToken_ (consumes<vector<reco::PFTau>>    (iConfig.getParameter<edm::InputTag>("recoPFTauToken")))
  , pfTauDiscrimination_byDecayModeToken_ (consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("taudecayModediscriminationToken")))
{
  return;
}


HLTTauTimeStudy::~HLTTauTimeStudy()
{
  for ( int particle_id=ParticleType::X; particle_id<=ParticleType::egamma_HF; particle_id++ )
  {
    for ( int match=GenMatch::kNoGenMatch; match<=GenMatch::kGenMatch; match++)
    {
      for ( std::map<std::string, TH1F*>::const_iterator itr = genMatched_histogram[match][particle_id].begin(); itr != genMatched_histogram[match][particle_id].end(); itr++)
      {
        delete itr->second;
      }
    }
  }
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
  void
HLTTauTimeStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(debug_){
    std::cout<<" Starting L1 nad HLTTau Analyzer ............     "<< std::endl;  
  }
  using namespace edm;
  h_analyzed->Fill(1);
  
  edm::Handle<reco::VertexCollection>  recvertexHandle;
  iEvent.getByToken(recoVertexToken_, recvertexHandle);
  h_vertex->Fill(recvertexHandle->size());

  edm::Handle<reco::PFCandidateCollection>  recpPFCandHandle;
  iEvent.getByToken(pfCandToken_, recpPFCandHandle);
  
  edm::Handle<reco::GenParticleCollection>  recogenPartCandHandle;
  iEvent.getByToken(genParticleToken_, recogenPartCandHandle);

  edm::Handle<vector<pat::Tau>>  recoPatTauHandle;
  iEvent.getByToken(recoPatTauToken_, recoPatTauHandle);

  edm::Handle<vector<reco::PFTau>>  recoPFTauHandle;
  iEvent.getByToken(recoPFTauToken_, recoPFTauHandle);

  edm::Handle<vector<reco::PFJet>>  recoJetHandle;
  iEvent.getByToken(recoJetToken_, recoJetHandle);
  
  edm::Handle<reco::PFTauDiscriminator> pfTauDiscrimination_byDecayModeHandle;
  iEvent.getByToken(pfTauDiscrimination_byDecayModeToken_, pfTauDiscrimination_byDecayModeHandle);

  //  edm::Handle<vector<double>>  rhoHandle_;
  //iEvent.getByToken(recoJetToken_, rhoHandle);

  reco::GenParticle tau1gen, tau2gen, maxptgen;
  int ntau(0);
  float maxpt(0);
  math::XYZPoint tauvertex(0,0,0);
  bool taufound(false);
  for(auto genpart : *recogenPartCandHandle)
  {
    if( genpart.status() == 1 && genpart.pt() > maxpt && genpart.charge() !=0)
    {
      maxpt = genpart.pt();
      maxptgen = const_cast<reco::GenParticle&>(genpart);
    }
    if (fabs(genpart.pdgId()) ==15 && !taufound)
    {
      taufound = true;
      tauvertex = genpart.vertex();
    }
  }
  math::XYZPoint vertex = maxptgen.vertex();
  h_tau_rec_vertex->Fill(tauvertex.z(), vertex.z());
  //  std::cout << "pt= " << maxptgen.pt() << "\t" << maxptgen.vertex().nTracks() << std::endl;

  int genmatched_time_count(0), vertexmatched_time_count(0);
  double genmatched_sum_time(0), vertexmatched_sum_time(0);
  int nrecoPFTau(recoPFTauHandle->size());

  /*if (0) 
  {
    if ( !usegenJets_ )
    {
      for(const auto recoGenCand : *recogenPartCandHandle)
      {
        if (recoGenCand.pt() > 20 && fabs(recoGenCand.eta()) < 2.3 && fabs(recoGenCand.pdgId()) == 15)
          fillWOverflow(h_deno, recoGenCand.pt());
      }
      for(const auto recoPatTau : *recoPatTauHandle)
      {
        if (recoPatTau.pt() >20 && std::fabs(recoPatTau.eta()) <2.3 && recoPatTau.tauID("chargedIsoPtSumdR03")/recoPatTau.pt() <0.20 && (recoPatTau.tauID("decayModeFindingNewDMs") == 0 || recoPatTau.tauID("decayModeFindingNewDMs") == 1 || recoPatTau.tauID("decayModeFindingNewDMs") == 2 || recoPatTau.tauID("decayModeFindingNewDMs") == 10 || recoPatTau.tauID("decayModeFindingNewDMs") == 11) )
        {
          for(const auto recoGenCand : *recogenPartCandHandle)
          {
            if ( fabs(recoGenCand.pdgId()) != 15 ) continue;
            if ( deltaR(recoPatTau.eta(), recoPatTau.phi(), recoGenCand.eta(), recoGenCand.phi()) < 0.3)
            {
              fillWOverflow(h_neo, recoGenCand.pt());
              break;
            }
          }
        }
      }
    }
    else
    {
      //    std::cout << "new event " << std::endl;
      
      for(const auto recoJet : *recoJetHandle)
      {
        if(recoJet.pt() >20) fillWOverflow(h_jetpt, recoJet.pt());
        if (recoJet.pt() > 20 && fabs(recoJet.eta()) < 2.3 )
        {
          fillWOverflow(h_jeteta, recoJet.pt());
          std::vector<reco::PFCandidatePtr> jetcands = recoJet.getPFConstituents();
          float chargedpt = 0;
          float neutralpt=0;
          for (auto part: jetcands)
          {
            fillWOverflow(h_pdgid, fabs(part->pdgId()));
            if(fabs(part->pdgId()) == 211) fillWOverflow(h_pionpt, part->pt());
            if(fabs(part->pdgId()) == 130) fillWOverflow(h_kaonpt, part->pt());
            if(fabs(part->pdgId()) == 22) fillWOverflow(h_gamapt, part->pt());
            if(fabs(part->pdgId()) == 11) fillWOverflow(h_elept, part->pt());
            if(fabs(part->pdgId()) == 13) fillWOverflow(h_muonpt, part->pt());
            if(fabs(part->pdgId()) == 15) fillWOverflow(h_jettaupt, part->pt());
            if(part->charge())
            {
              fillWOverflow(h_jetchargedpt, part->pt());
              chargedpt += part->pt();
            }
            else
            {
              neutralpt += part->pt();
            }
          }
          fillWOverflow(h_jetchargedsumpt, chargedpt);
          fillWOverflow(h_jetneutralsumpt, neutralpt);
        }
        //      if(recoJet.pt() > 20 && recoJet.pt() < 150)std::cout << "jet pt: " << recoJet.pt() << "\tjet eta: " << recoJet.eta() << std::endl;
        if (recoJet.pt() > 20 && fabs(recoJet.eta()) < 2.3 )
          fillWOverflow(h_deno, recoJet.pt());
      }
      std::vector<reco::PFJet *> matchedjets;
      for(const auto recoPatTau : *recoPatTauHandle)
      {
        if(recoPatTau.pt() >20 && recoPatTau.pt() <100)
        {
          float taupt = printpt(recoPatTau, "signal");
          fillWOverflow(h_tausignal, taupt);
          //std::cout << "new tau " << std::endl;
          //std::cout << "signal pt sum: " << taupt << std::endl;
          taupt = printpt(recoPatTau, "isolation");
          //std::cout << "isolation pt sum: " << printpt(recoPatTau, "isolation") << std::endl;
          //std::cout << "reco tau pt: " << recoPatTau.pt() << std::endl;
          //std::cout << "reco tau eta: " << std::fabs(recoPatTau.eta()) << "\tcharged rel iso: " << recoPatTau.tauID("chargedIsoPtSumdR03")/recoPatTau.pt() << "\tdecay mode: " << recoPatTau.tauID("decayModeFindingNewDMs") << std::endl;
        }
        if (recoPatTau.pt() >20 && std::fabs(recoPatTau.eta()) <2.3  && recoPatTau.tauID("chargedIsoPtSumdR03")/recoPatTau.pt() <0.20 && (recoPatTau.tauID("decayModeFindingNewDMs") == 0 || recoPatTau.tauID("decayModeFindingNewDMs") == 1 || recoPatTau.tauID("decayModeFindingNewDMs") == 2 || recoPatTau.tauID("decayModeFindingNewDMs") == 10 || recoPatTau.tauID("decayModeFindingNewDMs") == 11) )
        {
          double tau_total = printpt(recoPatTau, "isolation")+printpt(recoPatTau, "signal");
          //std::cout << "entering" << std::endl;
          double mindr = 99.;
          reco::PFJet * matchedtau = 0;
          for(const auto &recoJet : *recoJetHandle)
          {
            if ( !(recoJet.pt() > 20 && fabs(recoJet.eta()) < 2.3) ) continue;
            double deltar =  deltaR(recoPatTau.eta(), recoPatTau.phi(), recoJet.eta(), recoJet.phi());
            //          if (tau_total >20 && tau_total<150)std::cout << "deltar: " << deltar << "\tjet eta, jet phi " << recoJet.eta() << "\t" <<  recoJet.phi() << "\t" << recoJet.pt() << std::endl;
            if ( deltar<0.3  && deltar<mindr)
            {
              mindr = deltar;
              matchedtau = const_cast<reco::PFJet*>(&recoJet);
            }
        }
          if (mindr != 99. && std::count(matchedjets.begin(), matchedjets.end(), matchedtau) == 0)
          {
            fillWOverflow(h_neo, matchedtau->pt());
            matchedjets.push_back(matchedtau);
            fillWOverflow(h_match, matchedtau->pt());
            //std::cout << "matched jet pt: " << matchedtau->pt() << std::endl;
            if(tau_total > 20 && tau_total <100)std::cout << "matched jet pt: "   << matchedtau->pt() << "\t tau iso+signal pt / matched jet pt " << tau_total/matchedtau->pt() << std::endl; 
            if ( tau_total > 20 && tau_total < 150)fillWOverflow(h_diff, (printpt(recoPatTau, "isolation")+printpt(recoPatTau, "signal"))/matchedtau->pt());                               
            fillWOverflow(h_diff_inclu, (printpt(recoPatTau, "isolation")+printpt(recoPatTau, "signal"))/matchedtau->pt());
            if (tau_total >=150)fillWOverflow(h_diff_high, (printpt(recoPatTau, "isolation")+printpt(recoPatTau, "signal"))/matchedtau->pt());
          }
        }
        /*     fillWOverflow(h_neo, recoJet.pt());
               if(recoPatTau.pt() > 20 && recoPatTau.pt() <100)std::cout << "matched jet pt: " << recoJet.pt() << "\t tau iso+signal pt / matched jet pt " << (printpt(recoPatTau, "isolation")+printpt(recoPatTau, "signal"))/recoJet.pt() << std::endl;
               if ( tau_total > 20 && tau_total < 150)fillWOverflow(h_diff, (printpt(recoPatTau, "isolation")+printpt(recoPatTau, "signal"))/recoJet.pt());
               fillWOverflow(h_diff_inclu, (printpt(recoPatTau, "isolation")+printpt(recoPatTau, "signal"))/recoJet.pt());
               if (tau_total >=150)fillWOverflow(h_diff_high, (printpt(recoPatTau, "isolation")+printpt(recoPatTau, "signal"))/recoJet.pt());
               break;
               }
               }
               }
               double mindr = 99.;
               reco::PFJet * matchedtau = 0;
               std::vector<reco::PFJet *> matchedjets1;
               for(const auto &recoJet : *recoJetHandle)
               {
               if ( !(recoJet.pt() > 20 && fabs(recoJet.eta()) < 2.3) ) continue;
               double deltar =  deltaR(recoPatTau.eta(), recoPatTau.phi(), recoJet.eta(), recoJet.phi());
               if ( deltar<0.3 ){/* && deltar<mindr)
               {
               mindr = deltar;
               matchedtau = const_cast<reco::PFJet*>(&recoJet);
               }
               }
               if (mindr != 99. && std::count(matchedjets1.begin(), matchedjets1.end(), matchedtau) == 0)
               {
               if(recoPatTau.pt() >20) fillWOverflow(h_taupt, recoJet.pt());
               if(recoPatTau.pt() >20 && std::fabs(recoPatTau.eta()) <2.3) fillWOverflow(h_taueta, recoJet.pt());
               if(recoPatTau.pt() >20 && std::fabs(recoPatTau.eta()) <2.3 && recoPatTau.tauID("chargedIsoPtSumdR03")/recoPatTau.pt() <0.20 ) fillWOverflow(h_tauiso, recoJet.pt());
               if(recoPatTau.pt() >20 && std::fabs(recoPatTau.eta()) <2.3 && recoPatTau.tauID("chargedIsoPtSumdR03")/recoPatTau.pt() <0.20 && (recoPatTau.tauID("decayModeFindingNewDMs") ==  0 || recoPatTau.tauID("decayModeFindingNewDMs") == 1 || recoPatTau.tauID("decayModeFindingNewDMs") == 2 || recoPatTau.tauID("decayModeFindingNewDMs") == 10 || recoPatTau.tauID("decayModeFindingNewDMs") == 11)) fillWOverflow(h_taudecay, recoJet.pt());
               break;
               //        matchedjets1.push_back(matchedtau);
               ############}
      }
    }
  }*/

  for(const auto recoPFCand : *recpPFCandHandle) 
  {
    
    bool genmatched(0);
    double recpt = recoPFCand.pt();
    double receta = recoPFCand.eta();
    int recparticleid = recoPFCand.particleId();

    for(const auto recoGenCand : *recogenPartCandHandle)
    {
      genmatched = do_genmatch<reco::PFCandidate>(recoGenCand, recoPFCand);
      if (genmatched) break;
    }
    double time = recoPFCand.time();
    bool hastimeinfo = recoPFCand.isTimeValid();
    if ( hastimeinfo )
    {
      fillWOverflow(h_eta_wtime,  abs(receta));
      fillWOverflow(h_pt_wtime, recpt);
    }
    else
    {
      fillWOverflow(h_eta_wotime, abs(receta));
      fillWOverflow(h_pt_wotime,  recpt);
    }
    if ( genmatched ) 
    {
      h_genmatched_timeinfo->Fill(hastimeinfo);
    }
    else
    {
      h_nogenmatched_timeinfo->Fill(hastimeinfo);
    }
    for ( int match=GenMatch::kNoGenMatch; match<=GenMatch::kGenMatch; match++)
    {
      if ( ! (genmatched == match) ) continue;
      for ( int particle_id=ParticleType::X; particle_id<=ParticleType::egamma_HF; particle_id++)
      {
        if ( recparticleid != particle_id ) continue;
        if ( hastimeinfo )
        {
          fillWOverflow(genMatched_histogram[match][particle_id]["time_distribution"], time);
          fillWOverflow(genMatched_histogram[match][particle_id]["pt_distribution"], recpt);
          if ( recparticleid == ParticleType::h && recpt > 2 && abs(receta) < 2.4 )
          {
            if ( genmatched )
            {
              genmatched_time_count += 1;
              genmatched_sum_time += time;
            }
            math::XYZPoint recvertex = recoPFCand.vertex();
            math::XYZPoint diff(recvertex.x() - vertex.x(), recvertex.y() - vertex.y(), recvertex.z() - vertex.z());
            float dxy = TMath::Sqrt(square(diff.x()) + square(diff.y()));
            if ( fabs(diff.z()) < 0.1 && dxy < 0.05 )
            {
              vertexmatched_time_count += recpt;
              vertexmatched_sum_time += (recpt * time);
            }
          }
        }
        if ( recpt > 2 && abs(receta) < 2.4 )fillWOverflow(genMatched_histogram[match][particle_id]["time_info"], hastimeinfo);
      }
    }
  }
  float genmatched_avgtime = ( genmatched_time_count ) ? genmatched_sum_time / genmatched_time_count : -100;
  if ( genmatched_avgtime > -100 ) fillWOverflow(h_genmatchedavgtime, genmatched_avgtime );
  float vertexmatched_avgtime = ( vertexmatched_time_count ) ? vertexmatched_sum_time / vertexmatched_time_count : -100;
  if ( vertexmatched_avgtime > -100 ) fillWOverflow(h_vertexmatchedavgtime, vertexmatched_avgtime );
  h_avgtime->Fill(genmatched_avgtime, vertexmatched_avgtime);
  for(const auto recoPFCand : *recpPFCandHandle)
  {
    
    bool genmatched(0);

    for(const auto recoGenCand : *recogenPartCandHandle)
    {
      genmatched = do_genmatch<reco::PFCandidate>(recoGenCand, recoPFCand);
      if (genmatched) break;
    }
    int recparticleid = recoPFCand.particleId();
    for ( int match=GenMatch::kNoGenMatch; match<=GenMatch::kGenMatch; match++)
    {
      if ( ! (genmatched == match) ) continue;
      for ( int particle_id=ParticleType::X; particle_id<=ParticleType::egamma_HF; particle_id++)
      {
        if ( particle_id != ParticleType::h ) continue;
        if ( recparticleid != particle_id ) continue;
        if ( recoPFCand.isTimeValid() )
        {
          if ( genmatched_avgtime > -100 )
          {
            if ( recoPFCand.pt() > 2 && abs(recoPFCand.eta()) < 2.4 && genmatched )
              hprofile->Fill(recoPFCand.eta(), recoPFCand.phi(), recoPFCand.time() - genmatched_avgtime);
            fillWOverflow(genMatched_histogram[match][particle_id]["time_shift"], recoPFCand.time() - genmatched_avgtime );
            if ( genmatched_time_count >=3 && genmatched ) fillWOverflow(genMatched_histogram[match][particle_id]["time_shift_with5"], recoPFCand.time() - genmatched_avgtime );
          }
          if ( vertexmatched_avgtime > -100 )
          {
            math::XYZPoint diff(recoPFCand.vertex().x() - vertex.x(), recoPFCand.vertex().y() - vertex.y(), recoPFCand.vertex().z() - vertex.z());
            bool dxy = TMath::Sqrt(square(diff.x()) + square(diff.y())) < 0.05;
            bool dz = fabs(diff.z()) < 0.1;
            fillWOverflow(genMatched_histogram[dz][particle_id]["time_shift_vertexmatched"], recoPFCand.time() - vertexmatched_avgtime );
          }
        }
      }
    }
  }

  for (size_t tau_idx=0; tau_idx<recoPFTauHandle->size(); ++tau_idx)
    {
      auto tau = (*recoPFTauHandle)[tau_idx];
      if (!(tau.pt() >20 && fabs(tau.eta())<2.4 && (tau.decayMode()==0|| tau.decayMode()==1|| tau.decayMode()==2|| tau.decayMode()==10|| tau.decayMode()==11)))continue;
      //std::cout << tau.pt() << "\t" << tau.decayMode() << std::endl;
      for(auto sig : tau.signalPFCands()){
        //      if(sig->pdgId() !=22)std::cout << "time" << sig->time() << "\t" << sig->timeError() << std::endl;
       }
    }
  for (auto tau : *recoPatTauHandle) {
    int cand(0);
    float leadtime(0);
    bool leadgenmatched(0);
    reco::CandidatePtrVector signal = tau.signalCands();
    //std::cout << "new " << tau.pt() << std::endl;
    for (auto sig : signal){
      if (sig->pdgId() == 22) continue;
      bool genmatch(false);
      for (auto genpart : *recogenPartCandHandle)
      {
        genmatch = do_genmatch<reco::Candidate>(genpart, *sig);
        if (genmatch) break;
      }
      cand +=1;
      auto packedCand = dynamic_cast<const pat::PackedCandidate*>(sig.get());
      math::XYZPoint diff(packedCand->vertex().x() - vertex.x(), packedCand->vertex().y() - vertex.y(), packedCand->vertex().z() - vertex.z());
      bool dxy = TMath::Sqrt(square(diff.x()) + square(diff.y())) < 0.05;
      bool dz = fabs(diff.z()) < 0.1;
      //std::cout << "time " << packedCand->time() << "\t" << packedCand->vertexRef()->t() << "\t" << packedCand->vertexRef()->position().x() << "\t" << packedCand->vertexRef()->position().y() << "\t" << packedCand->vertexRef()->position().z() << "\t" << packedCand->vertexRef()->nTracks() << "\t" << genmatch << std::endl;
      if(cand ==1) {
        leadtime = packedCand->time();
        leadgenmatched = genmatch;
        fillWOverflow(genMatched_histogram[genmatch][ParticleType::h]["tauleadcandidate"], 1);
        fillWOverflow(genMatched_histogram[genmatch][ParticleType::h]["time_shift_leadcandidate"], leadtime-genmatched_avgtime);
        fillWOverflow(genMatched_histogram[dz][ParticleType::h]["time_shift_leadcandidate_vertexmatched"], leadtime-vertexmatched_avgtime);
        continue;
      }
      else if (cand==2) {
        fillWOverflow(genMatched_histogram[genmatch && leadgenmatched][ParticleType::h]["time_shift_of_1stsubleadCand_wrt_lead_chargedCand"], leadtime - packedCand->time());
        fillWOverflow(genMatched_histogram[genmatch][ParticleType::h]["time_shift_subleadcandidate"], packedCand->time()-genmatched_avgtime);
        fillWOverflow(genMatched_histogram[dz][ParticleType::h]["time_shift_subleadcandidate_vertexmatched"], packedCand->time()-vertexmatched_avgtime);
        continue;
      }
      else if (cand==3) {
        fillWOverflow(genMatched_histogram[genmatch && leadgenmatched][ParticleType::h]["time_shift_of_2ndsubleadCand_wrt_lead_chargedCand"], leadtime - packedCand->time());
        fillWOverflow(genMatched_histogram[genmatch][ParticleType::h]["time_shift_2ndsubleadcandidate"], packedCand->time()-genmatched_avgtime);
        fillWOverflow(genMatched_histogram[dz][ParticleType::h]["time_shift_2ndsubleadcandidate_vertexmatched"], packedCand->time()-vertexmatched_avgtime);
        break;
      }
    }
  }
}

void
HLTTauTimeStudy::beginJob()
{
  std::map<int, string> map_particlename = {
    {ParticleType::X, "undefined"},
    {ParticleType::h, "charged_hadron"},
    {ParticleType::e, "electron"},
    {ParticleType::mu, "muon"},
    {ParticleType::gamma, "photon"},
    {ParticleType::h0, "neutral_hadron"},
    {ParticleType::h_HF, "HF_hadron"},
    {ParticleType::egamma_HF, "HF_EM"},
  };
  edm::Service<TFileService> fs;
  h_analyzed = fs->make<TH1F>("analyzed", "analyzed", 1, 0.5, 1.5);
  h_analyzed->Sumw2();
  h_deno = fs->make<TH1F>("deno", "deno", 50, 20, 300.);
  h_deno->Sumw2();
  h_neo = fs->make<TH1F>("neo", "neo", 50, 20, 300.);
  h_neo->Sumw2();
  h_taupt = fs->make<TH1F>("taupt", "taupt", 50, 20, 300);
  h_taupt->Sumw2();
  h_tauptnocut = fs->make<TH1F>("tauptnocut", "taupt no cut", 50, 0, 300);
  h_tauptnocut->Sumw2();
  h_tausignal = fs->make<TH1F>("tausignal", "taupt with signal candidate, no cut", 50, 0, 300);
  h_tausignal->Sumw2();
  h_tauisolation = fs->make<TH1F>("tauisolation", "taupt with isolation candidate, no cut", 50, 0, 300);
  h_tauisolation->Sumw2();
  h_diff = fs->make<TH1F>("diff", "ratio between matched jet pt to pt of tau signal and isolation candidate with jet pt<150", 50, 0.20, 1.50);
  h_diff->Sumw2();
  h_diff_inclu = fs->make<TH1F>("diff_inclusive", "ratio between matched jet pt to pt of tau signal and isolation candidate", 50, 0.20, 1.50);
  h_diff_inclu->Sumw2();
  h_diff_high = fs->make<TH1F>("diff_high", "ratio between matched jet pt to pt of tau signal and isolation candidate with jet pt>150", 50, 0.20, 1.50);
  h_diff_high->Sumw2();
  h_jetpt = fs->make<TH1F>("jetpt", "jetpt", 50, 20, 300);
  h_jetpt->Sumw2();
  h_jetchargedpt = fs->make<TH1F>("jetchargedpt", "jet charged pt", 50, 0, 10);
  h_jetchargedpt->Sumw2();
  h_jetchargedsumpt = fs->make<TH1F>("jetchargedsumpt", "jet charged sum pt", 50, 0, 200);
  h_jetchargedsumpt->Sumw2();
  h_pdgid = fs->make<TH1F>("pdgid", "pdgid of jet constituent", 211, -0.5,  211.5);
  h_pdgid->Sumw2();
  h_pionpt = fs->make<TH1F>("pion_pt", "pion pt", 50, 0,  10);
  h_pionpt->Sumw2();
  h_kaonpt = fs->make<TH1F>("kaon_pt", "kaon pt", 100, 0,  50);
  h_kaonpt->Sumw2();
  h_gamapt = fs->make<TH1F>("gama_pt", "gama pt", 50, 0,  10);
  h_gamapt->Sumw2();
  h_elept = fs->make<TH1F>("ele_pt", "ele pt", 50, 0,  10);
  h_elept->Sumw2();
  h_muonpt = fs->make<TH1F>("muon_pt", "muon pt", 50, 0,  10);
  h_muonpt->Sumw2();
  h_jettaupt = fs->make<TH1F>("tau_pt", "tau pt", 50, 0,  10);
  h_jettaupt->Sumw2();
  h_jetneutralsumpt = fs->make<TH1F>("jetneutralsumpt", "jet neutral sum pt", 50, 0, 200);
  h_jetneutralsumpt->Sumw2();
  h_taueta = fs->make<TH1F>("taueta", "taueta", 50, 20, 300.);
  h_taueta->Sumw2();
  h_taudecay = fs->make<TH1F>("taudecay", "taudecay", 50, 20, 300.);
  h_taudecay->Sumw2();
  h_match = fs->make<TH1F>("match", "match", 50, 20, 300.);
  h_match->Sumw2();
  h_tauiso = fs->make<TH1F>("tauiso", "tauiso", 50, 20, 300.);
  h_tauiso->Sumw2();
  h_jeteta = fs->make<TH1F>("jeteta", "jeteta", 50, 20, 300.);
  h_jeteta->Sumw2();
  h_iso = fs->make<TH1F>("iso", "iso", 20, 0, 100.);
  h_iso->Sumw2();
  h_reliso = fs->make<TH1F>("reliso", "reliso", 20, 0, 100.);
  h_reliso->Sumw2();
  int largebin=101;
  int smallbin=1;
  float * a = new float[2*largebin+smallbin];
  float small_binw=(0.0025+0.00025)/(smallbin-1);
  float large_binw = (0.5-0.)/(largebin-1);
  for(int i=0; i<largebin; i++) a[i]=-0.5 + large_binw*i;
  int count =0;
  /*for(int i=largebin; i<largebin+smallbin;i++) {
    count++;
    a[i]= a[largebin-1]+small_binw*count;
    }*/
  //  a[largebin] = -0.0002/2;
  //a[largebin+1] = 0;
  //a[largebin+2] = 0.0002/2.;
  a[largebin] = 0.00001;
  count =0;
  for(int i=largebin+smallbin; i<largebin+largebin+smallbin; i++) {
    count++;
    a[i] = a[largebin+smallbin-1] +  large_binw*count;
  }
  for ( int match=GenMatch::kNoGenMatch; match<=GenMatch::kGenMatch; match++)
  {
    for ( int particle_id=ParticleType::X; particle_id<=ParticleType::egamma_HF; particle_id++ )
    {
      std::string prefix = Form("%s",
         (match==GenMatch::kGenMatch) ? "genmatched" : "nogenmatched");
      TFileDirectory dir = fs->mkdir(map_particlename[particle_id]);
      gdir[particle_id] = dir;
      std::string name = Form("%s_time_distribution", prefix.data());
      genMatched_histogram[match][particle_id]["time_distribution"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
      name = Form("%s_pt_distribution", prefix.data());
      genMatched_histogram[match][particle_id]["pt_distribution"] = dir.make<TH1F>(name.data(), name.data(), 100, 0, 20);
      name = Form("%s_time_shift", prefix.data());
      genMatched_histogram[match][particle_id]["time_shift"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
      if ( match==GenMatch::kGenMatch )
      {
        genMatched_histogram[match][particle_id]["time_shift_with5"] = dir.make<TH1F>((name+"_with5").data(), (name+"_with5").data(), 200, -0.5, 0.5);
      }
      if ( particle_id == ParticleType::h )
      {
        name = Form("%s_time_shift_vertexmatched", prefix.data());
        genMatched_histogram[match][particle_id]["time_shift_vertexmatched"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_leadcandidate", prefix.data());
        genMatched_histogram[match][particle_id]["tauleadcandidate"] = dir.make<TH1F>(name.data(), name.data(), 1, 0.5, 1.5);
        name = Form("%s_time_shift_of_1stsubleadCand_wrt_lead_chargedCand", prefix.data());
        genMatched_histogram[match][particle_id]["time_shift_of_1stsubleadCand_wrt_lead_chargedCand"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_of_2ndsubleadCand_wrt_lead_chargedCand", prefix.data());
        genMatched_histogram[match][particle_id]["time_shift_of_2ndsubleadCand_wrt_lead_chargedCand"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_leadcandidate", prefix.data());
        genMatched_histogram[match][ParticleType::h]["time_shift_leadcandidate"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_leadcandidate_vertexmatched", prefix.data());
        genMatched_histogram[match][ParticleType::h]["time_shift_leadcandidate_vertexmatched"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_1stleadcandidate", prefix.data());
        genMatched_histogram[match][ParticleType::h]["time_shift_subleadcandidate"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_1stleadcandidate_vertexmatched", prefix.data());
        genMatched_histogram[match][ParticleType::h]["time_shift_subleadcandidate_vertexmatched"]= dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_2ndleadcandidate", prefix.data());
        genMatched_histogram[match][ParticleType::h]["time_shift_2ndsubleadcandidate"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_2ndleadcandidate_vertexmatched", prefix.data());
        genMatched_histogram[match][ParticleType::h]["time_shift_2ndsubleadcandidate_vertexmatched"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
      }
      name = Form("%s_time_info", prefix.data());
      genMatched_histogram[match][particle_id]["time_info"] = dir.make<TH1F>(name.data(), name.data(), 2, -0.5, 1.5);
    }
  }
  delete [] a;
  gdir[-1] = fs->mkdir("inclusive");
  TFileDirectory dir = gdir[-1];
  h_genmatched_timeinfo = dir.make<TH1F>("genmatched_timeinfo", "timeinfo", 2, -0.5, 1.5);
  h_genmatchedavgtime = new TH1F("genmatchedavg_time", "average time of event from gen matching", 200, -0.5, 0.5);
  h_vertexmatchedavgtime = new TH1F("vertexmatchedavg_time", "average time of event from vertex matching", 200, -0.5, 0.5);
  h_avgtime = new TH2F("avg_time", "average time of event", 200, -0.5, 0.5, 200, -0.5, 0.5);
  h_avgtime->Sumw2();
  hprofile = new TProfile2D("profile", "time vs eta and phi", 12, -2.4, 2.4, 9, -3.14, 3.14);
  hprofile->Sumw2();
  h_genmatchedavgtime->Sumw2();
  h_vertexmatchedavgtime->Sumw2();
  h_nogenmatched_timeinfo = new TH1F("nogenmatched_timeinfo", "timeinfo", 2, -0.5, 1.5);
  h_eta_wtime = new TH1F("eta_wtime", "eta info with time", 50, 0., 3.);
  h_eta_wtime->Sumw2();
  h_vertex = new TH1F("vertex", "vertex size", 250, 0., 250.);
  h_vertex->Sumw2();
  h_tau_rec_vertex = new TH2F("tau_rec_vertex", "gen vs rec vertex", 50, 0., 10., 50, 0., 10.);
  h_tau_rec_vertex->Sumw2();
  h_pt_wtime = new TH1F("pt_wtime", "pt info with time", 100, 0, 50.);
  h_pt_wtime->Sumw2();
  h_eta_wotime = new TH1F("eta_wotime", "eta info w/o time", 50, 0., 3.);
  h_eta_wotime->Sumw2();
  h_pt_wotime = new TH1F("pt_wotime", "pt info w/o time", 100, 0, 50.);
  h_pt_wotime->Sumw2();
  h_nogenmatched_timeinfo->GetXaxis()->SetBinLabel(1, "without_timeinfo");
  h_nogenmatched_timeinfo->GetXaxis()->SetBinLabel(2, "with_timeinfo");
  h_genmatched_timeinfo->GetXaxis()->SetBinLabel(1, "without_timeinfo");
  h_genmatched_timeinfo->GetXaxis()->SetBinLabel(2, "with_timeinfo");
}
// ------------ method called once each job just after ending the event loop  ------------
void
HLTTauTimeStudy::endJob()
{
  h_genmatchedavgtime->Write();
  h_vertexmatchedavgtime->Write();
  h_avgtime->Write();
  h_vertex->Write();
  hprofile->Write();
  h_eta_wtime->Write();
  h_pt_wtime->Write();
  h_eta_wotime->Write();
  h_pt_wotime->Write();
  h_genmatched_timeinfo->Write();
  h_nogenmatched_timeinfo->Write();
  h_tau_rec_vertex->Write();
  for ( int particle_id=ParticleType::X; particle_id<=ParticleType::egamma_HF; particle_id++ )
  {
    gdir[particle_id].cd();
    for ( int match=GenMatch::kNoGenMatch; match<=GenMatch::kGenMatch; match++)
    {
      for ( std::map<std::string, TH1F*>::const_iterator itr = genMatched_histogram[match][particle_id].begin(); itr != genMatched_histogram[match][particle_id].end(); itr++)
      {
        itr->second->Write();
      }
    }
  }
}
  
//define this as a plug-in
DEFINE_FWK_MODULE(HLTTauTimeStudy);
