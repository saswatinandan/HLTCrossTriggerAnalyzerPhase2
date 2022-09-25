// Class:      HLTTauAnalyzer
//
// Original Author:  Sandeep Bhowmik
//         Created:  Tue, 12 Mar 2019 18:38:39 GMT
//
#include "FWCore/Framework/interface/one/EDAnalyzer.h" 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TH2F.h>
#include <TProfile2D.h>

#include "DataFormats/VertexReco/interface/Vertex.h"  
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <DataFormats/PatCandidates/interface/Tau.h> 

#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"   // reco::PFTauDiscriminator

#include "DataFormats/L1TParticleFlow/interface/PFTau.h"       // l1t::PFTauCollection for NNTau
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Math/interface/Point3D.h"
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
  std::map<int, std::map<int, std::map<std::string, TH1F*> >> histograms;
  TProfile2D* hprofile;
  TH1F* h_analyzed;
  TH1F* h_vertex;
  TH1F* h_genmatched_timeinfo;
  TH1F* h_nogenmatched_timeinfo;
  TH1F* h_eta_wtime;
  TH1F* h_pt_wtime;
  TH1F* h_eta_wotime;
  TH1F* h_pt_wotime;
  TH1F* h_genmatchedavgtime;
  TH1F* h_vertexmatchedavgtime;
  TH2F* h_avgtime;
  std::map<int, TFileDirectory> gdir;
  ULong64_t       indexevents_;
  Int_t           runNumber_;
  Int_t           lumi_;
  std::map<TDirectory *, std::vector<TH1 *>> gHistograms_;

  // ----------member data ---------------------------
  bool debug_;
  bool usegenJets_;
  edm::EDGetTokenT<std::vector<reco::Vertex>>      recoVertexToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection>      pfCandToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<vector<pat::Tau>>      recoPatTauToken_;
  edm::EDGetTokenT<std::vector<reco::PFJet>> recoJetToken_;
  edm::EDGetTokenT<math::XYZPointF>  genvertexToken_;
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
  , genvertexToken_ (consumes<math::XYZPointF>    (iConfig.getParameter<edm::InputTag>("genvertexToken")))
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
      for ( std::map<std::string, TH1F*>::const_iterator itr = histograms[match][particle_id].begin(); itr != histograms[match][particle_id].end(); itr++)
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
  
  const auto& genVtxPosition = iEvent.get(genvertexToken_);
  math::XYZPoint vertex(genVtxPosition.x(), genVtxPosition.y(), genVtxPosition.z());

  int genmatched_time_count(0), vertexmatched_time_count(0);
  double genmatched_sum_time(0), vertexmatched_sum_time(0);

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
          fillWOverflow(histograms[match][particle_id]["time_distribution"], time);
          fillWOverflow(histograms[match][particle_id]["pt_distribution"], recpt);
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
            if ( fabs(diff.z()) < 0.1)
            {
              vertexmatched_time_count += recpt;
              vertexmatched_sum_time += (recpt*time);
            }
          }
        }
        if ( recpt > 2 && abs(receta) < 2.4 )fillWOverflow(histograms[match][particle_id]["time_info"], hastimeinfo);
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
            fillWOverflow(histograms[match][particle_id]["time_shift"], recoPFCand.time() - genmatched_avgtime );
            if ( genmatched_time_count >=3 && genmatched ) fillWOverflow(histograms[match][particle_id]["time_shift_with5"], recoPFCand.time() - genmatched_avgtime );
          }
          if ( vertexmatched_avgtime > -100 )
          {
            math::XYZPoint diff(recoPFCand.vertex().x() - vertex.x(), recoPFCand.vertex().y() - vertex.y(), recoPFCand.vertex().z() - vertex.z());
            bool dxy = TMath::Sqrt(square(diff.x()) + square(diff.y())) < 0.05;
            bool dz = fabs(diff.z()) < 0.1;
            fillWOverflow(histograms[dz][particle_id]["time_shift_vertexmatched"], recoPFCand.time() - vertexmatched_avgtime );
          }
        }
      }
    }
  }

  for (auto tau : *recoPatTauHandle) {
    int cand(0);
    float leadtime(0);
    bool leadgenmatched(0);
    reco::CandidatePtrVector signal = tau.signalCands();
    for (auto sig : signal)
    {
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
      if(cand ==1)
      {
        leadtime = packedCand->time();
        leadgenmatched = genmatch;
        fillWOverflow(histograms[genmatch][ParticleType::h]["tauleadcandidate"], 1);
        fillWOverflow(histograms[genmatch][ParticleType::h]["time_shift_leadcandidate"], leadtime-genmatched_avgtime);
        fillWOverflow(histograms[dz][ParticleType::h]["time_shift_leadcandidate_vertexmatched"], leadtime-vertexmatched_avgtime);
        continue;
      }
      else if (cand==2) {
        fillWOverflow(histograms[genmatch && leadgenmatched][ParticleType::h]["time_shift_of_1stsubleadCand_wrt_lead_chargedCand"], leadtime - packedCand->time());
        fillWOverflow(histograms[genmatch][ParticleType::h]["time_shift_subleadcandidate"], packedCand->time()-genmatched_avgtime);
        fillWOverflow(histograms[dz][ParticleType::h]["time_shift_subleadcandidate_vertexmatched"], packedCand->time()-vertexmatched_avgtime);
        continue;
      }
      else if (cand==3) {
        fillWOverflow(histograms[genmatch && leadgenmatched][ParticleType::h]["time_shift_of_2ndsubleadCand_wrt_lead_chargedCand"], leadtime - packedCand->time());
        fillWOverflow(histograms[genmatch][ParticleType::h]["time_shift_2ndsubleadcandidate"], packedCand->time()-genmatched_avgtime);
        fillWOverflow(histograms[dz][ParticleType::h]["time_shift_2ndsubleadcandidate_vertexmatched"], packedCand->time()-vertexmatched_avgtime);
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
  for ( int match=GenMatch::kNoGenMatch; match<=GenMatch::kGenMatch; match++)
  {
    for ( int particle_id=ParticleType::X; particle_id<=ParticleType::egamma_HF; particle_id++ )
    {
      std::string prefix = Form("%s",
         (match==GenMatch::kGenMatch) ? "genmatched" : "nogenmatched");
      TFileDirectory dir = fs->mkdir(map_particlename[particle_id]);
      gdir[particle_id] = dir;
      std::string name = Form("%s_time_distribution", prefix.data());
      histograms[match][particle_id]["time_distribution"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
      name = Form("%s_pt_distribution", prefix.data());
      histograms[match][particle_id]["pt_distribution"] = dir.make<TH1F>(name.data(), name.data(), 100, 0, 20);
      name = Form("%s_time_shift", prefix.data());
      histograms[match][particle_id]["time_shift"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
      if ( match==GenMatch::kGenMatch )
      {
        histograms[match][particle_id]["time_shift_with5"] = dir.make<TH1F>((name+"_with5").data(), (name+"_with5").data(), 200, -0.5, 0.5);
      }
      if ( particle_id == ParticleType::h )
      {
        name = Form("%s_time_shift_vertexmatched", prefix.data());
        histograms[match][particle_id]["time_shift_vertexmatched"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_leadcandidate", prefix.data());
        histograms[match][particle_id]["tauleadcandidate"] = dir.make<TH1F>(name.data(), name.data(), 1, 0.5, 1.5);
        name = Form("%s_time_shift_of_1stsubleadCand_wrt_lead_chargedCand", prefix.data());
        histograms[match][particle_id]["time_shift_of_1stsubleadCand_wrt_lead_chargedCand"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_of_2ndsubleadCand_wrt_lead_chargedCand", prefix.data());
        histograms[match][particle_id]["time_shift_of_2ndsubleadCand_wrt_lead_chargedCand"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_leadcandidate", prefix.data());
        histograms[match][ParticleType::h]["time_shift_leadcandidate"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_leadcandidate_vertexmatched", prefix.data());
        histograms[match][ParticleType::h]["time_shift_leadcandidate_vertexmatched"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_1stleadcandidate", prefix.data());
        histograms[match][ParticleType::h]["time_shift_subleadcandidate"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_1stleadcandidate_vertexmatched", prefix.data());
        histograms[match][ParticleType::h]["time_shift_subleadcandidate_vertexmatched"]= dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_2ndleadcandidate", prefix.data());
        histograms[match][ParticleType::h]["time_shift_2ndsubleadcandidate"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
        name = Form("%s_time_shift_2ndleadcandidate_vertexmatched", prefix.data());
        histograms[match][ParticleType::h]["time_shift_2ndsubleadcandidate_vertexmatched"] = dir.make<TH1F>(name.data(), name.data(), 200, -0.5, 0.5);
      }
      name = Form("%s_time_info", prefix.data());
      histograms[match][particle_id]["time_info"] = dir.make<TH1F>(name.data(), name.data(), 2, -0.5, 1.5);
    }
  }
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
  for ( int particle_id=ParticleType::X; particle_id<=ParticleType::egamma_HF; particle_id++ )
  {
    gdir[particle_id].cd();
    for ( int match=GenMatch::kNoGenMatch; match<=GenMatch::kGenMatch; match++)
    {
      for ( std::map<std::string, TH1F*>::const_iterator itr = histograms[match][particle_id].begin(); itr != histograms[match][particle_id].end(); itr++)
      {
        itr->second->Write();
      }
    }
  }
}
  
//define this as a plug-in
DEFINE_FWK_MODULE(HLTTauTimeStudy);
