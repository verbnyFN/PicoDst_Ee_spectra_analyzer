// C++ headers
#include <iostream>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include <cmath>
#include "TLorentzVector.h"
#include <deque>
#include <vector>
#include <tuple>

// PicoDst headers
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"
#include "StPicoBTofHit.h"
#include "StPicoBTowHit.h"
#include "StPicoEmcTrigger.h"
#include "StPicoBTofPidTraits.h"
#include "StPicoTrackCovMatrix.h"
#include "StPicoFmsHit.h"
#include "StPicoETofHit.h"
#include "StPicoEpdHit.h"

const int nVzBins = 10;
const int nCentBins = 9;
const int poolDepth = 10;


int getVzBin(float vz)
{
    float vzMin = -40.0;
    float vzMax = 40.0;
    float binWidth = (vzMax - vzMin) / nVzBins;

    int bin = (vz - vzMin) / binWidth;

    if (bin < 0) bin = 0;
    if (bin >= nVzBins) bin = nVzBins - 1;

    return bin;
}

int getCentBin(int refMult, const std::vector<int>& centBins) {
    for (int i = 0; i < (int)centBins.size(); i++) {
        if (refMult >= centBins[i]) return i;
    }
    return centBins.size(); 
}

//_________________
int main(int argc, char* argv[]) {

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  R__LOAD_LIBRARY(libStPicoDst);
#else
  gSystem->Load("../libs/libStPicoDst.so");
#endif

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  const char* fileName;
  const char* oFileName;

  switch (argc) {
  case 3:
    fileName = argv[1];
    oFileName = argv[2];
    break;
  default:
    std::cout << "Usage: picoAnalyzerStandalone inputFileName outputFileName.root" << std::endl;
    return -1;
  }
  std::cout << " inputFileName : " << fileName << std::endl;
  std::cout << " outputFileName: " << oFileName << std::endl;
  
  StPicoDstReader* picoReader = new StPicoDstReader(fileName);
  picoReader->Init();

  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event*", 1);
  picoReader->SetStatus("Track*", 1);
  picoReader->SetStatus("BTofHit*", 1);
  picoReader->SetStatus("BTofPidTraits*", 1);
  picoReader->SetStatus("BTowHit*", 1);
  picoReader->SetStatus("ETofHit*", 1);
  picoReader->SetStatus("EpdHit*", 1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !picoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;


  TFile *oFile = new TFile(oFileName, "recreate");
  
  // Histogramming
  // Event
  TH1D *hRefMult = new TH1D("hRefMult",
			    "Reference multiplicity;refMult",
			    500, 0, 500);
  TH2D *hVtxXvsY = new TH2D("hVtxXvsY",
			    "hVtxXvsY",
			    200,-10.,10.,200,-10.,10.);
  TH1D *hVtxZ = new TH1D("hVtxZ","hVtxZ",
			 140, -70., 70.);

  // Track
  TH1D *hGlobalPtot = new TH1D("hGlobalPtot",
			       "Global track momentum;p (GeV/c)",
			       100, 0., 1. );
  TH1D *hGlobalPtotCut = new TH1D("hGlobalPtotCut",
				  "Global track momentum after cut;p (GeV/c)",
				  100, 0., 1. );
  TH1D *hPrimaryPtot = new TH1D("hPrimaryPtot",
				"Primary track momentum;p (GeV/c)",
			       100, 0., 1. );
  TH1D *hPrimaryPtotCut = new TH1D("hPrimaryPtotCut",
				   "Primary track momentum after cut;p (GeV/c)",
				  100, 0., 1. );
  TH1D *hTransvMomentum = new TH1D("hTransvMomentum",
				   "Track transverse momentum;p_{T} (GeV/c)",
				   200, 0., 2.);
  TH2D *hGlobalPhiVsPt[2];
  for(int i=0; i<2; i++) {
    hGlobalPhiVsPt[i] = new TH2D(Form("hGlobalPhiVsPt_%d",i),
				 Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
				 300, 0., 3.,
				 630, -3.15, 3.15);
  }
  TH1D *hNSigmaPion = new TH1D("hNSigmaPion",
			       "n#sigma(#pi);n#sigma(#pi)",
			       400, -10., 10.);
  TH1D *hNSigmaElectron = new TH1D("hNSigmaElectron",
				   "n#sigma(e);n#sigma(e)",
				   400,-10.,10.);
  TH1D *hNSigmaKaon = new TH1D("hNSigmaKaon",
			       "n#sigma(K);n#sigma(K)",
			       400, -10., 10.);
  TH1D *hNSigmaProton = new TH1D("hNSigmaProton",
				 "n#sigma(p);n#sigma(p)",
				 400, -10., 10.);
    
  // BTof pid traits
  TH1D *hTofBeta = new TH1D("hTofBeta", "BTofPidTraits #beta;#beta",
			    2000, 0., 2.);

  // BTOF hit
  TH1D *hBTofTrayHit = new TH1D("hBTofTrayHit","BTof tray number with the hit",
				120, -0.5, 119.5);

  // BTOW hit
  TH1D *hBTowAdc = new TH1D("hBTowAdc","Barrel tower ADC;ADC",500,0.,500);

  // FMS hit
  TH1D *hFmsAdc = new TH1D("hFmsAdc","ADC in FMS modules;ADC",1000, 0.,5000);

  // ETOF hit
  TH1D *hETofToT = new TH1D("hETofToT","eTOF TOT;Time over threshold (ns)",300, 0.,150);

  // EPD hit
  TH1D *hEpdAdc = new TH1D("hEpdAdc","ADC in EPD;ADC",4095, 0., 4095);

  // Vx vs Vy
  TH2D *hVxVy_before = new TH2D("hVxVy_before", "Vx vs Vy (before event cut);V_{x} (cm);V_{y} (cm)", 2000, -10, 10, 2000, -10, 10);
  TH2D *hVxVy_after  = new TH2D("hVxVy_after",  "Vx vs Vy (after event cut);V_{x} (cm);V_{y} (cm)", 500, -2.5, 2.5, 500, -2.5, 2.5);
  // Vz
  TH1D *hVz_before = new TH1D("hVz_before", "V_{z} (before event cut);V_{z} (cm)", 400, -100, 100);
  TH1D *hVz_after  = new TH1D("hVz_after",  "V_{z} (after event cut);V_{z} (cm)", 200, -50, 50);
  // RefMult
  TH1D *hRefMult_after = new TH1D("hRefMult_after", "Reference multiplicity (after event selection);refMult", 1500, 0, 1500);

  // nHitsFit
  TH1D *hNfit_before = new TH1D("hNfit_before", "nHitsFit (before track cut);nHitsFit", 60, 0, 60);
  TH1D *hNfit_after  = new TH1D("hNfit_after",  "nHitsFit (after track cut);nHitsFit", 60, 0, 60);

  // nHitsDedx 
  TH1D *hNdEdx_before = new TH1D("hNdEdx_before", "nHitsDedx (before track cut);nHitsDedx", 40, 0, 40);
  TH1D *hNdEdx_after  = new TH1D("hNdEdx_after",  "nHitsDedx (after track cut);nHitsDedx", 40, 0, 40);

  // DCA 
  TH1D *hDCA_before = new TH1D("hDCA_before", "|DCA| (before track cut);|DCA| (cm)", 400, 0, 10);
  TH1D *hDCA_after  = new TH1D("hDCA_after",  "|DCA| (after track cut);|DCA| (cm)", 200, 0, 5);

  // pT
  TH1D *hPt_before = new TH1D("hPt_before", "p_{T} (before track cut);p_{T} (GeV/c)", 200, 0, 5);
  TH1D *hPt_after  = new TH1D("hPt_after",  "p_{T} (after track cut);p_{T} (GeV/c)", 200, 0, 5);

  // eta
  TH1D *hEta_before = new TH1D("hEta_before", "#eta (before track cut);#eta", 200, -2, 2);
  TH1D *hEta_after  = new TH1D("hEta_after",  "#eta (after track cut);#eta", 200, -2, 2);

  // pT vs eta
  TH2D *hPtEta_before = new TH2D("hPtEta_before","p_{T} vs #eta (before track cut);#eta;p_{T} (GeV/c)", 400, -2, 2, 10000, 0, 50);
  TH2D *hPtEta_after  = new TH2D("hPtEta_after", "p_{T} vs #eta (after track cut);#eta;p_{T} (GeV/c)", 400, -2, 2, 10000, 0, 50);

  // TPC dE/dx vs p/q
  TH2D *hDedxVsPq_before = new TH2D("hDedxVsPq_before", "TPC dE/dx vs p/q (before PID cut);p/q (GeV/c);dE/dx (keV/cm)", 3000, 0, 30, 500, 0, 50);
  TH2D *hDedxVsPq_after  = new TH2D("hDedxVsPq_after", "TPC dE/dx vs p/q (after PID cut);p/q (GeV/c);dE/dx (keV/cm)", 300, 0, 10, 600, 0, 6);

  // nSigma(e) vs p/q 
  TH2D *hNsigmaE_vs_Pq_before = new TH2D("hNsigmaE_vs_Pq_before", "n#sigma(e) vs p/q (before PID);p/q;n#sigma(e)", 300, 0, 10, 750, -20, 50);
  TH2D *hNsigmaE_vs_Pq_after  = new TH2D("hNsigmaE_vs_Pq_after", "n#sigma(e) vs p/q (after PID);p/q;n#sigma(e)", 300, 0, 10, 250, -10, 15);

  // TOF 1/beta vs p/q 
  TH2D *hInvBeta_vs_Pq_before = new TH2D("hInvBeta_vs_Pq_before", "1/#beta vs p/q (before PID);p/q;1/#beta", 1000, 0, 10, 1000, 0, 10);
  TH2D *hInvBeta_vs_Pq_after  = new TH2D("hInvBeta_vs_Pq_after", "1/#beta vs p/q (after PID);p/q;1/#beta", 1000, 0, 10, 100, 0.5, 1.5);

  // BEMC E/p vs p/q
  TH2D *hEoverP_vs_Pq_before = new TH2D("hEoverP_vs_Pq_before", "P/E vs p/q (before PID);p/q;P/E", 300, 0, 4, 300, 0, 3);
  TH2D *hEoverP_vs_Pq_after  = new TH2D("hEoverP_vs_Pq_after", "P/E vs p/q (after PID);p/q;P/E", 300, 0, 4, 300, 0, 3);

  // inv mass
  TH1D* hInvMassEE = new TH1D("hInvMassEE", "Invariant mass of e+ e- pairs from same events;M_{inv} (GeV/c^{2});Counts", 1000, 0, 10.0);
  TH1D* hInvMassEEMixed = new TH1D("hInvMassEEMixed", "Invariant mass of e+ e- pairs from mixed events;M_{inv} (GeV/c^{2});Counts", 1000, 0, 10.0);
  TH1D* hInvMassEEFinal = new TH1D("hInvMassEEFinal", "Invariant mass of e+ e- pairs without background;M_{inv} (GeV/c^{2});Counts", 1000, 0, 10.0);
  
  // Пул только для электронов — миксуем только в одну сторону
  std::deque< std::vector<TLorentzVector> > poolElectrons[nVzBins][nCentBins];
  
  float mass_e = 0.000511;
  std::vector<TLorentzVector> electrons;
  std::vector<TLorentzVector> positrons;

  // Вектор кандидатов для conversion cut: {заряд, индекс трека, 4-вектор}
  std::vector<std::tuple<int, int, TLorentzVector>> eCandidates;

  // =====================================================================
  // Первый проход — определение centrality
  // =====================================================================
  for (Long64_t iEvent = 0; iEvent < events2read; iEvent++) {
    picoReader->readPicoEvent(iEvent);
    StPicoEvent *event = picoReader->picoDst()->event();
    if (!event) continue;    

    hRefMult->Fill( event->refMult() );
    if (abs(event->primaryVertex().Z()) >= 40) continue;
    if (sqrt(pow(event->primaryVertex().X(),2) + pow(event->primaryVertex().Y(),2)) >= 2) continue;
    hRefMult_after->Fill(event->refMult());
    std::cout << "First pass, event #[" << (iEvent+1) << "/" << events2read << "]" << std::endl;
  }
  
  double centEdges[] = {0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80};
  std::vector<int> centBins;
  double totalEvents = hRefMult_after->Integral(1, hRefMult_after->GetNbinsX());

  for (int i = 0; i < nCentBins; i++) {
      double target = centEdges[i] * totalEvents;
      double multSum = 0;
      for (int bin = hRefMult_after->GetNbinsX(); bin >= 1; bin--) {
          multSum += hRefMult_after->GetBinContent(bin);
          if (multSum >= target) {
              centBins.push_back(bin);
              break;
          }
      }
  }

  std::cout << "Centrality bins: ";
  for (int i = 0; i < (int)centBins.size(); i++) { std::cout << centBins[i] << " "; }
  std::cout << std::endl;
  
  picoReader->Finish();

  // =====================================================================
  // Второй проход — основной анализ
  // =====================================================================
  StPicoDstReader* myReaderSecond = new StPicoDstReader(fileName);
  myReaderSecond->Init();
  std::cout << "Explicit read status for some branches" << std::endl;
  myReaderSecond->SetStatus("*",0);
  myReaderSecond->SetStatus("Event*", 1);
  myReaderSecond->SetStatus("Track*", 1);
  myReaderSecond->SetStatus("BTofHit*", 1);
  myReaderSecond->SetStatus("BTofPidTraits*", 1);
  myReaderSecond->SetStatus("BTowHit*", 1);
  myReaderSecond->SetStatus("ETofHit*", 1);
  myReaderSecond->SetStatus("EpdHit*", 1);

  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    std::cout << "Second pass, event #[" << (iEvent+1) << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = myReaderSecond->readPicoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      break;
    }

    StPicoDst *dst = myReaderSecond->picoDst();
    StPicoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    hVxVy_before->Fill(event->primaryVertex().X(), event->primaryVertex().Y());
    hVz_before->Fill(event->primaryVertex().Z());

    if (abs(event->primaryVertex().Z()) >= 40) {continue;}
    if (sqrt(event->primaryVertex().X()*event->primaryVertex().X() + event->primaryVertex().Y()*event->primaryVertex().Y()) >= 2) {continue;}

    hVxVy_after->Fill(event->primaryVertex().X(), event->primaryVertex().Y());
    hVz_after->Fill(event->primaryVertex().Z());

    TVector3 pVtx = event->primaryVertex();
    hVtxXvsY->Fill( event->primaryVertex().X(), event->primaryVertex().Y() );
    hVtxZ->Fill( event->primaryVertex().Z() );

    Int_t nTracks = dst->numberOfTracks();

    // =================================================================
    // Трековая петля — собираем кандидатов в eCandidates
    // =================================================================
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

      StPicoTrack *picoTrack = dst->track(iTrk);
      if(!picoTrack) continue;

      hGlobalPtot->Fill( picoTrack->gMom().Mag() );
      if( picoTrack->isPrimary() ) {
        hPrimaryPtot->Fill( picoTrack->pMom().Mag() );
      }
      
      if( picoTrack->gMom().Mag() < 0.1 ||
          picoTrack->gDCA(pVtx).Mag()>50. ) {
        continue;
      } 

      hNfit_before->Fill(picoTrack->nHitsFit());
      hNdEdx_before->Fill(picoTrack->nHitsDedx());
      hDCA_before->Fill(picoTrack->gDCA(pVtx).Mag());
      hPt_before->Fill(picoTrack->pPt());
      hEta_before->Fill(picoTrack->pMom().PseudoRapidity());
      hPtEta_before->Fill(picoTrack->pMom().PseudoRapidity(), picoTrack->pPt());

      if (!picoTrack->isPrimary()) {continue;}
      if ( TMath::Abs(picoTrack->pMom().PseudoRapidity()) >= 1.0 ) {continue;}
      if (picoTrack->pPt() >= 50 || picoTrack->pPt() <= 0.2) {continue;}
      if (picoTrack->gDCA(event->primaryVertex()).Mag() > 1.0) {continue;}
      if ((float)picoTrack->nHitsFit() / picoTrack->nHitsMax() <= 0.52 ) {continue;}
      if (picoTrack->nHitsFit() < 20) {continue;}
      if (picoTrack->nHitsDedx() < 11) {continue;}

      hNfit_after->Fill(picoTrack->nHitsFit());
      hNdEdx_after->Fill(picoTrack->nHitsDedx());
      hDCA_after->Fill(picoTrack->gDCA(pVtx).Mag());
      hPt_after->Fill(picoTrack->pPt());
      hEta_after->Fill(picoTrack->pMom().PseudoRapidity());
      hPtEta_after->Fill(picoTrack->pMom().PseudoRapidity(), picoTrack->pPt());

      float p = picoTrack->pMom().Mag();
      float pq = p / picoTrack->charge();
      float dedx = picoTrack->dEdx();

      hDedxVsPq_before->Fill(pq, dedx);
      hNsigmaE_vs_Pq_before->Fill(pq, picoTrack->nSigmaElectron());
          
      if (picoTrack->isTofTrack()) {
          StPicoBTofPidTraits* tof = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());
          if (tof && tof->btofBeta() > 0) {
              hInvBeta_vs_Pq_before->Fill(pq, 1.0 / tof->btofBeta());
          }
      }
      
      if (picoTrack->isBemcTrack()) {
          StPicoBTowHit* towHit = dst->btowHit(picoTrack->bemcTowerIndex());
          if (towHit && towHit->energy() > 0) {
              hEoverP_vs_Pq_before->Fill(pq, p / towHit->energy());
          }
      }

      // PID cuts — идентификация электронов
      if (picoTrack->nSigmaElectron() < -1.9 || picoTrack->nSigmaElectron() > 3) {continue;}

      if (picoTrack->isBemcTrack()) {
        StPicoBTowHit* towHit = dst->btowHit(picoTrack->bemcTowerIndex());
        if (towHit) { 
          // if(picoTrack->pPt() <= 3.5 || 0.3 >= (picoTrack->gMom().Mag() / towHit->energy()) || 1.5 <= (picoTrack->gMom().Mag() / towHit->energy())) {continue;}
        }
      }

      if (picoTrack->isTofTrack()) {
        StPicoBTofPidTraits* tof = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());
        if (tof) {
          if (fabs(1.0/tof->btofBeta() - 1.0) >= 0.03 || fabs(tof->btofYLocal()) >= 10) {continue;} 
        } else {continue;}
      } else {continue;}

      hDedxVsPq_after->Fill(pq, dedx);
      hNsigmaE_vs_Pq_after->Fill(pq, picoTrack->nSigmaElectron());
          
      if (picoTrack->isTofTrack()) {
          StPicoBTofPidTraits* tof = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());
          if (tof && tof->btofBeta() > 0) {
              hInvBeta_vs_Pq_after->Fill(pq, 1.0 / tof->btofBeta());
          }
      }
      
      if (picoTrack->isBemcTrack()) {
          StPicoBTowHit* towHit = dst->btowHit(picoTrack->bemcTowerIndex());
          if (towHit && towHit->energy() > 0) {
              hEoverP_vs_Pq_after->Fill(pq, p / towHit->energy());
          }
      }

      TLorentzVector lv;
      lv.SetPtEtaPhiM(picoTrack->pPt(), picoTrack->pMom().PseudoRapidity(), picoTrack->pMom().Phi(), mass_e);

      // Сохраняем кандидата: {заряд, индекс, 4-вектор}
      eCandidates.push_back({picoTrack->charge(), iTrk, lv});

      hGlobalPtotCut->Fill( picoTrack->gMom().Mag() );
      if( picoTrack->isPrimary() ) {
        hPrimaryPtotCut->Fill( picoTrack->pMom().Mag() );
      }
      if( picoTrack->charge() > 0 ) {
        hGlobalPhiVsPt[0]->Fill( picoTrack->gMom().Pt(), picoTrack->gMom().Phi() );
      } else {
        hGlobalPhiVsPt[1]->Fill( picoTrack->gMom().Pt(), picoTrack->gMom().Phi() );	
      }
      hNSigmaElectron->Fill( picoTrack->nSigmaElectron() );
      hNSigmaPion->Fill( picoTrack->nSigmaPion() );
      hNSigmaKaon->Fill( picoTrack->nSigmaKaon() );
      hNSigmaProton->Fill( picoTrack->nSigmaProton() );
      hTransvMomentum->Fill( picoTrack->gMom().Pt() );

      if( picoTrack->isTofTrack() ) {
        StPicoBTofPidTraits *trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
        if( !trait ) {
          std::cout << "O-oh... No BTofPidTrait # " << picoTrack->bTofPidTraitsIndex()
                    << " for track # " << iTrk << std::endl;
          continue;
        }
        hTofBeta->Fill( trait->btofBeta() );
      }
      
    } // for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    // =================================================================
    // Conversion cut — проверяем кандидатов попарно
    // =================================================================
    for (auto& [charge1, idx1, lv1] : eCandidates) {
        bool isConversion = false;
        for (auto& [charge2, idx2, lv2] : eCandidates) {
            if (idx1 == idx2) continue;
            if (charge1 * charge2 > 0) continue; // только противоположные заряды
            if ((lv1 + lv2).M() < 0.05) {
                isConversion = true;
                break;
            }
        }
        if (isConversion) continue;

        if (charge1 > 0) positrons.push_back(lv1);
        else             electrons.push_back(lv1);
    }
    eCandidates.clear();

    // =================================================================
    // Hit analysis
    // =================================================================
    Int_t nBTofHits = dst->numberOfBTofHits();
    for(Int_t iHit=0; iHit<nBTofHits; iHit++) {
      StPicoBTofHit *btofHit = dst->btofHit(iHit);
      if( !btofHit ) continue;
      hBTofTrayHit->Fill( btofHit->tray() );
    }

    Int_t nBTowHits = dst->numberOfBTowHits();
    for(Int_t iHit=0; iHit<nBTowHits; iHit++) {
      StPicoBTowHit *btowHit = dst->btowHit(iHit);
      if( !btowHit ) continue;
      hBTowAdc->Fill( btowHit->adc() );
    }

    Int_t nFmsHits = dst->numberOfFmsHits();
    for(Int_t iHit=0; iHit<nFmsHits; iHit++) {
      StPicoFmsHit *fmsHit = dst->fmsHit(iHit);
      if( !fmsHit ) continue;
      hFmsAdc->Fill( fmsHit->adc() );
    }

    Int_t nETofHits = dst->numberOfETofHits();
    for(Int_t iHit=0; iHit<nETofHits; iHit++) {
      StPicoETofHit *etofHit = dst->etofHit(iHit);
      if( !etofHit ) continue;
      hETofToT->Fill( etofHit->timeOverThreshold() );
    }
    
    Int_t nEpdHits = dst->numberOfEpdHits();
    for(Int_t iHit=0; iHit<nEpdHits; iHit++) {
      StPicoEpdHit *epdHit = dst->epdHit(iHit);
      if( !epdHit ) continue;
      hEpdAdc->Fill( epdHit->adc() );
    }

    // =================================================================
    // Инвариантная масса — same events
    // =================================================================
    for (auto& eplus : positrons) {
      for (auto& eminus : electrons) {
        hInvMassEE->Fill((eplus + eminus).M());
      }
    }
    
    // =================================================================
    // Centrality и Vz бины
    // =================================================================
    int vzBin  = getVzBin(event->primaryVertex().Z());
    int centBin = getCentBin(event->refMult(), centBins);
    if (centBin >= nCentBins) {
      positrons.clear();
      electrons.clear();
      continue;
    }

    // =================================================================
    // Mixed event background — только в одну сторону:
    // новые позитроны × старые электроны
    // =================================================================
    for (auto &oldElectrons : poolElectrons[vzBin][centBin]) {
      for (auto &eplus : positrons) {
          for (auto &eminus : oldElectrons) {
              hInvMassEEMixed->Fill((eplus + eminus).M());
          }
      }
    }

    // Обновляем пул электронами текущего события
    poolElectrons[vzBin][centBin].push_back(electrons);
    if (poolElectrons[vzBin][centBin].size() > poolDepth) {
      poolElectrons[vzBin][centBin].pop_front();
    }

    positrons.clear();
    electrons.clear();

  } // for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  // =================================================================
  // Нормировка и вычитание фона
  // =================================================================
  double mLow  = 1.2, mHigh = 2.5;
  int bLow  = hInvMassEE->FindBin(mLow);
  int bHigh = hInvMassEE->FindBin(mHigh);

  double sameSB  = hInvMassEE->Integral(bLow, bHigh);
  double mixedSB = hInvMassEEMixed->Integral(bLow, bHigh);
  double norm_factor = sameSB / mixedSB;

  std::cout << "Normalization factor: " << norm_factor << std::endl;

  hInvMassEEMixed->Scale(norm_factor);
  hInvMassEEFinal->Add(hInvMassEE, 1);
  hInvMassEEFinal->Add(hInvMassEEMixed, -1); 

  myReaderSecond->Finish();

  hVxVy_after->SetOption("colz");
  hVxVy_before->SetOption("colz");
  hPtEta_before->SetOption("colz"); 
  hPtEta_after->SetOption("colz");  
  hDedxVsPq_before->SetOption("colz");
  hDedxVsPq_after->SetOption("colz");
  hNsigmaE_vs_Pq_before->SetOption("colz"); 
  hNsigmaE_vs_Pq_after->SetOption("colz");
  hInvBeta_vs_Pq_before->SetOption("colz");
  hInvBeta_vs_Pq_after->SetOption("colz");

  oFile->Write();
  oFile->Close();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!" << std::endl;

  return 0;
}
