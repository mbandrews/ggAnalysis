#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "DQM/HcalCommon/interface/Constants.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
//#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

// Fill ECAL rechits at HCAL granularity /////////////////
// Project ECAL event rechits into a vector of length
// equal to number of *towers* in HBHE (iphi:72, ieta:56).
//
// NOTE: We do not decrease the iphi granularity as
// happens in the real HE. The findDetIdCalo() enforces
// the even iphi numbering in the coarse region so must
// resort to intermediate helper histograms. Since helper
// histograms are binned by eta,phi some approx. involved.

TH2F *hEvt_HBHE_EMenergy;
TProfile2D *hHBHE_EMenergy;
std::vector<float> vHBHE_EMenergy_;

// HBHE eta bin edges
double eta_bins_HBHE[2*(hcaldqm::constants::IETA_MAX_HE-1)+1] =
                  {-3.000, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
                   -1.218, -1.131, -1.044, -0.957, -0.870, -0.783, -0.695, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000,
                    0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,  0.695,  0.783,  0.870,  0.957,  1.044,  1.131,  1.218,
                    1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,  1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  3.000}; // 57

// Initialize branches _____________________________________________________________//
void ggNtuplizer::branchesECALatHCAL ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("HBHE_EMenergy",    &vHBHE_EMenergy_);
  // Intermediate helper histogram (single event only)
  hEvt_HBHE_EMenergy = new TH2F("evt_HBHE_EMenergy", "E(#phi,#eta);#phi;#eta",
      hcaldqm::constants::IPHI_NUM,         -TMath::Pi(),     TMath::Pi(),
      2*(hcaldqm::constants::IETA_MAX_HE-1), eta_bins_HBHE );

  // Histograms for monitoring
  hHBHE_EMenergy = fs->make<TProfile2D>("HBHE_EMenergy", "E(i#phi,i#eta);i#phi;i#eta",
      hcaldqm::constants::IPHI_NUM,           hcaldqm::constants::IPHI_MIN-1,    hcaldqm::constants::IPHI_MAX,
      2*(hcaldqm::constants::IETA_MAX_HE-1),-(hcaldqm::constants::IETA_MAX_HE-1),hcaldqm::constants::IETA_MAX_HE-1 );

} // branchesECALatHCAL

// Fill ECAL rechits at HBHE granularity ___________________________________________________//
void ggNtuplizer::fillECALatHCAL ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int ieta_, iphi_, idx_;
  float eta,  phi, energy_;
  const CaloGeometry* caloGeom;
  GlobalPoint pos;

  vHBHE_EMenergy_.assign( 2*hcaldqm::constants::IPHI_NUM*(hcaldqm::constants::IETA_MAX_HE-1),0. );
  hEvt_HBHE_EMenergy->Reset();

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByToken( EBRecHitCollectionT_, EBRecHitsH_ );
  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByToken( EERecHitCollectionT_, EERecHitsH_ );
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  caloGeom = caloGeomH_.product();

  // Fill EB rechits
  for ( EcalRecHitCollection::const_iterator iRHit = EBRecHitsH_->begin();
        iRHit != EBRecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ == 0. ) continue;
    // Get position of cell centers
    pos = caloGeom->getPosition( iRHit->id() );
    eta = pos.eta();
    phi = pos.phi();
    // Fill intermediate helper histogram by eta,phi
    hEvt_HBHE_EMenergy->Fill( phi, eta, energy_ );

    //HcalDetId hId( spr::findDetIdHCAL( caloGeom, eta, phi, false ) );

  } // EB rechits

  // Fill EE rechits
  for ( EcalRecHitCollection::const_iterator iRHit = EERecHitsH_->begin();
        iRHit != EERecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ == 0. ) continue;
    // Get position of cell centers
    pos = caloGeom->getPosition( iRHit->id() );
    eta = pos.eta();
    phi = pos.phi();
    // Fill intermediate helper histogram by eta,phi
    hEvt_HBHE_EMenergy->Fill( phi, eta, energy_ );

    //HcalDetId hId( spr::findDetIdHCAL( caloGeom, eta, phi, false ) );

  } // EE rechits

  // Fill vector for full ECAL@HCAL image using helper histograms
  for (int ieta = 1; ieta < hEvt_HBHE_EMenergy->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    for (int iphi = 1; iphi < hEvt_HBHE_EMenergy->GetNbinsX()+1; iphi++) {

      energy_ = hEvt_HBHE_EMenergy->GetBinContent( iphi, ieta );
      if ( energy_ == 0. ) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 38; // shift
      iphi_ = iphi_ > hcaldqm::constants::IPHI_MAX ? iphi_-hcaldqm::constants::IPHI_MAX : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_*hcaldqm::constants::IPHI_NUM + iphi_;
      // Fill vector for image
      vHBHE_EMenergy_[idx_] = energy_;
      // Fill histogram for monitoring
      hHBHE_EMenergy->Fill( iphi_, ieta_-(hcaldqm::constants::IETA_MAX_HE-1), energy_ );

    } // iphi
  } // ieta

} // fillECALatHCAL()
