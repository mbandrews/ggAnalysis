#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DQM/HcalCommon/interface/Constants.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

// Fill stitched EEm_EB_EEp image /////////////////////............/
// Store all ECAL event rechits into a stitched EEm_EB_EEp image 
// segmented in iphi,ieta spannning the full -3 < eta < 3. 
// Use EB-like granularity giving an extended range ieta=[-140,140].
//
// For endcaps, project EE hits into a helper histogram binned by 
// phi,eta before filling the full extended ECAL(iphi,eta) image.
// For barrel, fill EB hits directly since geometries are 1:1. 
//
// 'ieta_global' keeps track of the global ieta index count used
// for filling the extended image vector vECAL_energy.
// 'ieta_signed' keeps track of the position along [-140,140] used
// for filling the monitoring histogram hECAL_energy.

static const int ECAL_IETA_MAX_EXT = 140;
static const int HBHE_IETA_MAX_EB = hcaldqm::constants::IETA_MAX_HB + 1; // 17
static const int nEE = 2;

TH2F *hEvt_EE_energy[nEE];
TProfile2D *hECAL_energy;
std::vector<float> vECAL_energy_;

// EE-(phi,eta) projection eta edges
// These are generated by requiring 5 fictional crystals 
// to uniformly span each HCAL tower in eta (as in EB).
double eta_bins_EEm[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                  {-3.    , -2.93  , -2.86  , -2.79  , -2.72  , -2.65  , -2.62  ,
                   -2.59  , -2.56  , -2.53  , -2.5   , -2.4644, -2.4288, -2.3932,
                   -2.3576, -2.322 , -2.292 , -2.262 , -2.232 , -2.202 , -2.172 ,
                   -2.1462, -2.1204, -2.0946, -2.0688, -2.043 , -2.0204, -1.9978,
                   -1.9752, -1.9526, -1.93  , -1.91  , -1.89  , -1.87  , -1.85  ,
                   -1.83  , -1.812 , -1.794 , -1.776 , -1.758 , -1.74  , -1.7226,
                   -1.7052, -1.6878, -1.6704, -1.653 , -1.6356, -1.6182, -1.6008,
                   -1.5834, -1.566 , -1.5486, -1.5312, -1.5138, -1.4964, -1.479 }; // 56
// EE+(phi,eta) projection eta edges
double eta_bins_EEp[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                   {1.479 ,  1.4964,  1.5138,  1.5312,  1.5486,  1.566 ,  1.5834,
                    1.6008,  1.6182,  1.6356,  1.653 ,  1.6704,  1.6878,  1.7052,
                    1.7226,  1.74  ,  1.758 ,  1.776 ,  1.794 ,  1.812 ,  1.83  ,
                    1.85  ,  1.87  ,  1.89  ,  1.91  ,  1.93  ,  1.9526,  1.9752,
                    1.9978,  2.0204,  2.043 ,  2.0688,  2.0946,  2.1204,  2.1462,
                    2.172 ,  2.202 ,  2.232 ,  2.262 ,  2.292 ,  2.322 ,  2.3576,
                    2.3932,  2.4288,  2.4644,  2.5   ,  2.53  ,  2.56  ,  2.59  ,
                    2.62  ,  2.65  ,  2.72  ,  2.79  ,  2.86  ,  2.93  ,  3.    }; // 56

// Initialize branches _______________________________________________________________//
void ggNtuplizer::branchesECALstitched ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("ECAL_energy",    &vECAL_energy_);
  // Intermediate helper histogram (single event only)
  hEvt_EE_energy[0] = new TH2F("evt_EEm_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI, -TMath::Pi(), TMath::Pi(),
      5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_energy[1] = new TH2F("evt_EEp_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI, -TMath::Pi(), TMath::Pi(),
      5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  // Histograms for monitoring
  hECAL_energy = fs->make<TProfile2D>("ECAL_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EBDetId::MAX_IPHI,    EBDetId::MIN_IPHI-1, EBDetId::MAX_IPHI,
      2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );

} // branchesECALstitched()

// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
void fillECAL_with_EEproj ( TH2F *hEvt_EE_energy_, int ieta_global_offset, int ieta_signed_offset ) {

  int ieta_global_, ieta_signed_;
  int ieta_, iphi_, idx_;
  float energy_;

  for (int ieta = 1; ieta < hEvt_EE_energy_->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    ieta_global_ = ieta_ + ieta_global_offset;
    ieta_signed_ = ieta_ + ieta_signed_offset;
    for (int iphi = 1; iphi < hEvt_EE_energy_->GetNbinsX()+1; iphi++) {

      energy_ = hEvt_EE_energy_->GetBinContent( iphi, ieta );
      if ( energy_ == 0. ) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 5*38; // shift
      iphi_ = iphi_ > EBDetId::MAX_IPHI ? iphi_-EBDetId::MAX_IPHI : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_global_*EBDetId::MAX_IPHI + iphi_;
      // Fill vector for image
      vECAL_energy_[idx_] = energy_;
      // Fill histogram for monitoring
      hECAL_energy->Fill( iphi_, ieta_signed_, energy_ );

    } // iphi_
  } // ieta_

} // fillECAL_with_EEproj

// Fill stitched EE-, EB, EE+ rechits ________________________________________________________//
void ggNtuplizer::fillECALstitched ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int iphi_, ieta_, iz_, idx_;
  int ieta_global, ieta_signed;
  int ieta_global_offset, ieta_signed_offset;
  float eta, phi, energy_;
  const CaloGeometry* caloGeom;
  GlobalPoint pos;

  vECAL_energy_.assign( 2*ECAL_IETA_MAX_EXT*EBDetId::MAX_IPHI, 0. );
  for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_energy[iz]->Reset();

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByToken( EBRecHitCollectionT_, EBRecHitsH_ );
  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByToken( EERecHitCollectionT_, EERecHitsH_ );
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  caloGeom = caloGeomH_.product();

  // Fill the EE+/-(phi,eta) projection with the EE hits.
  for ( EcalRecHitCollection::const_iterator iRHit = EERecHitsH_->begin();
        iRHit != EERecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ == 0. ) continue;
    // Get detector id
    EEDetId eeId( iRHit->id() );
    iz_ = ( eeId.zside() > 0 ) ? 1 : 0;
    // Get position of cell centers
    pos  = caloGeom->getPosition( eeId );
    eta = pos.eta();
    phi = pos.phi();
    // Fill intermediate helper histogram by eta,phi
    hEvt_EE_energy[iz_]->Fill( phi, eta, iRHit->energy() );

  } // EE+/-

  // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
  ieta_global_offset = 0;
  ieta_signed_offset = -ECAL_IETA_MAX_EXT;
  fillECAL_with_EEproj( hEvt_EE_energy[0], ieta_global_offset, ieta_signed_offset );

  // Fill middle part of ECAL(iphi,ieta) with the EB rechits.
  ieta_global_offset = 55;
  for ( EcalRecHitCollection::const_iterator iRHit = EBRecHitsH_->begin();
        iRHit != EBRecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ == 0. ) continue;
    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi() - 1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    // Fill vector for image
    ieta_global = ieta_ + ieta_global_offset;
    idx_ = ieta_global*EBDetId::MAX_IPHI + iphi_; 
    vECAL_energy_[idx_] = energy_;
    // Fill histogram for monitoring
    ieta_signed = ieta_;
    hECAL_energy->Fill( iphi_, ieta_signed, energy_ );

  } // EB

  // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
  ieta_global_offset = ECAL_IETA_MAX_EXT + EBDetId::MAX_IETA;
  ieta_signed_offset = EBDetId::MAX_IETA;
  fillECAL_with_EEproj( hEvt_EE_energy[1], ieta_global_offset, ieta_signed_offset );

} // fillECALstitched()
