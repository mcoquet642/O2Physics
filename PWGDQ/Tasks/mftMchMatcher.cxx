// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copdyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file mftMchMatcher.cxx
/// \brief MFT-MCH matching tool for data preparation

#include "Common/DataModel/EventSelection.h"

#include "PWGDQ/Core/VarManager.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FwdTrackReAlignTables.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Field/MagneticField.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"


#include "GlobalTracking/MatchGlobalFwd.h"
#include "MFTTracking/Constants.h"
#include <TGeoGlobalMagField.h>

#include <TTree.h>

#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

using MyEvents = soa::Join<aod::Collisions, aod::EvSels, aod::MFTMults>;
using MyMuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMuonsMC = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels, aod::FwdTracksDCA, aod::FwdTrkCompColls>;
using MyMFTs = aod::MFTTracks;
using MyMFTCovariances = aod::MFTTracksCov;
using MyMFTsMC = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

using MyMuon = MyMuonsWithCov::iterator;
using MyMuonMC = MyMuonsMC::iterator;
using MyMFT = MyMFTs::iterator;
using MyMFTCovariance = MyMFTCovariances::iterator;


using namespace std;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix55Std = ROOT::Math::SMatrix<double, 5>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

namespace o2::aod
{
namespace fwdmatchcandidates
{
DECLARE_SOA_COLUMN(XMCH, xMCH, float);
DECLARE_SOA_COLUMN(YMCH, yMCH, float);
DECLARE_SOA_COLUMN(PhiMCH, phiMCH, float);
DECLARE_SOA_COLUMN(TanlMCH, tanlMCH, float);
DECLARE_SOA_COLUMN(InvQPtMCH, invQPtMCH, float);
DECLARE_SOA_COLUMN(TimeMCH, timeMCH, float);
DECLARE_SOA_COLUMN(TimeResMCH, timeResMCH, float);
DECLARE_SOA_COLUMN(Chi2MCH, chi2MCH, float);
DECLARE_SOA_COLUMN(Rabs, rabs, float);
DECLARE_SOA_COLUMN(PDCA, pDCA, float);
DECLARE_SOA_COLUMN(McMaskMCH, mcMaskMCH, int);

DECLARE_SOA_COLUMN(CXXMCH, cXXMCH, float);
DECLARE_SOA_COLUMN(CYYMCH, cYYMCH, float);
DECLARE_SOA_COLUMN(CPhiPhiMCH, cPhiPhiMCH, float);
DECLARE_SOA_COLUMN(CTglTglMCH, cTglTglMCH, float);
DECLARE_SOA_COLUMN(C1Pt1PtMCH, c1Pt1PtMCH, float);
DECLARE_SOA_COLUMN(CXYMCH, cXYMCH, float);
DECLARE_SOA_COLUMN(CPhiXMCH, cPhiXMCH, float);
DECLARE_SOA_COLUMN(CPhiYMCH, cPhiYMCH, float);
DECLARE_SOA_COLUMN(CTglXMCH, cTglXMCH, float);
DECLARE_SOA_COLUMN(CTglYMCH, cTglYMCH, float);
DECLARE_SOA_COLUMN(CTglPhiMCH, cTglPhiMCH, float);
DECLARE_SOA_COLUMN(C1PtXMCH, c1PtXMCH, float);
DECLARE_SOA_COLUMN(C1PtYMCH, c1PtYMCH, float);
DECLARE_SOA_COLUMN(C1PtPhiMCH, c1PtPhiMCH, float);
DECLARE_SOA_COLUMN(C1PtTglMCH, c1PtTglMCH, float);

DECLARE_SOA_COLUMN(XMFT, xMFT, float);
DECLARE_SOA_COLUMN(YMFT, yMFT, float);
DECLARE_SOA_COLUMN(PhiMFT, phiMFT, float);
DECLARE_SOA_COLUMN(TanlMFT, tanlMFT, float);
DECLARE_SOA_COLUMN(InvQPtMFT, invQPtMFT, float);
DECLARE_SOA_COLUMN(TimeMFT, timeMFT, float);
DECLARE_SOA_COLUMN(TimeResMFT, timeResMFT, float);
DECLARE_SOA_COLUMN(Chi2MFT, chi2MFT, float);
DECLARE_SOA_COLUMN(McMaskMFT, mcMaskMFT, int);
DECLARE_SOA_COLUMN(MftClusterSizesAndTrackFlags, mftClusterSizesAndTrackFlags, uint64_t);

DECLARE_SOA_COLUMN(CXXMFT, cXXMFT, float);
DECLARE_SOA_COLUMN(CYYMFT, cYYMFT, float);
DECLARE_SOA_COLUMN(CPhiPhiMFT, cPhiPhiMFT, float);
DECLARE_SOA_COLUMN(CTglTglMFT, cTglTglMFT, float);
DECLARE_SOA_COLUMN(C1Pt1PtMFT, c1Pt1PtMFT, float);
DECLARE_SOA_COLUMN(CXYMFT, cXYMFT, float);
DECLARE_SOA_COLUMN(CPhiYMFT, cPhiYMFT, float);
DECLARE_SOA_COLUMN(CPhiXMFT, cPhiXMFT, float);
DECLARE_SOA_COLUMN(CTglXMFT, cTglXMFT, float);
DECLARE_SOA_COLUMN(CTglYMFT, cTglYMFT, float);
DECLARE_SOA_COLUMN(CTglPhiMFT, cTglPhiMFT, float);
DECLARE_SOA_COLUMN(C1PtXMFT, c1PtXMFT, float);
DECLARE_SOA_COLUMN(C1PtYMFT, c1PtYMFT, float);
DECLARE_SOA_COLUMN(C1PtPhiMFT, c1PtPhiMFT, float);
DECLARE_SOA_COLUMN(C1PtTglMFT, c1PtTglMFT, float);

DECLARE_SOA_COLUMN(Chi2Glob, chi2Glob, float);
DECLARE_SOA_COLUMN(Chi2Match, chi2Match, float);
DECLARE_SOA_COLUMN(IsAmbig, isAmbig, bool);
DECLARE_SOA_COLUMN(MFTMult, mftMult, int);
DECLARE_SOA_COLUMN(DCAX, dcaX, float);
DECLARE_SOA_COLUMN(DCAY, dcaY, float);
DECLARE_SOA_COLUMN(McMaskGlob, mcMaskGlob, int);
}

DECLARE_SOA_TABLE(FwdMatchMLCandidates, "AOD", "FWDMLCAND",
		fwdmatchcandidates::XMCH,
		fwdmatchcandidates::YMCH,
		fwdmatchcandidates::PhiMCH,
		fwdmatchcandidates::TanlMCH,
		fwdmatchcandidates::InvQPtMCH,
		fwdmatchcandidates::TimeMCH,
		fwdmatchcandidates::TimeResMCH,
		fwdmatchcandidates::Chi2MCH,
		fwdmatchcandidates::PDCA,
		fwdmatchcandidates::Rabs,
		fwdmatchcandidates::CXXMCH,
		fwdmatchcandidates::CYYMCH,
		fwdmatchcandidates::CPhiPhiMCH, 
		fwdmatchcandidates::CTglTglMCH, 
		fwdmatchcandidates::C1Pt1PtMCH, 
		fwdmatchcandidates::CXYMCH,
		fwdmatchcandidates::CPhiYMCH,
		fwdmatchcandidates::CPhiXMCH,
		fwdmatchcandidates::CTglXMCH,
		fwdmatchcandidates::CTglYMCH,
		fwdmatchcandidates::CTglPhiMCH,
		fwdmatchcandidates::C1PtXMCH,
		fwdmatchcandidates::C1PtYMCH,
		fwdmatchcandidates::C1PtPhiMCH,
		fwdmatchcandidates::C1PtTglMCH,
		fwdmatchcandidates::XMFT,
		fwdmatchcandidates::YMFT,
		fwdmatchcandidates::PhiMFT,
		fwdmatchcandidates::TanlMFT,
		fwdmatchcandidates::InvQPtMFT,
		fwdmatchcandidates::TimeMFT,
		fwdmatchcandidates::TimeResMFT,
		fwdmatchcandidates::Chi2MFT,
		fwdmatchcandidates::MftClusterSizesAndTrackFlags,
		fwdmatchcandidates::CXXMFT,
		fwdmatchcandidates::CYYMFT,
		fwdmatchcandidates::CPhiPhiMFT, 
		fwdmatchcandidates::CTglTglMFT, 
		fwdmatchcandidates::C1Pt1PtMFT, 
		fwdmatchcandidates::CXYMFT,
		fwdmatchcandidates::CPhiYMFT,
		fwdmatchcandidates::CPhiXMFT,
		fwdmatchcandidates::CTglXMFT,
		fwdmatchcandidates::CTglYMFT,
		fwdmatchcandidates::CTglPhiMFT,
		fwdmatchcandidates::C1PtXMFT,
		fwdmatchcandidates::C1PtYMFT,
		fwdmatchcandidates::C1PtPhiMFT,
		fwdmatchcandidates::C1PtTglMFT,
		fwdmatchcandidates::Chi2Glob,
		fwdmatchcandidates::Chi2Match,
		fwdmatchcandidates::DCAX,
		fwdmatchcandidates::DCAY,
		fwdmatchcandidates::IsAmbig,
		fwdmatchcandidates::MFTMult,
		fwdmatchcandidates::McMaskMCH,
		fwdmatchcandidates::McMaskMFT,
		fwdmatchcandidates::McMaskGlob);
}

struct mftMchMatcher {
  Produces<o2::aod::FwdMatchMLCandidates> fwdMatchMLCandidates;
  ////   Variables for selecting muon tracks
  Configurable<float> fPMchLow{"cfgPMchLow", 0.0f, ""};
  Configurable<float> fPtMchLow{"cfgPtMchLow", 0.7f, ""};
  Configurable<float> fEtaMchLow{"cfgEtaMchLow", -4.0f, ""};
  Configurable<float> fEtaMchUp{"cfgEtaMchUp", -2.5f, ""};
  Configurable<float> fRabsLow{"cfgRabsLow", 17.6f, ""};
  Configurable<float> fRabsUp{"cfgRabsUp", 89.5f, ""};
  Configurable<float> fSigmaPdcaUp{"cfgPdcaUp", 6.f, ""};
  Configurable<float> fTrackChi2MchUp{"cfgTrackChi2MchUp", 5.f, ""};
  Configurable<float> fMatchingChi2MchMidUp{"cfgMatchingChi2MchMidUp", 999.f, ""};

  ////   Variables for selecting mft tracks
  Configurable<float> fEtaMFTLow{"cfgEtaMFTlow", -3.6f, ""};
  Configurable<float> fEtaMFTUp{"cfgEtaMFTup", -2.5f, ""};
  Configurable<float> fTrackChi2MFTUp{"cfgTrackChi2MFTUp", 10.f, ""};
  Configurable<float> fPtMFTLow{"cfgPtMFTLow", 0.1f, ""};

  ////   Variables for matching configuration
  Configurable<int> fMaxCandidates{"cfgMaxCandidates", 0, ""};

  Configurable<bool> fKeepBestMatch{"cfgKeepBestMatch", false, "Keep only the best match global muons in the skimming"};
  Configurable<float> fzMatching{"cfgzMatching", -77.5f, "Plane for MFT-MCH matching"};

  ////   Variables for ccdb
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  double mBzAtMftCenter{0};
  o2::globaltracking::MatchGlobalFwd mExtrap;

  int mRunNumber{0}; // needed to detect if the run changed and trigger update of magnetic field
  Service<o2::ccdb::BasicCCDBManager> ccdbManager;
  o2::ccdb::CcdbApi fCCDBApi;

  o2::parameters::GRPMagField* fGrpMag = nullptr;

  o2::globaltracking::MatchGlobalFwd mMatching;

  std::unordered_map<int64_t, int32_t> mftCovIndexes;

  std::map<uint32_t, bool> fBestMatch;

  HistogramRegistry registry{"registry", {}};

  template <class T, class C>
  bool pDCACut(const T& mchTrack, const C& collision, double nSigmaPDCA)
  {
    static const double sigmaPDCA23 = 80.;
    static const double sigmaPDCA310 = 54.;
    static const double relPRes = 0.0004;
    static const double slopeRes = 0.0005;

    double thetaAbs = TMath::ATan(mchTrack.rAtAbsorberEnd() / 505.) * TMath::RadToDeg();

    // propagate muon track to vertex
    auto trackConv = VarManager::FwdToTrackPar(mchTrack, mchTrack);
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(trackConv.getParameters());
    track.setZ(trackConv.getZ());
    track.setCovariances(trackConv.getCovariances());
    auto mchTrackAtVertex = mMatching.FwdtoMCH(track);
    o2::mch::TrackExtrap::extrapToVertex(mchTrackAtVertex, collision.posX(), collision.posY(), collision.posZ(), collision.covXX(), collision.covYY());

    // double pUncorr = mchTrack.p();
    double p = mchTrackAtVertex.p();

    double pDCA = mchTrack.pDca();
    double sigmaPDCA = (thetaAbs < 3) ? sigmaPDCA23 : sigmaPDCA310;
    double nrp = nSigmaPDCA * relPRes * p;
    double pResEffect = sigmaPDCA / (1. - nrp / (1. + nrp));
    double slopeResEffect = 535. * slopeRes * p;
    double sigmaPDCAWithRes = TMath::Sqrt(pResEffect * pResEffect + slopeResEffect * slopeResEffect);
    if (pDCA > nSigmaPDCA * sigmaPDCAWithRes) {
      return false;
    }

    return true;
  }

  template <class T, class C>
  bool IsGoodMuon(const T& mchTrack, const C& collision,
                  double chi2Cut,
                  double pCut,
                  double pTCut,
                  std::array<double, 2> etaCut,
                  std::array<double, 2> rAbsCut,
                  double nSigmaPdcaCut)
  {
    // chi2 cut
    if (mchTrack.chi2() > chi2Cut)
      return false;

    // momentum cut
    if (mchTrack.p() < pCut) {
      return false; // skip low-momentum tracks
    }

    // transverse momentum cut
    if (mchTrack.pt() < pTCut) {
      return false; // skip low-momentum tracks
    }

    // Eta cut
    double eta = mchTrack.eta();
    if ((eta < etaCut[0] || eta > etaCut[1])) {
      return false;
    }

    // RAbs cut
    double rAbs = mchTrack.rAtAbsorberEnd();
    if ((rAbs < rAbsCut[0] || rAbs > rAbsCut[1])) {
      return false;
    }

    // pDCA cut
    if (!pDCACut(mchTrack, collision, nSigmaPdcaCut)) {
      return false;
    }

    return true;
  }

  template <class T>
  bool IsGoodMFT(const T& mftTrack,
                  double chi2Cut,
                  double pTCut,
                  std::array<double, 2> etaCut)
  {
    // chi2 cut
    if (mftTrack.chi2() > chi2Cut)
      return false;

    // transverse momentum cut
    if (mftTrack.pt() < pTCut) {
      return false; // skip low-momentum tracks
    }

    // Eta cut
    double eta = mftTrack.eta();
    if ((eta < etaCut[0] || eta > etaCut[1])) {
      return false;
    }

    return true;
  }


  template <typename BC>
  void initCCDB(BC const& bc)
  {
    if (mRunNumber == bc.runNumber())
      return;

    fGrpMag = ccdbManager->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());

        if (fGrpMag != nullptr) {
          o2::base::Propagator::initFieldFromGRP(fGrpMag);
          VarManager::SetMagneticField(fGrpMag->getNominalL3Field());
          VarManager::SetupMuonMagField();
        }
    mRunNumber = bc.runNumber();
  }

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    ccdbManager->setURL(ccdburl);
    ccdbManager->setCaching(true);
    ccdbManager->setLocalObjectValidityChecking();
    fCCDBApi.init(ccdburl);
    mRunNumber = 0;

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      ccdbManager->get<TGeoManager>(geoPath);
    }
  }

  template <typename TMuons>
  void skimBestMuonMatches(TMuons const& muons)
  {
    std::unordered_map<int, std::pair<float, int>> mCandidates;
    for (const auto& muon : muons) {
      if (static_cast<int>(muon.trackType()) < 2) {
        auto muonID = muon.matchMCHTrackId();
        auto chi2 = muon.chi2MatchMCHMFT();
        if (mCandidates.find(muonID) == mCandidates.end()) {
          mCandidates[muonID] = {chi2, muon.globalIndex()};
        } else {
          if (chi2 < mCandidates[muonID].first) {
            mCandidates[muonID] = {chi2, muon.globalIndex()};
          }
        }
      }
    }
    for (auto& pairCand : mCandidates) {
      fBestMatch[pairCand.second.second] = true;
    }
  }

  void processMC(MyEvents const& collisions,
                 aod::BCsWithTimestamps const& bcs,
                 MyMuonsMC const& muonTracks,
                 MyMFTsMC const& mftTracks,
                 MyMFTCovariances const& mftCovs,
                 aod::McParticles const& /*mcParticles*/)
  {
    if (bcs.size() > 0) {
      auto bc = bcs.begin();
      initCCDB(bc);
      VarManager::SetMatchingPlane(fzMatching.value);
    }
    if (fKeepBestMatch) {
	    skimBestMuonMatches(muonTracks);
    }

    mftCovIndexes.clear();
    for (auto& mftTrackCov : mftCovs) {
      mftCovIndexes[mftTrackCov.matchMFTTrackId()] = mftTrackCov.globalIndex();
    }

    fwdMatchMLCandidates.reserve(muonTracks.size());

    for (auto muon : muonTracks) {
      // only consider global MFT-MCH-MID matches
      if (static_cast<int>(muon.trackType()) != 0) {
        continue;
      }

      if (!muon.has_collision()) {
        continue;
      }

      if (fKeepBestMatch) {
        if (fBestMatch.find(muon.globalIndex()) == fBestMatch.end()) {
          continue;
        }
      }

      const auto& collision = collisions.rawIteratorAt(muon.collisionId());
      auto bc_coll = collision.bc_as<aod::BCsWithTimestamps>();

      auto muontrack = muon.template matchMCHTrack_as<MyMuonsMC>();
      auto mfttrack = muon.template matchMFTTrack_as<MyMFTsMC>();
      auto const& mfttrackcov = mftCovs.rawIteratorAt(mftCovIndexes[mfttrack.globalIndex()]);

      auto muonTime = muontrack.trackTime() + bc_coll.globalBC()*o2::constants::lhc::LHCBunchSpacingNS;
      auto mftTime = mfttrack.trackTime() + bc_coll.globalBC()*o2::constants::lhc::LHCBunchSpacingNS;

      o2::track::TrackParCovFwd mftprop = VarManager::FwdToTrackPar(mfttrack, mfttrackcov);
      o2::track::TrackParCovFwd muonprop = VarManager::FwdToTrackPar(muontrack, muontrack);
      if (fzMatching.value < 0.) {
        mftprop = VarManager::PropagateFwd(mfttrack, mfttrackcov, fzMatching.value);
        muonprop = VarManager::PropagateMuon(muontrack, collision, VarManager::kToMatching);
      }
      auto muonpropCov = muonprop.getCovariances();
      auto mftpropCov = mftprop.getCovariances();

      if (!IsGoodMuon(muontrack, collision, fTrackChi2MchUp, fPMchLow, fPtMchLow, {fEtaMchLow, fEtaMchUp}, {fRabsLow, fRabsUp}, fSigmaPdcaUp)){
	      continue;
      }


      if (!IsGoodMFT(mfttrack, fTrackChi2MFTUp, fPtMFTLow, {fEtaMFTLow, fEtaMFTUp})){
	      continue;
      }

      bool IsAmbig = (muon.compatibleCollIds().size() != 1);
      int MFTMult = collision.mftNtracks();

    fwdMatchMLCandidates(
		    muonprop.getX(),
		    muonprop.getY(),
		    muonprop.getPhi(),
		    muonprop.getTgl(),
		    muonprop.getInvQPt(),
		    muonTime,
		    muontrack.trackTimeRes(),
		    muontrack.chi2(),
		    muontrack.pDca(),
		    muontrack.rAtAbsorberEnd(),
		    muonpropCov(0, 0),
		    muonpropCov(1, 1),
		    muonpropCov(2, 2),
		    muonpropCov(3, 3),
		    muonpropCov(4, 4),
		    muonpropCov(1, 0),
		    muonpropCov(2, 0),
		    muonpropCov(2, 1),
		    muonpropCov(3, 0),
		    muonpropCov(3, 1),
		    muonpropCov(3, 2),
		    muonpropCov(4, 0),
		    muonpropCov(4, 1),
		    muonpropCov(4, 2),
		    muonpropCov(4, 3),
		    mftprop.getX(),
		    mftprop.getY(),
		    mftprop.getPhi(),
		    mftprop.getTgl(),
		    mftprop.getInvQPt(),
		    mftTime,
		    mfttrack.trackTimeRes(),
		    mfttrack.chi2(),
		    mfttrack.mftClusterSizesAndTrackFlags(),
		    mftpropCov(0, 0),
		    mftpropCov(1, 1),
		    mftpropCov(2, 2),
		    mftpropCov(3, 3),
		    mftpropCov(4, 4),
		    mftpropCov(1, 0),
		    mftpropCov(2, 0),
		    mftpropCov(2, 1),
		    mftpropCov(3, 0),
		    mftpropCov(3, 1),
		    mftpropCov(3, 2),
		    mftpropCov(4, 0),
		    mftpropCov(4, 1),
		    mftpropCov(4, 2),
		    mftpropCov(4, 3),
		    muon.chi2(),
		    muon.chi2MatchMCHMFT(),
		    muon.fwdDcaX(),
		    muon.fwdDcaY(),
		    IsAmbig,
		    MFTMult,
		    muontrack.mcMask(),
		    mfttrack.mcMask(),
		    muon.mcMask()
		    );



    }

  }

  PROCESS_SWITCH(mftMchMatcher, processMC, "process_MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<mftMchMatcher>(cfgc)};
};

