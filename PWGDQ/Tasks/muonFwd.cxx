// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief This task shows how to access the Muons belong to a collision.
///        The association is made through the BC column (and in Run 3 may not be unique!).
///        Note that one has to subscribe to aod::Collisions const& to load
///        the relevant data even if you access the data itself through m.collision().
/// \author
/// \since

#include "CCDB/BasicCCDBManager.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"
#include <TH1F.h>
#include <THashList.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <algorithm>


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;



using MyDileptons = soa::Join<aod::Dileptons, aod::DileptonsVtx>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using MyMuonTracksWithCov = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;

// Iterate on muon using the collision iterator in the dq-analysis style
struct IterateMuons {
    HistogramRegistry registry{
    "registry",
    {
	    {"pT", "; pT; dileptons", {HistType::kTH1F, {{200, 0, 20}}}},
	    {"Tauz", "; Tauz; dileptons", {HistType::kTH1F, {{200, -0.01, 0.01}}}},
	    {"Rap", "; y; dileptons", {HistType::kTH1F, {{200, 2.5, 4.0}}}},
	    {"pT_Tauz", "; pT; Tauz; dileptons", {HistType::kTH2F, {{200, 0, 20}, {200, -0.01, 0.01}}}},
	    {"Rap_Tauz", "; y; Tauz; dileptons", {HistType::kTH2F, {{200, 2.5, 4.0}, {200, -0.01, 0.01}}}}
    }
    };

  void process(MyDileptons const& dileptons)
  {
//    LOGF(info, "Vertex = %f has %d muons", collision.posZ(), muons.size());
    for (auto& dilepton : dileptons) {
      LOGF(info, "  pT = %.2f", dilepton.pt());
      LOGF(info, "  Tauz = %.5f", dilepton.tauz());
      registry.fill(HIST("pT"), dilepton.pt());
      registry.fill(HIST("Rap"), dilepton.rap());
      registry.fill(HIST("Tauz"), dilepton.tauz());
      registry.fill(HIST("pT_Tauz"), dilepton.pt(), dilepton.tauz());
      registry.fill(HIST("Rap_Tauz"), dilepton.rap(), dilepton.tauz());
    }
  }
};

struct MuonsDCA {
    HistogramRegistry registry{
    "registry",
    {
	    {"pT", "; pT; muon", {HistType::kTH1F, {{200, 0, 20}}}},
	    {"Eta", "; Eta; muon", {HistType::kTH1F, {{500, -5.0, 5.0}}}},
	    {"Dca", "; dca; muon", {HistType::kTH1F, {{100, 0, 2}}}},
	    {"pT_dca", "; pT; dca; muon", {HistType::kTH2F, {{200, 0, 20}, {100, 0, 2}}}},
	    {"Eta_dca", "; Eta; dca; muon", {HistType::kTH2F, {{500, -5.0, 5.0}, {100, 0, 2}}}}
    }
    };

	void process(MyEventsVtxCov::iterator const& collision, MyMuonTracksWithCov const& muons){
		LOGF(info, "Vertex = %f has %d muons", collision.posZ(), muons.size());
	        for (auto& muon : muons) {
			    LOGF(info, "  pT = %.2f", muon.pt());
			    double chi2 = muon.chi2();
			    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
			    std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                          	    muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
	                            muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
			    SMatrix55 tcovs(v1.begin(), v1.end());
			    o2::track::TrackParCovFwd pars1{muon.z(), tpars, tcovs, chi2};
			    pars1.propagateToZlinear(collision.posZ());
			    auto dca = std::sqrt((pars1.getX()-collision.posX())*(pars1.getX()-collision.posX())+(pars1.getY()-collision.posY())*(pars1.getY()-collision.posY()));

			    registry.fill(HIST("pT"), muon.pt());
			    registry.fill(HIST("Eta"), muon.eta());
			    registry.fill(HIST("Dca"), dca);
			    registry.fill(HIST("pT_dca"), muon.pt(), dca);
			    registry.fill(HIST("Eta_dca"), muon.eta(), dca);



    		}
	}
};

// This uses the sparse matcher, so you also get BCs without a collision.
// You need to check with m.has_collision()

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    //adaptAnalysisTask<IterateMuonsExclusives>(cfgc), // currently does not work
    adaptAnalysisTask<IterateMuons>(cfgc),
    adaptAnalysisTask<MuonsDCA>(cfgc)
  };
}
