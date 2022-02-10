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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "TDatabasePDG.h"
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

constexpr int each = 1;



namespace
{
inline auto normalize(float a)
{
  if (a > M_PI) {
    a -= 2 * M_PI;
  }
  if (a < -M_PI) {
    a += 2 * M_PI;
  }
  return a;
}
} // namespace



struct analysemfttracks {
  
   HistogramRegistry registry{
    "registry",
    {
      {"Zvtx_Colln", "; Z_{vtx}; events", {HistType::kTH1F, {{402, -20.1, 20.1}}}},       //
      {"hNMFTTracks", ";Ntracks_{MFT}; events", {HistType::kTH1F, {{301, -0.5, 300.5}}}},       //
       {"hNMFTTracks_1", ";Ntracks_{MFT}; events", {HistType::kTH1F, {{301, -0.5, 300.5}}}},       //
       {"hNMFTTracks_2", ";Ntracks_{MFT}; events", {HistType::kTH1F, {{301, -0.5, 300.5}}}},       //
       {"hNMFTTracks_pos", ";Ntracks_{MFT}; events", {HistType::kTH1F, {{301, -0.5, 300.5}}}},       //
       {"hNMFTTracks_neg", ";Ntracks_{MFT}; events", {HistType::kTH1F, {{301, -0.5, 300.5}}}},       //

      {"NtrkZvtx", "; Z_{vtx};  N_{trk};events", {HistType::kTH2F, {{402, -20.1, 20.1}, {101, -0.5, 100.5}}}},       //
      {"hEtaZvtx", "; Z_{vtx};  #eta;events", {HistType::kTH2F, {{402, -20.1, 20.1},  {20, -4., -2.}}}},       //
       

      {"hEtaMFT", "#eta; N_{trk}; events", {HistType::kTH1F, {{200, -4., -2.}}}},     // -2.5 to -3.6
      {"hPhiMFT", "#phi; N_{trk}; events", {HistType::kTH1F, {{300, -1.5*M_PI, 1.5 *M_PI}}}},                           
      {"hPhiEtaMFT", "; #varphi; #eta; tracks", {HistType::kTH2F, {{150, -1.5*M_PI, 1.5 *M_PI}, {20, -4., -2.}}}}, //

      {"hXYMFT", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hZMFT", "; #z; tracks", {HistType::kTH1F, {{1000, -100, 0}}}},

      {"hXYMFT_Disk5_b", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk5_a", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk4_b", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk4_a", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk3_b", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk3_a", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk2_b", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk2_a", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk1_b", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hXYMFT_Disk1_a", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
     		
      {"hClustersMFT", "#eta; N_{trk}; events", {HistType::kTH1F, {{20, -0.5, 19.5}}}},     // -2.5 to -3.6
      
      {"Zvtx_MC_Generated_MFT", "; Z_{vtx}; events", {HistType::kTH1F, {{402, -20.1, 20.1}}}},       //
      {"Size_MC_Generated_MFT", "; Number; events", {HistType::kTH1F, {{301, -0.5, 300.5}}}},
      {"PDG_MC_Generated_MFT", ";#phi_{MC}", {HistType::kTH1F, {{6000,-0.5, 5999.5}}}},
      {"eta_MC_Generated_MFT", ";#eta_{MC}", {HistType::kTH1F, {{600, -4, 2}}}},
      {"phi_MC_Generated_MFT", ";#phi_{MC}", {HistType::kTH1F, {{100, 0, 2 * M_PI}}}},

      {"Zvtx_MC_MFT_JoinedTable", "; Z_{vtx}; events", {HistType::kTH1F, {{402, -20.1, 20.1}}}},       //
      {"Zvtx_MC_MFT_JoinedTable_mcCollision", "; Z_{vtx}; events", {HistType::kTH1F, {{402, -20.1, 20.1}}}},       //
      {"Size_MC_MFT_JoinedTable", "; Number; events", {HistType::kTH1F, {{301, -0.5, 300.5}}}},
      {"eta_MC_MFT_JoinedTable_mcParticle", ";#eta_{MC} (mcParticle)", {HistType::kTH1F,  {{200, -4., -2.}}}},
      {"eta_MC_MFT_JoinedTable", ";#eta_{MC} ", {HistType::kTH1F,  {{200, -4., -2.}}}},
      {"etaDiff_MC_MFT_JoinedTable", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
      {"phiDiff_MC_MFT_JoinedTable", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH1F, {{100, -M_PI, M_PI}}}},
      {"PDG_MC_MFT_JoinedTable", ";#phi_{MC}", {HistType::kTH1F, {{6000,-0.5, 5999.5}}}}
          }                                                                          //
  };


  void processMFT(o2::aod::Collision const& collision, o2::aod::MFTTracks const& tracks)//, o2::aod::SingleTrack const& stracks)
  {
      auto z = collision.posZ();
      registry.fill(HIST("Zvtx_Colln"), z);
      registry.fill(HIST("hNMFTTracks"), tracks.size());
      if(z> 10 || z< -10) registry.fill(HIST("hNMFTTracks_1"), tracks.size());
      if(z<= 10 && z>= -10) registry.fill(HIST("hNMFTTracks_2"), tracks.size());
	
      if(z<= 0 ) registry.fill(HIST("hNMFTTracks_neg"), tracks.size());
      if( z>= 0) registry.fill(HIST("hNMFTTracks_pos"), tracks.size());
     
      registry.fill(HIST("NtrkZvtx"), z, tracks.size());

	
    for (auto& track : tracks) {
      registry.fill(HIST("hPhiEtaMFT"), track.phi(), track.eta());
      registry.fill(HIST("hEtaMFT"), track.eta());
      registry.fill(HIST("hPhiMFT"), track.phi());

      registry.fill(HIST("hXYMFT"), track.x(), track.y());
      registry.fill(HIST("hZMFT"), track.z());

      registry.fill(HIST("hClustersMFT"), track.nClusters());

      registry.fill(HIST("hEtaZvtx"), z,  track.eta());
	
      if(track.z()< -77 && track.z() > -78)  registry.fill(HIST("hXYMFT_Disk5_b"), track.x(), track.y());
      if(track.z()< -75 && track.z() > -77)  registry.fill(HIST("hXYMFT_Disk5_a"), track.x(), track.y());

      if(track.z()< -69 && track.z() > -70)  registry.fill(HIST("hXYMFT_Disk4_b"), track.x(), track.y());
      if(track.z()< -67 && track.z() > -68)  registry.fill(HIST("hXYMFT_Disk4_a"), track.x(), track.y());

      if(track.z()< -53 && track.z() > -54)  registry.fill(HIST("hXYMFT_Disk3_b"), track.x(), track.y());
      if(track.z()< -52 && track.z() > -53)  registry.fill(HIST("hXYMFT_Disk3_a"), track.x(), track.y());

      if(track.z()< -49 && track.z() > -50.5)  registry.fill(HIST("hXYMFT_Disk2_b"), track.x(), track.y());
      if(track.z()< -48 && track.z() > -49)  registry.fill(HIST("hXYMFT_Disk2_a"), track.x(), track.y());

      if(track.z()< -46 && track.z() > -47)  registry.fill(HIST("hXYMFT_Disk1_b"), track.x(), track.y());
      if(track.z()< -45 && track.z() > -46)  registry.fill(HIST("hXYMFT_Disk1_a"), track.x(), track.y());
 
     
      // LOGP(info, "Track {} has x = {}, y = {}, z = {}", track.index(), track.x(),  track.y(),  track.z());
    }
    // LOGP(info, "============================");

  }
  PROCESS_SWITCH(analysemfttracks, processMFT, "Process MFT info", true);


   void processGen(aod::McCollisions::iterator const& collisionMC, aod::McParticles const& mcParticles)
  {
     registry.fill(HIST("Zvtx_MC_Generated_MFT"), collisionMC.posZ());
     registry.fill(HIST("Size_MC_Generated_MFT"), mcParticles.size());

    for (auto& particle : mcParticles) {
      //LOGF(info, " Particle MC PDG =  %d, ", abs(particle.pdgCode()) );
      registry.fill(HIST("PDG_MC_Generated_MFT"), abs(particle.pdgCode()) );
      registry.fill(HIST("eta_MC_Generated_MFT"), particle.eta());
      registry.fill(HIST("phi_MC_Generated_MFT"), particle.phi());
    }
  }

  PROCESS_SWITCH(analysemfttracks, processGen, "Process gen level MFT", true);


  void processResolutionMFT(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::MFTTracks, aod::McMFTTrackLabels> const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
     //  LOGF(info, "vtx-z (data) = %f | vtx-z (MC) = %f", collision.posZ(), collision.mcCollision().posZ());

       registry.fill(HIST("Zvtx_MC_MFT_JoinedTable"), collision.posZ());
       registry.fill(HIST("Zvtx_MC_MFT_JoinedTable_mcCollision"),collision.mcCollision().posZ());
   
       registry.fill(HIST("Size_MC_MFT_JoinedTable"), tracks.size());

     for (auto& track : tracks) {
      auto particle = track.mcParticle();
      if (particle.isPhysicalPrimary()){
       // LOGF(info, " Particle Joined PDG =  %d, Eta Gen = %f, Eta Reco = %f, Eta Diff = %f", abs(particle.pdgCode()),  track.eta(), track.mcParticle().eta() , track.mcParticle().eta() - track.eta());
        registry.fill(HIST("PDG_MC_MFT_JoinedTable"), abs(particle.pdgCode()) );
         registry.fill(HIST("eta_MC_MFT_JoinedTable_mcParticle"), track.mcParticle().eta());
      	 registry.fill(HIST("eta_MC_MFT_JoinedTable"),  track.eta());
        registry.fill(HIST("etaDiff_MC_MFT_JoinedTable"), track.mcParticle().eta() - track.eta());
      	registry.fill(HIST("phiDiff_MC_MFT_JoinedTable"), normalize(track.mcParticle().phi() - track.phi()));
      }
    }
  }
  PROCESS_SWITCH(analysemfttracks, processResolutionMFT, "Process reco/gen matching", true);


};


struct IterateMuons {

    HistogramRegistry registry{
    "registry",
    {
       {"Zvtx_FwdTracks_Muons", "; Z_{vtx}; events", {HistType::kTH1F, {{402, -20.1, 20.1}}}}, 
       {"XYvtx_FwdTracks_Muons", "; X_{vtx}; Y_{vtx};  events", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}}, 
      {"Size_FwdTracks_Muons", "; muons.size(); events", {HistType::kTH1F, {{51, -0.5, 50.5}}}},

      {"Size_FwdTracks_Muons_Type0", "; muons.size(); events", {HistType::kTH1F, {{51, -0.5, 50.5}}}},
      {"Size_FwdTracks_Muons_Type2", "; muons.size(); events", {HistType::kTH1F, {{51, -0.5, 50.5}}}},
      {"Size_FwdTracks_Muons_Type4", "; muons.size(); events", {HistType::kTH1F, {{51, -0.5, 50.5}}}},

      {"Pt_FwdTracks_Muons", "; Pt; events", {HistType::kTH1F, {{400, 0, 40}}}},  
      {"MuonTrackType", "; Type; events", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
      {"Chi2MCHMIDMatch", "; Type; events", {HistType::kTH1F, {{100, -2.5, 7.5 }}}},
      {"Chi2MCHMFTMatch", "; Type; events", {HistType::kTH1F, {{100, -2.5, 7.5}}}},


      {"hXY_FwdTracks_Type0", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hZ_FwdTracks_Type0", "; #z; tracks", {HistType::kTH1F, {{1000, -100, 0}}}},
      {"hX_FwdTracks_Type0", "; #x; tracks", {HistType::kTH1F, {{400, -20, 20}}}},
      {"hY_FwdTracks_Type0", "; #y; tracks", {HistType::kTH1F, {{400, -20, 20}}}},
      {"hTime_FwdTracks_Type0", "; #Time; tracks", {HistType::kTH1F, {{520, -2, 50}}}},
      {"Pt_FwdTracks_Type0", "; Pt; events", {HistType::kTH1F, {{400, 0, 40}}}},
      {"XYDiff_Trial_FwdTracks_Type0", "; XYDiff; tracks", {HistType::kTH1F, {{400, 0, 40}}}}, 

      {"hXY_FwdTracks_Type2", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hZ_FwdTracks_Type2", "; #z; tracks", {HistType::kTH1F, {{1000, -100, 0}}}},
      {"hX_FwdTracks_Type2", "; #x; tracks", {HistType::kTH1F, {{400, -20, 20}}}},
      {"hY_FwdTracks_Type2", "; #y; tracks", {HistType::kTH1F, {{400, -20, 20}}}},
      {"hTime_FwdTracks_Type2", "; #Time; tracks", {HistType::kTH1F, {{520, -2, 50}}}},
      {"Pt_FwdTracks_Type2", "; Pt; events", {HistType::kTH1F, {{400, 0, 40}}}},
      {"XYDiff_Trial_FwdTracks_Type2", "; XYDiff; tracks", {HistType::kTH1F, {{400, 0, 40}}}}, 


      {"hXY_FwdTracks_Type4", "; #x; #y; tracks", {HistType::kTH2F, {{400, -20, 20}, {400, -20, 20}}}},
      {"hZ_FwdTracks_Type4", "; #z; tracks", {HistType::kTH1F, {{1500, -100, 50}}}},
      {"hX_FwdTracks_Type4", "; #x; tracks", {HistType::kTH1F, {{400, -20, 20}}}},
      {"hY_FwdTracks_Type4", "; #y; tracks", {HistType::kTH1F, {{400, -20, 20}}}},
      {"hTime_FwdTracks_Type4", "; #Time; tracks", {HistType::kTH1F, {{400, -20, 20}}}},
      {"Pt_FwdTracks_Type4", "; Pt; tracks", {HistType::kTH1F, {{400, 0, 40}}}},
      {"XYDiff_Trial_FwdTracks_Type4", "; XYDiff; tracks", {HistType::kTH1F, {{400, 0, 40}}}}, 	
	
      {"Zvtx_Generated_MC", "; Z_{vtx}; events", {HistType::kTH1F, {{402, -20.1, 20.1}}}}, 
      {"Size_Generated_MC", "; mcParticles.size(); events", {HistType::kTH1F, {{51, -0.5, 50.5}}}},
      {"eta_Generated_MC", ";#eta_{MC}", {HistType::kTH1F, {{100, -2, 2}}}},
      {"phi_Generated_MC", ";#phi_{MC}", {HistType::kTH1F, {{100, 0, 2 * M_PI}}}}, ///	
      {"PDG_Generated_MC", ";#phi_{MC}", {HistType::kTH1F, {{6000,-0.5, 5999.5}}}},


      {"Zvtx_MC_Fwd_JoinedTable", "; Z_{vtx}; events", {HistType::kTH1F, {{402, -20.1, 20.1}}}},       //
      {"Zvtx_MC_Fwd_JoinedTable_mcCollision", "; Z_{vtx}; events", {HistType::kTH1F, {{402, -20.1, 20.1}}}},       //
      {"Size_MC_Fwd_JoinedTable", "; Number; events", {HistType::kTH1F, {{51, -0.5, 50.5}}}},
      {"MuonTrackType_MC_Fwd_JoinedTable", "; Type; events", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
      // {"MuonTrackType_Generated_MC_Fwd_JoinedTable", "; Type; events", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
       {"PDG_MC_Fwd_JoinedTable", ";#phi_{MC}", {HistType::kTH1F, {{6000,-0.5, 5999.5}}}},
	{"PDG_MC_Fwd_JoinedTable_Type0", ";#phi_{MC}", {HistType::kTH1F, {{6000,-0.5, 5999.5}}}},
	{"PDG_MC_Fwd_JoinedTable_Type2", ";#phi_{MC}", {HistType::kTH1F, {{6000,-0.5, 5999.5}}}},
	{"PDG_MC_Fwd_JoinedTable_Type4", ";#phi_{MC}", {HistType::kTH1F, {{6000,-0.5, 5999.5}}}},

	{"PDG_Type_2D_MC_Fwd_JoinedTable", ";PDG, Type", {HistType::kTH2F, {{6000,-0.5, 5999.5},{6, -0.5, 5.5}}}},

        {"etaDiff_MC_Fwd_JoinedTable", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
        {"phiDiff_MC_Fwd_JoinedTable", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH1F, {{100, -M_PI, M_PI}}}},
        {"pTDiff_MC_Fwd_JoinedTable", ";#pT_{MC} - #pT_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
	 {"etaDiff_MC_Fwd_JoinedTable_Type0", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
        {"phiDiff_MC_Fwd_JoinedTable_Type0", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH1F, {{100, -M_PI, M_PI}}}},
        {"pTDiff_MC_Fwd_JoinedTable_Type0", ";#pT_{MC} - #pT_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
	 {"etaDiff_MC_Fwd_JoinedTable_Type2", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
        {"phiDiff_MC_Fwd_JoinedTable_Type2", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH1F, {{100, -M_PI, M_PI}}}},
        {"pTDiff_MC_Fwd_JoinedTable_Type2", ";#pT_{MC} - #pT_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
	{"etaDiff_MC_Fwd_JoinedTable_Type4", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
        {"phiDiff_MC_Fwd_JoinedTable_Type4", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH1F, {{100, -M_PI, M_PI}}}},
        {"pTDiff_MC_Fwd_JoinedTable_Type4", ";#pT_{MC} - #pT_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
     // {"etaRec", ";#eta_{Rec}", {HistType::kTH1F, {{100, -2, 2}}}},
     // {"phiRec", ";#phi_{Rec}", {HistType::kTH1F, {{100, 0, 2 * M_PI}}}},

	{"eta_2D_MC_Fwd_JoinedTable", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH2F, {{300, -5, -2},{300, -5, -2}}}},
        {"phi_2D_MC_Fwd_JoinedTable", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH2F, {{300, -1.5*M_PI, 2.5*M_PI},{300, -1.5*M_PI, 2.5*M_PI}}}},
        {"pT_2D_MC_Fwd_JoinedTable", ";#pT_{MC} - #pT_{Rec}", {HistType::kTH2F, {{200, 0, 20},{200, 0, 20}}}},

	{"eta_2D_MC_Fwd_JoinedTable_Type0", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH2F, {{300, -5, -2},{300, -5, -2}}}},
        {"phi_2D_MC_Fwd_JoinedTable_Type0", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH2F, {{300, -1.5*M_PI, 2.5*M_PI},{300, -1.5*M_PI, 2.5*M_PI}}}},
        {"pT_2D_MC_Fwd_JoinedTable_Type0", ";#pT_{MC} - #pT_{Rec}", {HistType::kTH2F, {{200, 0, 20},{200, 0, 20}}}},

	{"eta_2D_MC_Fwd_JoinedTable_Type2", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH2F, {{300, -5, -2},{300, -5, -2}}}},
        {"phi_2D_MC_Fwd_JoinedTable_Type2", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH2F, {{300, -1.5*M_PI, 2.5*M_PI},{300, -1.5*M_PI, 2.5*M_PI}}}},
        {"pT_2D_MC_Fwd_JoinedTable_Type2", ";#pT_{MC} - #pT_{Rec}", {HistType::kTH2F, {{200, 0, 20},{200, 0, 20}}}},

	{"eta_2D_MC_Fwd_JoinedTable_Type4", ";#eta_{MC} - #eta_{Rec}", {HistType::kTH2F, {{300, -5, -2},{300, -5, -2}}}},
        {"phi_2D_MC_Fwd_JoinedTable_Type4", ";#phi_{MC} - #phi_{Rec}", {HistType::kTH2F, {{300, -1.5*M_PI, 2.5*M_PI},{300, -1.5*M_PI, 2.5*M_PI}}}},
        {"pT_2D_MC_Fwd_JoinedTable_Type4", ";#pT_{MC} - #pT_{Rec}", {HistType::kTH2F, {{200, 0, 20},{200, 0, 20}}}}

    }
    };


  void process(aod::Collisions::iterator const& collision, aod::FwdTracks const& muons)
  {
    //LOGF(info, "Vertex = %f has %d muons", collision.posZ(), muons.size());
     registry.fill(HIST("Zvtx_FwdTracks_Muons"),  collision.posZ());
     registry.fill(HIST("XYvtx_FwdTracks_Muons"),  collision.posX(), collision.posY());
     registry.fill(HIST("Size_FwdTracks_Muons"), muons.size());

     auto iSize0=0, iSize2=0, iSize4=0;
     auto xyDiff=0.0;

    for (auto& muon : muons) {
       //LOGF(info, "  pT = %.2f, chi2 MCH-MID = %.2f, chi2 MCH-MFT = %.2f", muon.pt(), muon.chi2MatchMCHMID(), muon.chi2MatchMCHMFT()); // chi2 MCH-MID = -1.00, chi2 MCH-MFT = -1.00
       //LOGF(info, "  Track Type = %u", muon.trackType()); 
       registry.fill(HIST("Pt_FwdTracks_Muons"), muon.pt());
       registry.fill(HIST("MuonTrackType"), muon.trackType()); // 4 = o2::aod::fwdtrack::MCHStandaloneTrack
       registry.fill(HIST("Chi2MCHMIDMatch"), muon.chi2MatchMCHMID());	
       registry.fill(HIST("Chi2MCHMFTMatch"), muon.chi2MatchMCHMFT());
       if(muon.trackType() ==4) iSize4++;
       if(muon.trackType() ==2) iSize2++;
       if(muon.trackType() ==0) iSize0++;
//http://aliceo2group.github.io/AliceO2/d2/d12/DataTypes_8h.html
//http://aliceo2group.github.io/AliceO2/d2/d12/DataTypes_8h_source.html#l00057

	xyDiff =  sqrt( pow((collision.posX() - muon.x()), 2.0) +  pow((collision.posY() - muon.y()), 2.0)   );
	 
 	if(muon.trackType() ==0){
		registry.fill(HIST("hXY_FwdTracks_Type0"), muon.x(), muon.y());
		registry.fill(HIST("hZ_FwdTracks_Type0"), muon.z());
		registry.fill(HIST("hX_FwdTracks_Type0"), muon.x());
		registry.fill(HIST("hY_FwdTracks_Type0"), muon.y());
		registry.fill(HIST("hTime_FwdTracks_Type0"), muon.trackTime());
		registry.fill(HIST("Pt_FwdTracks_Type0"), muon.pt());
		registry.fill(HIST("XYDiff_Trial_FwdTracks_Type0"), xyDiff);
		//LOGF(info, "  ColX  = %.2f, muon x = %.2f,            ColY  = %.2f, muon y = %.2f,   xyDiff = %.2f", collision.posX() , muon.x(), collision.posY() , muon.y(),   xyDiff); // chi2 MCH-MID = -1.00, chi2 MCH-MFT = -1.00             
        
	}
	if(muon.trackType() ==2){
		registry.fill(HIST("hXY_FwdTracks_Type2"), muon.x(), muon.y());
		registry.fill(HIST("hZ_FwdTracks_Type2"), muon.z());
		registry.fill(HIST("hX_FwdTracks_Type2"), muon.x());
		registry.fill(HIST("hY_FwdTracks_Type2"), muon.y());
		registry.fill(HIST("hTime_FwdTracks_Type2"), muon.trackTime());
		registry.fill(HIST("Pt_FwdTracks_Type2"), muon.pt());
		registry.fill(HIST("XYDiff_Trial_FwdTracks_Type2"), xyDiff);
		//LOGF(info, "  ColX  = %.2f, muon x = %.2f,            ColY  = %.2f, muon y = %.2f,   xyDiff = %.2f", collision.posX() , muon.x(), collision.posY() , muon.y(),   xyDiff); // chi2 MCH-MID = -1.00, chi2 MCH-MFT = -1.00
               
        
	}
	if(muon.trackType() ==4){
		registry.fill(HIST("hXY_FwdTracks_Type4"), muon.x(), muon.y());
		registry.fill(HIST("hZ_FwdTracks_Type4"), muon.z());
		registry.fill(HIST("hX_FwdTracks_Type4"), muon.x());
		registry.fill(HIST("hY_FwdTracks_Type4"), muon.y());
		registry.fill(HIST("hTime_FwdTracks_Type4"), muon.trackTime());
 		registry.fill(HIST("Pt_FwdTracks_Type4"), muon.pt());
                registry.fill(HIST("XYDiff_Trial_FwdTracks_Type4"), xyDiff);
	}		
       	
       	
    }
    // LOGF(info, " iSize = %d , muons.size() = %d ",iSize0, muons.size());
     registry.fill(HIST("Size_FwdTracks_Muons_Type0"), iSize0);
     registry.fill(HIST("Size_FwdTracks_Muons_Type2"), iSize2);
     registry.fill(HIST("Size_FwdTracks_Muons_Type4"), iSize4);

  }
   
  void processGeneratedMC(aod::McCollisions::iterator const& collisionMC, aod::McParticles const& mcParticles)
  {
     registry.fill(HIST("Zvtx_Generated_MC"), collisionMC.posZ());
     registry.fill(HIST("Size_Generated_MC"), mcParticles.size());

    for (auto& particle : mcParticles) {
      registry.fill(HIST("PDG_Generated_MC"), abs(particle.pdgCode()) );
      registry.fill(HIST("eta_Generated_MC"), particle.eta());
      registry.fill(HIST("phi_Generated_MC"), particle.phi());
    }
  }

  PROCESS_SWITCH(IterateMuons, processGeneratedMC, "Process gen level MC ", true);



   void processResolutionFwdTacks(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::FwdTracks, aod::McFwdTrackLabels> const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
     //LOGF(info, "vtx-z (data) = %f | vtx-z (MC) = %f", collision.posZ(), collision.mcCollision().posZ());
    //  LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisions.size());
  
     registry.fill(HIST("Zvtx_MC_Fwd_JoinedTable"), collision.posZ());
     registry.fill(HIST("Zvtx_MC_Fwd_JoinedTable_mcCollision"),collision.mcCollision().posZ());
     registry.fill(HIST("Size_MC_Fwd_JoinedTable"), tracks.size());

    for (auto& track : tracks) {
	
      auto particle = track.mcParticle();
      if (particle.isPhysicalPrimary()){
        //std::cout<<"  Is Primary "<<endl;
       // LOGF(info, " Particle pid =  %d, ", particle.pdgCode());
        registry.fill(HIST("PDG_MC_Fwd_JoinedTable"), abs(particle.pdgCode()) );
        registry.fill(HIST("PDG_Type_2D_MC_Fwd_JoinedTable"), abs(particle.pdgCode()), track.trackType());
        registry.fill(HIST("MuonTrackType_MC_Fwd_JoinedTable"), track.trackType());
        // This is wrong => registry.fill(HIST("MuonTrackType_Generated_MC_Fwd_JoinedTable"), particle.eta());
        if(track.trackType() ==0){
		registry.fill(HIST("PDG_MC_Fwd_JoinedTable_Type0"), abs(particle.pdgCode()) );
		registry.fill(HIST("etaDiff_MC_Fwd_JoinedTable_Type0"), particle.eta() - track.eta());
	      	registry.fill(HIST("pTDiff_MC_Fwd_JoinedTable_Type0"), particle.pt() - track.pt());
	      	registry.fill(HIST("phiDiff_MC_Fwd_JoinedTable_Type0"), normalize(particle.phi() - track.phi()));

		registry.fill(HIST("eta_2D_MC_Fwd_JoinedTable_Type0"), particle.eta(), track.eta());
	      	registry.fill(HIST("pT_2D_MC_Fwd_JoinedTable_Type0"), particle.pt(), track.pt());
	      	registry.fill(HIST("phi_2D_MC_Fwd_JoinedTable_Type0"), normalize(particle.phi()), track.phi());
	}
        if(track.trackType() == 2){
		registry.fill(HIST("PDG_MC_Fwd_JoinedTable_Type2"), abs(particle.pdgCode()) );
		registry.fill(HIST("etaDiff_MC_Fwd_JoinedTable_Type2"), particle.eta() - track.eta());
	      	registry.fill(HIST("pTDiff_MC_Fwd_JoinedTable_Type2"), particle.pt() - track.pt());
	      	registry.fill(HIST("phiDiff_MC_Fwd_JoinedTable_Type2"), normalize(particle.phi() - track.phi()));

		registry.fill(HIST("eta_2D_MC_Fwd_JoinedTable_Type2"), particle.eta(), track.eta());
	      	registry.fill(HIST("pT_2D_MC_Fwd_JoinedTable_Type2"), particle.pt(), track.pt());
	      	registry.fill(HIST("phi_2D_MC_Fwd_JoinedTable_Type2"), normalize(particle.phi()), track.phi());
	}
        if(track.trackType() == 4){
		registry.fill(HIST("PDG_MC_Fwd_JoinedTable_Type4"), abs(particle.pdgCode()) );
		registry.fill(HIST("etaDiff_MC_Fwd_JoinedTable_Type4"), particle.eta() - track.eta());
	      	registry.fill(HIST("pTDiff_MC_Fwd_JoinedTable_Type4"), particle.pt() - track.pt());
	      	registry.fill(HIST("phiDiff_MC_Fwd_JoinedTable_Type4"), normalize(particle.phi() - track.phi()));

		registry.fill(HIST("eta_2D_MC_Fwd_JoinedTable_Type4"), particle.eta(), track.eta());
	      	registry.fill(HIST("pT_2D_MC_Fwd_JoinedTable_Type4"), particle.pt(), track.pt());
	      	registry.fill(HIST("phi_2D_MC_Fwd_JoinedTable_Type4"), normalize(particle.phi()), track.phi());
	}
        registry.fill(HIST("etaDiff_MC_Fwd_JoinedTable"), particle.eta() - track.eta());
      	registry.fill(HIST("pTDiff_MC_Fwd_JoinedTable"), particle.pt() - track.pt());
      	registry.fill(HIST("phiDiff_MC_Fwd_JoinedTable"), normalize(particle.phi() - track.phi()));

	registry.fill(HIST("eta_2D_MC_Fwd_JoinedTable"), particle.eta(), track.eta());
	registry.fill(HIST("pT_2D_MC_Fwd_JoinedTable"), particle.pt(), track.pt());
	registry.fill(HIST("phi_2D_MC_Fwd_JoinedTable"),normalize(particle.phi()), track.phi());
      }
    }

  }

  PROCESS_SWITCH(IterateMuons, processResolutionFwdTacks, "Process reco/gen matching FwdTracks", true);
	

  void processOffset(aod::Collisions::iterator const& collision, aod::FwdTracks const& muons, aod::Tracks const& tracks)
  {
  }
  PROCESS_SWITCH(IterateMuons, processOffset, "Process Offset FwdTracks", true);
	

  
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
	adaptAnalysisTask<analysemfttracks>(cfgc),
	adaptAnalysisTask<IterateMuons>(cfgc)
  };
}

/*struct analysemfttracks {
  

  void process(o2::aod::Collision const& collision, o2::aod::Tracks const& tracks)
  {
      
  }

  void processMFT(o2::aod::Collision const& collision, o2::aod::MFTTracks const& tracks)
  {
      
  }
  PROCESS_SWITCH(analysemfttracks, processMFT, "Process MFT info", true);


};
*/
