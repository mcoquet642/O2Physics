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

/// \file taskMFTTrkEfficiency.cxx
///

#include "DataFormatsITSMFT/TimeDeadMap.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include <CCDB/BasicCCDBManager.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MFTTrackLabeled = soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>;

/// default histogram output binning
namespace mft_trk_eff_bins
{
static constexpr int nBinsPt = 24;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double binsPt[nBinsPt + 1] = {
  0., 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0,
  2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0,
  15.0, 18.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0};
auto vecBinsPt = std::vector<double>{binsPt, binsPt + nBinsPt + 1};

constexpr float mX[936] = {
  -8.8, -8.8, 8.2, 8.2, 9.9, 9.9, -9.9, -9.9, -8.2, -8.2,
  8.8, 8.8, -7.1, -7.1, -7.1, -5.4, -5.4, -5.4, -3.7, -3.7,
  -3.7, -2, -2, -2, -0.3, -0.3, -0.3, 1.4, 1.4, 1.4,
  3.1, 3.1, 3.1, 4.8, 4.8, 4.8, 6.5, 6.5, 6.5, -6.5,
  -6.5, -6.5, -4.8, -4.8, -4.8, -3.1, -3.1, -3.1, -1.4, -1.4,
  -1.4, 0.3, 0.3, 0.3, 2, 2, 2, 3.7, 3.7, 3.7,
  5.4, 5.4, 5.4, 7.1, 7.1, 7.1, -8.8, -8.8, 8.2, 8.2,
  9.9, 9.9, -9.9, -9.9, -8.2, -8.2, 8.8, 8.8, -7.1, -7.1,
  -7.1, -5.4, -5.4, -5.4, -3.7, -3.7, -3.7, -2, -2, -2,
  -0.3, -0.3, -0.3, 1.4, 1.4, 1.4, 3.1, 3.1, 3.1, 4.8,
  4.8, 4.8, 6.5, 6.5, 6.5, -6.5, -6.5, -6.5, -4.8, -4.8,
  -4.8, -3.1, -3.1, -3.1, -1.4, -1.4, -1.4, 0.3, 0.3, 0.3,
  2, 2, 2, 3.7, 3.7, 3.7, 5.4, 5.4, 5.4, 7.1,
  7.1, 7.1, -10.5, -10.5, 9.9, 9.9, -9.9, -9.9, 10.5, 10.5,
  -8.8, -8.8, -8.8, -7.1, -7.1, -7.1, -2, -2, -2, -0.3,
  -0.3, -0.3, 1.4, 1.4, 1.4, 6.5, 6.5, 6.5, 8.2, 8.2,
  8.2, -8.2, -8.2, -8.2, -6.5, -6.5, -6.5, -1.4, -1.4, -1.4,
  0.3, 0.3, 0.3, 2, 2, 2, 7.1, 7.1, 7.1, 8.8,
  8.8, 8.8, -5.4, -5.4, -5.4, -5.4, -3.7, -3.7, -3.7, -3.7,
  3.1, 3.1, 3.1, 3.1, 4.8, 4.8, 4.8, 4.8, -4.8, -4.8,
  -4.8, -4.8, -3.1, -3.1, -3.1, -3.1, 3.7, 3.7, 3.7, 3.7,
  5.4, 5.4, 5.4, 5.4, -12.2, -12.2, -12.2, -10.5, -10.5, -10.5,
  9.9, 9.9, 9.9, 11.6, 11.6, 11.6, 13.3, 13.3, 13.3, -13.3,
  -13.3, -13.3, -11.6, -11.6, -11.6, -9.9, -9.9, -9.9, 10.5, 10.5,
  10.5, 12.2, 12.2, 12.2, -8.8, -8.8, -8.8, -8.8, -7.1, -7.1,
  -7.1, -7.1, -5.4, -5.4, -5.4, -5.4, -3.7, -3.7, -3.7, -3.7,
  -2, -2, -2, -2, -0.3, -0.3, -0.3, -0.3, 1.4, 1.4,
  1.4, 1.4, 3.1, 3.1, 3.1, 3.1, 4.8, 4.8, 4.8, 4.8,
  6.5, 6.5, 6.5, 6.5, 8.2, 8.2, 8.2, 8.2, -8.2, -8.2,
  -8.2, -8.2, -6.5, -6.5, -6.5, -6.5, -4.8, -4.8, -4.8, -4.8,
  -3.1, -3.1, -3.1, -3.1, -1.4, -1.4, -1.4, -1.4, 0.3, 0.3,
  0.3, 0.3, 2, 2, 2, 2, 3.7, 3.7, 3.7, 3.7,
  5.4, 5.4, 5.4, 5.4, 7.1, 7.1, 7.1, 7.1, 8.8, 8.8,
  8.8, 8.8, -13.9, -13.9, -13.9, -12.2, -12.2, -12.2, 11.6, 11.6,
  11.6, 13.3, 13.3, 13.3, -13.3, -13.3, -13.3, -11.6, -11.6, -11.6,
  12.2, 12.2, 12.2, 13.9, 13.9, 13.9, -10.5, -10.5, -10.5, -10.5,
  -8.8, -8.8, -8.8, -8.8, -3.7, -3.7, -3.7, -3.7, -2, -2,
  -2, -2, -0.3, -0.3, -0.3, -0.3, 1.4, 1.4, 1.4, 1.4,
  3.1, 3.1, 3.1, 3.1, 8.2, 8.2, 8.2, 8.2, 9.9, 9.9,
  9.9, 9.9, -9.9, -9.9, -9.9, -9.9, -8.2, -8.2, -8.2, -8.2,
  -3.1, -3.1, -3.1, -3.1, -1.4, -1.4, -1.4, -1.4, 0.3, 0.3,
  0.3, 0.3, 2, 2, 2, 2, 3.7, 3.7, 3.7, 3.7,
  8.8, 8.8, 8.8, 8.8, 10.5, 10.5, 10.5, 10.5, -7.1, -7.1,
  -7.1, -7.1, -7.1, -5.4, -5.4, -5.4, -5.4, -5.4, 4.8, 4.8,
  4.8, 4.8, 4.8, 6.5, 6.5, 6.5, 6.5, 6.5, -6.5, -6.5,
  -6.5, -6.5, -6.5, -4.8, -4.8, -4.8, -4.8, -4.8, 5.4, 5.4,
  5.4, 5.4, 5.4, 7.1, 7.1, 7.1, 7.1, 7.1, 8.8, 8.8,
  -8.2, -8.2, -9.9, -9.9, 9.9, 9.9, 8.2, 8.2, -8.8, -8.8,
  7.1, 7.1, 7.1, 5.4, 5.4, 5.4, 3.7, 3.7, 3.7, 2,
  2, 2, 0.3, 0.3, 0.3, -1.4, -1.4, -1.4, -3.1, -3.1,
  -3.1, -4.8, -4.8, -4.8, -6.5, -6.5, -6.5, 6.5, 6.5, 6.5,
  4.8, 4.8, 4.8, 3.1, 3.1, 3.1, 1.4, 1.4, 1.4, -0.3,
  -0.3, -0.3, -2, -2, -2, -3.7, -3.7, -3.7, -5.4, -5.4,
  -5.4, -7.1, -7.1, -7.1, 8.8, 8.8, -8.2, -8.2, -9.9, -9.9,
  9.9, 9.9, 8.2, 8.2, -8.8, -8.8, 7.1, 7.1, 7.1, 5.4,
  5.4, 5.4, 3.7, 3.7, 3.7, 2, 2, 2, 0.3, 0.3,
  0.3, -1.4, -1.4, -1.4, -3.1, -3.1, -3.1, -4.8, -4.8, -4.8,
  -6.5, -6.5, -6.5, 6.5, 6.5, 6.5, 4.8, 4.8, 4.8, 3.1,
  3.1, 3.1, 1.4, 1.4, 1.4, -0.3, -0.3, -0.3, -2, -2,
  -2, -3.7, -3.7, -3.7, -5.4, -5.4, -5.4, -7.1, -7.1, -7.1,
  10.5, 10.5, -9.9, -9.9, 9.9, 9.9, -10.5, -10.5, 8.8, 8.8,
  8.8, 7.1, 7.1, 7.1, 2, 2, 2, 0.3, 0.3, 0.3,
  -1.4, -1.4, -1.4, -6.5, -6.5, -6.5, -8.2, -8.2, -8.2, 8.2,
  8.2, 8.2, 6.5, 6.5, 6.5, 1.4, 1.4, 1.4, -0.3, -0.3,
  -0.3, -2, -2, -2, -7.1, -7.1, -7.1, -8.8, -8.8, -8.8,
  5.4, 5.4, 5.4, 5.4, 3.7, 3.7, 3.7, 3.7, -3.1, -3.1,
  -3.1, -3.1, -4.8, -4.8, -4.8, -4.8, 4.8, 4.8, 4.8, 4.8,
  3.1, 3.1, 3.1, 3.1, -3.7, -3.7, -3.7, -3.7, -5.4, -5.4,
  -5.4, -5.4, 12.2, 12.2, 12.2, 10.5, 10.5, 10.5, -9.9, -9.9,
  -9.9, -11.6, -11.6, -11.6, -13.3, -13.3, -13.3, 13.3, 13.3, 13.3,
  11.6, 11.6, 11.6, 9.9, 9.9, 9.9, -10.5, -10.5, -10.5, -12.2,
  -12.2, -12.2, 8.8, 8.8, 8.8, 8.8, 7.1, 7.1, 7.1, 7.1,
  5.4, 5.4, 5.4, 5.4, 3.7, 3.7, 3.7, 3.7, 2, 2,
  2, 2, 0.3, 0.3, 0.3, 0.3, -1.4, -1.4, -1.4, -1.4,
  -3.1, -3.1, -3.1, -3.1, -4.8, -4.8, -4.8, -4.8, -6.5, -6.5,
  -6.5, -6.5, -8.2, -8.2, -8.2, -8.2, 8.2, 8.2, 8.2, 8.2,
  6.5, 6.5, 6.5, 6.5, 4.8, 4.8, 4.8, 4.8, 3.1, 3.1,
  3.1, 3.1, 1.4, 1.4, 1.4, 1.4, -0.3, -0.3, -0.3, -0.3,
  -2, -2, -2, -2, -3.7, -3.7, -3.7, -3.7, -5.4, -5.4,
  -5.4, -5.4, -7.1, -7.1, -7.1, -7.1, -8.8, -8.8, -8.8, -8.8,
  13.9, 13.9, 13.9, 12.2, 12.2, 12.2, -11.6, -11.6, -11.6, -13.3,
  -13.3, -13.3, 13.3, 13.3, 13.3, 11.6, 11.6, 11.6, -12.2, -12.2,
  -12.2, -13.9, -13.9, -13.9, 10.5, 10.5, 10.5, 10.5, 8.8, 8.8,
  8.8, 8.8, 3.7, 3.7, 3.7, 3.7, 2, 2, 2, 2,
  0.3, 0.3, 0.3, 0.3, -1.4, -1.4, -1.4, -1.4, -3.1, -3.1,
  -3.1, -3.1, -8.2, -8.2, -8.2, -8.2, -9.9, -9.9, -9.9, -9.9,
  9.9, 9.9, 9.9, 9.9, 8.2, 8.2, 8.2, 8.2, 3.1, 3.1,
  3.1, 3.1, 1.4, 1.4, 1.4, 1.4, -0.3, -0.3, -0.3, -0.3,
  -2, -2, -2, -2, -3.7, -3.7, -3.7, -3.7, -8.8, -8.8,
  -8.8, -8.8, -10.5, -10.5, -10.5, -10.5, 7.1, 7.1, 7.1, 7.1,
  7.1, 5.4, 5.4, 5.4, 5.4, 5.4, -4.8, -4.8, -4.8, -4.8,
  -4.8, -6.5, -6.5, -6.5, -6.5, -6.5, 6.5, 6.5, 6.5, 6.5,
  6.5, 4.8, 4.8, 4.8, 4.8, 4.8, -5.4, -5.4, -5.4, -5.4,
  -5.4, -7.1, -7.1, -7.1, -7.1, -7.1};

constexpr float mY[936] = {
  -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715,
  -1.7, -4.715, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -3.546, -6.561, -9.576, -3.8, -6.815, -9.83, -3.706, -6.721, -9.736,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -3.706, -6.721,
  -9.736, -3.8, -6.815, -9.83, -3.546, -6.561, -9.576, -1.7, -4.715, -7.73,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -1.7, -4.715,
  -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -3.546, -6.561, -9.576,
  -3.8, -6.815, -9.83, -3.706, -6.721, -9.736, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -3.706, -6.721, -9.736, -3.8, -6.815, -9.83,
  -3.546, -6.561, -9.576, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -3.55, -6.565, -9.58, -3.8,
  -6.815, -9.83, -3.71, -6.725, -9.74, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -3.71, -6.725, -9.74,
  -3.8, -6.815, -9.83, -3.55, -6.565, -9.58, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -3.5, -6.515, -9.53, -12.545,
  -4.735, -7.75, -10.765, -13.78, -4.9, -7.915, -10.93, -13.945, -4.84, -7.855,
  -10.87, -13.885, -3.97, -6.985, -10, -13.015, -1.7, -4.715, -7.73, -10.745,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745,
  -3.97, -6.985, -10, -13.015, -4.84, -7.855, -10.87, -13.885, -4.9, -7.915,
  -10.93, -13.945, -4.735, -7.75, -10.765, -13.78, -3.5, -6.515, -9.53, -12.545,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -10.745,
  -1.7, -4.715, -7.73, -10.745, -4.125, -7.14, -10.155, -13.17, -5.155, -8.17,
  -11.185, -14.2, -5.3, -8.315, -11.33, -14.345, -5.245, -8.26, -11.275, -14.29,
  -4.5, -7.515, -10.53, -13.545, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745,
  -4.5, -7.515, -10.53, -13.545, -5.245, -8.26, -11.275, -14.29, -5.3, -8.315,
  -11.33, -14.345, -5.155, -8.17, -11.185, -14.2, -4.125, -7.14, -10.155, -13.17,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -13.76, -1.7, -4.715, -7.73, -10.745, -13.76, -1.7, -4.715,
  -7.73, -10.745, -13.76, -1.7, -4.715, -7.73, -10.745, -13.76, -1.7, -4.715,
  -7.73, -10.745, -13.76, -1.7, -4.715, -7.73, -10.745, -13.76, -1.7, -4.715,
  -7.73, -10.745, -13.76, -1.7, -4.715, -7.73, -10.745, -13.76, 1.7, 4.715,
  1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 3.546,
  6.561, 9.576, 3.8, 6.815, 9.83, 3.706, 6.721, 9.736, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 3.706, 6.721, 9.736, 3.8,
  6.815, 9.83, 3.546, 6.561, 9.576, 1.7, 4.715, 7.73, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715,
  1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 1.7, 4.715, 7.73, 3.546, 6.561, 9.576, 3.8, 6.815,
  9.83, 3.706, 6.721, 9.736, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 3.706, 6.721, 9.736, 3.8, 6.815, 9.83, 3.546, 6.561,
  9.576, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 3.55, 6.565, 9.58, 3.8, 6.815, 9.83,
  3.71, 6.725, 9.74, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 1.7, 4.715, 7.73, 3.71, 6.725, 9.74, 3.8, 6.815,
  9.83, 3.55, 6.565, 9.58, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 10.745, 3.5, 6.515, 9.53, 12.545, 4.735, 7.75,
  10.765, 13.78, 4.9, 7.915, 10.93, 13.945, 4.84, 7.855, 10.87, 13.885,
  3.97, 6.985, 10, 13.015, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 3.97, 6.985,
  10, 13.015, 4.84, 7.855, 10.87, 13.885, 4.9, 7.915, 10.93, 13.945,
  4.735, 7.75, 10.765, 13.78, 3.5, 6.515, 9.53, 12.545, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715,
  7.73, 10.745, 4.125, 7.14, 10.155, 13.17, 5.155, 8.17, 11.185, 14.2,
  5.3, 8.315, 11.33, 14.345, 5.245, 8.26, 11.275, 14.29, 4.5, 7.515,
  10.53, 13.545, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 4.5, 7.515,
  10.53, 13.545, 5.245, 8.26, 11.275, 14.29, 5.3, 8.315, 11.33, 14.345,
  5.155, 8.17, 11.185, 14.2, 4.125, 7.14, 10.155, 13.17, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  13.76, 1.7, 4.715, 7.73, 10.745, 13.76, 1.7, 4.715, 7.73, 10.745,
  13.76, 1.7, 4.715, 7.73, 10.745, 13.76, 1.7, 4.715, 7.73, 10.745,
  13.76, 1.7, 4.715, 7.73, 10.745, 13.76, 1.7, 4.715, 7.73, 10.745,
  13.76, 1.7, 4.715, 7.73, 10.745, 13.76};

} // namespace mft_trk_eff_bins

struct taskMFTTrkEfficiency {

  Configurable<int> nClusters{"nClusters", 5, "Minimum number of clusters per track"}; /// Muon track type to be selected if value >=0 (no selection by default)

  Configurable<double> ptMuonMin{"ptMin", 0.1, "Lower bound of pT"};                 /// Muon minimum pt to be studied
  Configurable<double> chi2ndfMax{"chi2ndfMax", 40, "cut on the chi2/ndf of track"}; /// Muon minimum pt to be studied
  Configurable<double> etaMax{"etaMax", -2.5, "Upper bound of eta"};                 /// Muon minimum |eta| to be studied
  Configurable<double> etaMin{"etaMin", -4.0, "Lower bound of eta"};                 /// Muon maximum |eta| to be studied

  Configurable<std::vector<double>> binsMuonPt{"binsPt", std::vector<double>{mft_trk_eff_bins::vecBinsPt}, "pT bin limits"}; /// Pt intervals for the histograms
  Configurable<int> nEtaBins{"nEtaBins", 400, "Number of Eta bins"};
  Configurable<int> nPhiBins{"nPhiBins", 400, "Number of Phi bins"};
  Configurable<int64_t> fConfigNoLaterThan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  const char* ccdbpath = "MFT/Calib/TimeDeadMap";
  const char* ccdburl = "http://alice-ccdb.cern.ch";

  int nBCsPerOrbit = 3564;
  int mRunNumber;
  unsigned long mOrbit;
  unsigned long mPrevOrbit;

  o2::itsmft::TimeDeadMap* deadmap = nullptr;

  std::vector<unsigned long> orbits = {};
  std::array<std::vector<int>, 10> chipsPerLayer;

  std::array<TH2I*, 10> layerMasks;

  const o2::itsmft::ChipMappingMFT maping;
  const std::array<o2::itsmft::MFTChipMappingData, 936> chipMap = maping.getChipMappingData();

  float dX = 1.7;
  float dY = 3.015;

  std::array<float, 10> layersZ = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};

  ///  Initialize: configure, create specifics
  void init(o2::framework::InitContext&)
  {

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(fConfigNoLaterThan.value);

    auto vbins = (std::vector<double>)binsMuonPt;
    const AxisSpec axisEta{nEtaBins, etaMin, etaMax, "#eta"};
    const AxisSpec axisPt{vbins, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPhi{nPhiBins, -3.14, 3.14, "#varphi"};

    const AxisSpec axisEtaGen{nEtaBins, etaMin, etaMax, "#eta Gen"};
    const AxisSpec axisPtGen{vbins, "#it{p}_{T} (GeV/#it{c}) Gen"};
    const AxisSpec axisPhiGen{nPhiBins, -3.14, 3.14, "#varphi Gen"};

    const AxisSpec axisNhits{15, -0.5, 15.5, ""};
    // Labels for the chambers hit per station, numbering starting from 1
    // i.e. (Nij, N0j, Ni0) correspond to hit on i-j, hit on j not on i, hit on i not on j
    const char* elabels[15] = {"N12", "N10", "N02", "N34", "N30", "N04", "N56", "N50", "N06", "N78", "N70", "N08", "N910", "N90", "N010"};

    HistogramConfigSpec defaultNhitsEtaPtPhi({HistType::kTHnF, {{axisNhits}, {axisEta}, {axisPt}, {axisPhi}}});

    registry.add("hPtRecPtGen", "hPtRecPtGen", {HistType::kTH2F, {{axisPt}, {axisPtGen}}});
    registry.add("hEtaRecEtaGen", "hEtaRecEtaGen", {HistType::kTH2F, {{axisEta}, {axisEtaGen}}});
    registry.add("hPhiRecPhiGen", "hPhiRecPhiGen", {HistType::kTH2F, {{axisPhi}, {axisPhiGen}}});

    registry.add("hHitsEtaPtPhi", "hHitsEtaPtPhi", defaultNhitsEtaPtPhi, false);
    auto hHitsEtaPtPhi = registry.get<THn>(HIST("hHitsEtaPtPhi"));
    for (int i = 0; i < 15; i++)
      hHitsEtaPtPhi->GetAxis(0)->SetBinLabel(i + 1, elabels[i]);

    for (int i = 0; i < 10; i++) {
      layerMasks[i] = new TH2I("", "", nEtaBins, -4.0, 2.0, nPhiBins, -TMath::Pi(), TMath::Pi());
    }

    mRunNumber = 0;
    mOrbit = 0;
    mPrevOrbit = 0;

  } //! end of Initialize: configure, create specifics

  void decodeChipVector(std::vector<uint16_t>& inChips, std::array<std::vector<int>, 10>& outChipsPerLayer)
  {
    bool prevIsDead = false;
    uint16_t lastDead = -1;
    for (auto& chip : inChips) {
      if (chip & (uint16_t)(0x8000)) {
        auto first_chip = chip - 0x8000;
        outChipsPerLayer[chipMap[first_chip].layer].push_back(first_chip);
        prevIsDead = true;
        lastDead = first_chip;
      } else {
        if (prevIsDead) {
          for (int i = 1; i < chip - lastDead; i++) {
            outChipsPerLayer[chipMap[i + lastDead].layer].push_back(i + lastDead);
          }
        }
        outChipsPerLayer[chipMap[chip].layer].push_back(chip);
        prevIsDead = false;
      }
    }
  }

  double normPhi(double phi)
  {
    while (phi <= -TMath::Pi())
      phi += 2.0 * TMath::Pi();
    while (phi > TMath::Pi())
      phi -= 2.0 * TMath::Pi();
    return phi;
  }

  std::tuple<double, double, double, double> computeEtaPhiCoverage(double posX, double posY, double dX, double dY, double z)
  {
    double x1 = posX - dX / 2.0;
    double x2 = posX + dX / 2.0;
    double y1 = posY - dY / 2.0;
    double y2 = posY + dY / 2.0;

    std::vector<double> etas;
    std::vector<double> phis;
    etas.reserve(4);
    phis.reserve(4);

    for (double x : {x1, x2}) {
      for (double y : {y1, y2}) {
        double r = std::sqrt(x * x + y * y);
        double theta = std::atan2(r, z);
        double eta = -std::log(std::tan(theta / 2.0));
        double phi = std::atan2(y, x);
        etas.push_back(eta);
        phis.push_back(phi);
      }
    }

    double eta_min = *std::min_element(etas.begin(), etas.end());
    double eta_max = *std::max_element(etas.begin(), etas.end());
    double phiMin = *std::min_element(phis.begin(), phis.end());
    double phiMax = *std::max_element(phis.begin(), phis.end());

    if (phiMax - phiMin > M_PI) {
      for (auto& phi : phis) {
        if (phi < 0)
          phi += 2 * M_PI;
      }
      phiMin = *std::min_element(phis.begin(), phis.end());
      phiMax = *std::max_element(phis.begin(), phis.end());
    }

    return std::make_tuple(eta_min, eta_max, phiMin, phiMax);
  }

  void applyChipToMask(TH2* mask, double etaMin, double etaMax, double phiMin, double phiMax)
  {
    if (!mask)
      return;

    int etaBinMin = mask->GetXaxis()->FindBin(etaMin);
    int etaBinMax = mask->GetXaxis()->FindBin(etaMax);
    int phiBinMin = mask->GetYaxis()->FindBin(phiMin);
    int phiBinMax = mask->GetYaxis()->FindBin(phiMax);

    for (int iEta = etaBinMin; iEta <= etaBinMax; ++iEta) {
      for (int iPhi = phiBinMin; iPhi <= phiBinMax; ++iPhi) {
        mask->SetBinContent(iEta, iPhi, 1.0);
      }
    }
  }

  bool isSelected(double eta, double pt, int nclusters, double chi2ndf)
  {
    return ((pt > ptMuonMin) && (eta > etaMin) && (eta < etaMax) && (nclusters >= nClusters) && (chi2ndf < chi2ndfMax));
  }

  void computeExclusionMap(int layer)
  {
    auto chips = chipsPerLayer[layer];
    float z = layersZ[layer];
    for (auto& chip : chips) {
      float posX = mft_trk_eff_bins::mX[chip];
      float posY = mft_trk_eff_bins::mY[chip];
      auto [eta_min, eta_max, phi_min, phi_max] = computeEtaPhiCoverage(posX, posY, dX, dY, z);
      applyChipToMask(layerMasks[layer], eta_min, eta_max, phi_min, phi_max);
    }
  }

  bool isExcluded(double eta, double phi, int layer)
  {
    int binEta = layerMasks[layer]->GetXaxis()->FindBin(eta);
    int binPhi = layerMasks[layer]->GetYaxis()->FindBin(phi);
    bool isBad = (layerMasks[layer]->GetBinContent(binEta, binPhi) > 0);
    return isBad;
  }

  void FillHistosWeight(double eta, double pt, double phi, uint64_t map, double /*etaGen*/, double /*ptGen*/, double /*phiGen*/)
  {
    double weighteta = 1, weightpt = 1, weightphi = 1; // default weight set to unity for now: no effect
    double etaw = eta * weighteta;
    double ptw = pt * weightpt;
    double phiw = phi * weightphi;
    FillHistos(etaw, ptw, phiw, map);
  }

  /// Filling histograms from generated & reconstructed information
  void FillHistosMC(double eta, double pt, double phi, uint64_t map, double etaGen, double ptGen, double phiGen)
  {
    registry.fill(HIST("hPtRecPtGen"), pt, ptGen);
    registry.fill(HIST("hEtaRecEtaGen"), eta, etaGen);
    registry.fill(HIST("hPhiRecPhiGen"), phi, phiGen);
    FillHistosWeight(eta, pt, phi, map, etaGen, ptGen, phiGen);
  }

  void FillHistos(double eta, double pt, double phi, uint64_t map)
  {
    bool iN[10];
    for (int ilayer = 0; ilayer < 10; ilayer++) {
      iN[ilayer] = (map >> (ilayer * 6)) & 0x3F;
    }
    for (int i = 0; i < 5; i++) {
      if (iN[2 * i] && iN[2 * i + 1]) {
        registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(3 * i, eta, pt, phi);
      }
      if (iN[2 * i] && (!iN[2 * i + 1])) {
        registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(3 * i + 1, eta, pt, phi);
      }
      if ((!iN[2 * i]) && iN[2 * i + 1]) {
        registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(3 * i + 2, eta, pt, phi);
      }
    }
  }

  /// Filling histograms from reconstructed quantities
  void FillHistosExclude(double eta, double pt, double phi, uint64_t map)
  {
    bool iN[10];
    for (int ilayer = 0; ilayer < 10; ilayer++) {
      iN[ilayer] = (map >> (ilayer * 6)) & 0x3F;
    }
    for (int i = 0; i < 5; i++) {
      if (iN[2 * i] && iN[2 * i + 1]) {
        registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(3 * i, eta, pt, phi);
      }
      if (iN[2 * i] && (!iN[2 * i + 1]) && !isExcluded(eta, phi, 2 * i + 1)) {
        registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(3 * i + 1, eta, pt, phi);
      }
      if ((!iN[2 * i]) && iN[2 * i + 1] && !isExcluded(eta, phi, 2 * i)) {
        registry.get<THn>(HIST("hHitsEtaPtPhi"))->Fill(3 * i + 2, eta, pt, phi);
      }
    }
  }

  void processReco(aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&)
  {
    for (auto& mfttrack : mfttracks) {
      auto eta = mfttrack.eta();
      auto pt = mfttrack.pt();
      auto phi = mfttrack.phi();
      auto nclusters = mfttrack.nClusters();
      auto chi2 = mfttrack.chi2();
      auto chi2ndf = chi2 / (2. * nclusters - 5);
      auto map = mfttrack.mftClusterSizesAndTrackFlags();
      if (isSelected(eta, pt, nclusters, chi2ndf)) {
        FillHistos(eta, pt, phi, map);
      }
    }
  }
  PROCESS_SWITCH(taskMFTTrkEfficiency, processReco, "process reconstructed information", true);

  void processRecoEffOnly(aod::MFTTracks const& mfttracks, aod::BCsWithTimestamps const&, aod::Collisions const&)
  {
    for (auto& mfttrack : mfttracks) {
      if (mfttrack.has_collision()) {
        auto bc = mfttrack.collision_as<aod::Collisions>().bc_as<aod::BCsWithTimestamps>();
        if (mRunNumber != bc.runNumber()) {
          deadmap = ccdb->getForTimeStamp<o2::itsmft::TimeDeadMap>(ccdbpath, bc.timestamp());
          if (deadmap != nullptr) {
            LOGF(info, "Using deadmap for run %d", bc.runNumber());
          } else {
            LOGF(fatal, "DeadMap is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
          }
          mRunNumber = bc.runNumber();
          orbits = deadmap->getEvolvingMapKeys();
        }
        if (mOrbit != (bc.globalBC() / nBCsPerOrbit)) {
          mOrbit = (bc.globalBC() / nBCsPerOrbit);
          std::vector<uint16_t> encodeChips;
          unsigned long lower_orbit = deadmap->getMapAtOrbit(mOrbit, encodeChips);
          if ((mOrbit - lower_orbit) > mPrevOrbit) {
            for (auto& v : chipsPerLayer) {
              v.clear();
            }
            for (auto& h : layerMasks) {
              if (h)
                h->Reset("ICES");
            }
            decodeChipVector(encodeChips, chipsPerLayer);
            for (int i = 0; i < 10; i++) {
              computeExclusionMap(i);
            }
            mPrevOrbit = mOrbit - lower_orbit;
          }
        }
      }
      auto eta = mfttrack.eta();
      auto pt = mfttrack.pt();
      auto phi = mfttrack.phi();
      auto nclusters = mfttrack.nClusters();
      auto chi2 = mfttrack.chi2();
      auto chi2ndf = chi2 / (2. * nclusters - 5);
      auto map = mfttrack.mftClusterSizesAndTrackFlags();
      if (isSelected(eta, pt, nclusters, chi2ndf)) {
        FillHistosExclude(eta, pt, phi, map);
      }
    }
  }
  PROCESS_SWITCH(taskMFTTrkEfficiency, processRecoEffOnly, "process reconstructed information", false);

  void processSim(MFTTrackLabeled const& mfttracks, aod::McParticles const& /*mcTracks*/)
  {
    for (auto& mfttrack : mfttracks) {
      auto eta = mfttrack.eta();
      auto pt = mfttrack.pt();
      auto phi = mfttrack.phi();
      auto nclusters = mfttrack.nClusters();
      auto chi2 = mfttrack.chi2();
      auto chi2ndf = chi2 / (2. * nclusters - 5);
      auto map = mfttrack.mftClusterSizesAndTrackFlags();
      if (isSelected(eta, pt, nclusters, chi2ndf)) {
        if (!mfttrack.has_mcParticle()) {
          continue;
        }
        auto mctrack = mfttrack.template mcParticle_as<aod::McParticles>();
        FillHistosMC(eta, pt, phi, map, mctrack.eta(), mctrack.pt(), mctrack.phi());
      }
    }
  }
  PROCESS_SWITCH(taskMFTTrkEfficiency, processSim, "process simulated information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<taskMFTTrkEfficiency>(cfgc)};
}
