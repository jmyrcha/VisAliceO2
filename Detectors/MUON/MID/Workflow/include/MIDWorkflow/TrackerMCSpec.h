// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   MIDWorkflow/TrackerSpec.h
/// \brief  Data processor spec for MID MC tracker device
/// \author Diego Stocco <Diego.Stocco at cern.ch>
/// \date   27 September 2019

#ifndef O2_MID_TRACKERMCSPEC_H
#define O2_MID_TRACKERMCSPEC_H

#include "Framework/DataProcessorSpec.h"

namespace o2
{
namespace mid
{
framework::DataProcessorSpec getTrackerMCSpec();
}
} // namespace o2

#endif //O2_MID_TRACKERMCSPEC_H
