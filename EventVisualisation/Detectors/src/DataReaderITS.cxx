// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file   DataReaderITS.cxx
/// \brief  ITS Detector-specific reading from file(s)
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#include "EventVisualisationDetectors/DataReaderITS.h"
#include <TTree.h>
#include <TVector2.h>
#include <TError.h>
#include "DataFormatsITSMFT/ROFRecord.h"

namespace o2  {
namespace event_visualisation {

    DataReaderITS::DataReaderITS() = default;

    void DataReaderITS::open() {
        TString clusterFile = "o2clus_its.root";
        TString trackFile = "o2trac_its.root";

        this->tracFile = TFile::Open(trackFile);
        this->clusFile = TFile::Open(clusterFile);

        TTree* roft = (TTree*)this->tracFile->Get("ITSTracksROF");
        TTree* rofc = (TTree*)this->clusFile->Get("ITSClustersROF");

        if(roft != nullptr && rofc != nullptr) {
            //TTree *tracks = (TTree*)this->tracFile->Get("o2sim");
            // temporary number of readout frames as number of events
            TTree *tracksRof = (TTree*)this->tracFile->Get("ITSTracksROF");

            //TTree *clusters = (TTree*)this->clusFile->Get("o2sim");
            TTree *clustersRof = (TTree*)this->clusFile->Get("ITSClustersROF");

            //Read all track RO frames to a buffer to count number of elements
            std::vector<o2::itsmft::ROFRecord> *trackROFrames = nullptr;
            tracksRof->SetBranchAddress("ITSTracksROF", &trackROFrames);
            tracksRof->GetEntry(0);

            //Read all cluster RO frames to a buffer
            std::vector<o2::itsmft::ROFRecord> *clusterROFrames = nullptr;
            clustersRof->SetBranchAddress("ITSClustersROF", &clusterROFrames);
            clustersRof->GetEntry(0);


            if(trackROFrames->size() == clusterROFrames->size()) {
                fMaxEv = trackROFrames->size();
            } else {
                Error("DataReaderITS", "Inconsistent number of readout frames in files");
                exit(1);
            }
        }
    }

    Int_t DataReaderITS::GetEventCount() {
        return fMaxEv;
    }

    TObject *DataReaderITS::getEventData(int no) {
        /// FIXME: Redesign the data reader class
        TList *list = new TList();
        list->Add(this->tracFile);
        list->Add(this->clusFile);
        TVector2 *v = new TVector2(no, 0);
        list->Add(v);
        return list;
    }
}
 }
