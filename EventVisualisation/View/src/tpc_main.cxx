//
// Created by jmy on 10.07.19.
//

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
/// \file    main.cxx
/// \author  Jeremi Niedziela
///

#include "EventVisualisationView/Initializer.h"
#include "EventVisualisationBase/ConfigurationManager.h"
#include "EventVisualisationBase/VisualisationConstants.h"

#include <TApplication.h>
#include <TEveBrowser.h>
#include <TEveManager.h>

#include <ctime>
#include <iostream>

using namespace std;
using namespace o2::event_visualisation;

#include <iostream>
#include <array>
#include <algorithm>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TEveManager.h>
#include <TEveBrowser.h>
#include <TGButton.h>
#include <TGNumberEntry.h>
#include <TGFrame.h>
#include <TGTab.h>
#include <TGLCameraOverlay.h>
#include <TEveFrameBox.h>
#include <TEveQuadSet.h>
#include <TEveTrans.h>
#include <TEvePointSet.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <Rtypes.h>
#include <gsl/span>
#include <DetectorsBase/GeometryManager.h>

#include "EventVisualisationView/MultiView.h"

#include "TPCSimulation/Point.h" // for o2::tpc::HitGroup
#include "TPCBase/Digit.h"
#include "TPCBase/Mapper.h"
#include "TPCReconstruction/RawReader.h"
//#include "DataFormatsTPC/ClusterNative.h"
//#include "DataFormatsTPC/ClusterNativeHelper.h"
#include "DataFormatsTPC/TrackTPC.h"
//#include "DetectorsCommonDataFormats/DetID.h"
//#include "CommonDataFormat/InteractionRecord.h"
#include "TGenericClassInfo.h"

using namespace o2::tpc;

extern TEveManager* gEve;

static TGNumberEntry* gEntry;

class Data
{
public:
    void loadData(int entry);
    void displayData(int entry);
    int getLastEvent() const { return mLastEvent; }
    void setHitsTree(TTree* t);
    void setDigiTree(TTree* t);
    //void setClusTree(TTree* t);
    void setTracTree(TTree* t);

    void setRawReader(std::string infile)
    {
      auto reader = new RawReader();
      if(!reader->addInputFile(-1, -1, -1, infile, -1)) {
        std::cout << "Could not add file for raw reader!" << std::endl;
        return;
      }
      mRawReader = reader;
      mRawReader->loadEvent(0);
    }
    bool getRawData() { return mRawData; }
    void setRawData(bool rawData) { mRawData = rawData; }

private:
    // Data loading members
    //ClusterNativeHelper::Reader reader;
    bool mRawData = false;
    RawReader* mRawReader = nullptr;
    int mLastEvent = 0;
    int kNumberOfSectors = 36;
    std::vector<std::vector<HitGroup>*> mHitsBuffer;
    std::vector<gsl::span<HitGroup>> mHits;
    std::vector<std::vector<Digit>*> mDigitBuffer;
    std::vector<gsl::span<Digit>> mDigits;
    //ClusterNativeAccess clusterIndex;
    //std::unique_ptr<ClusterNative[]> mClusterBuffer;
    //std::vector<ClusterNative>* mClusterBuffer = nullptr;
    //gsl::span<ClusterNative> mClusters;
    //std::unique_ptr<ClusterNativeAccess> mClusters;
    std::vector<TrackTPC>* mTrackBuffer = nullptr;
    gsl::span<TrackTPC> mTracks;
    void loadHits();
    void loadHits(int entry);
    void loadDigits();
    void loadDigits(int entry);
    void loadRawDigits();
    void loadRawDigits(int entry);
    //void loadClusters(int entry);
    void loadTracks(int entry);

    TTree* mHitsTree = nullptr;
    TTree* mDigiTree = nullptr;
    //TTree* mClusTree = nullptr;
    TTree* mTracTree = nullptr;

    // TEve-related members
    TEveElementList* mEvent = nullptr;
    TEveElement* getEveHits();
    TEveElement* getEveDigits();
    //TEveElement* getEveClusters();
    TEveElement* getEveTracks();
} evdata;

void Data::setHitsTree(TTree* tree)
{
  if (tree == nullptr) {
    std::cerr << "No tree for hits!\n";
    return;
  }
  mHitsBuffer.resize(kNumberOfSectors);
  for(int sector = 0; sector < kNumberOfSectors; sector++) {
    std::stringstream hitsStr;
    hitsStr << "TPCHitsShiftedSector" << sector;
    tree->SetBranchAddress(hitsStr.str().c_str(), &(mHitsBuffer[sector]));
  }

  mHitsTree = tree;
}

// Based on: macro/analyzeHits::analyseTPC()
void Data::loadHits(int entry)
{
  static int lastLoaded = -1;

  if (mHitsTree == nullptr)
    return;

  int eventsCount = mHitsTree->GetBranch("TPCHitsShiftedSector0")->GetEntries();
  std::cout << "Hits: Number of events: " << eventsCount << std::endl;
  if ((entry < 0) || (entry >= eventsCount)) {
    std::cerr << "Hits: Out of event range ! " << entry << '\n';
    return;
  }
  if (entry != lastLoaded) {
    for(int sector = 0; sector < kNumberOfSectors; sector++) {
      mHitsBuffer[sector]->clear();
    }
    mHitsTree->GetEntry(entry);
    lastLoaded = entry;
  }

  int size = 0;
  int first = 0;
  mHits.resize(kNumberOfSectors);
  for(int sector = 0; sector < kNumberOfSectors; sector++) {
    int last = mHitsBuffer[sector]->size();
    mHits[sector] = gsl::make_span(&(*mHitsBuffer[sector])[first], last - first);
    size += mHits[sector].size();
  }

  std::cout << "Number of TPC Hits: " << size << '\n';
}

void Data::setDigiTree(TTree* tree)
{
    if (tree == nullptr) {
        std::cerr << "No tree for digits!\n";
        return;
    }

    mDigitBuffer.resize(kNumberOfSectors);
    for(int sector = 0; sector < kNumberOfSectors; sector++) {
      std::stringstream digiStr;
      digiStr << "TPCDigit_" << sector;
      tree->SetBranchAddress(digiStr.str().c_str(), &(mDigitBuffer[sector]));
    }

    mDigiTree = tree;
}

// Based on: CalibRawBase::processEventRawReader()
void Data::loadRawDigits()
{
  mDigitBuffer.resize(1);

  int size = 0;
  o2::tpc::PadPos padPos;
  while (std::shared_ptr<std::vector<uint16_t>> data = mRawReader->getNextData(padPos)) {
    if (!data)
      continue;

    CRU cru(mRawReader->getRegion());

    // row is local in region (CRU)
    const int row = padPos.getRow();
    const int pad = padPos.getPad();
    if (row == 255 || pad == 255)
      continue;

    int timeBin = 0;
    for (const auto& signalI : *data) {
      const float signal = float(signalI);
      mDigitBuffer[0]->emplace_back(cru, signal, row, pad, timeBin);
      size++;
      ++timeBin;
    }
    std::cout << "timeBin counts: " << timeBin << std::endl;
  }

  int first = 0;
  int last = mDigitBuffer[0]->size();
  mDigits[0] = gsl::make_span(&(*mDigitBuffer[0])[first], last - first);

  std::cout << "Number of TPC Digits: " << size << '\n';
}


void Data::loadRawDigits(int entry)
{
  if (mRawReader == nullptr)
    return;

  if ((entry < 0) || (entry >= mRawReader->getEventNumber())) {
    std::cerr << "Raw digits: Out of event range ! " << entry << '\n';
    return;
  }

  int eventId = mRawReader->loadEvent(entry);

  mLastEvent++;
  loadRawDigits();
}

void Data::loadDigits(int entry)
{
  static int lastLoaded = -1;

  if (mDigiTree == nullptr)
    return;

  int eventsCount = 0;
  eventsCount = mDigiTree->GetBranch("TPCDigit_0")->GetEntries();
  std::cout << "Digits: Number of events: " << eventsCount << std::endl;
  if ((entry < 0) || (entry >= eventsCount)) {
    std::cerr << "Digits: Out of event range ! " << entry << '\n';
    return;
  }
  if (entry != lastLoaded) {
    for(int i = 0; i < kNumberOfSectors; i++) {
      mDigitBuffer[i]->clear();
    }
    mDigiTree->GetEntry(entry);
    lastLoaded = entry;
  }

  int size = 0;
  int first = 0;
  mDigits.resize(kNumberOfSectors);
  for(int i = 0; i < kNumberOfSectors; i++) {
    int last = mDigitBuffer[i]->size();
    mDigits[i] = gsl::make_span(&(*mDigitBuffer[i])[first], last - first);
    size += mDigits[i].size();
  }

  std::cout << "Number of TPC Digits: " << size << '\n';
}


//void Data::setClusTree(TTree* tree)
//{
//    if (tree == nullptr) {
//        std::cerr << "No tree for clusters !\n";
//        return;
//    }
//    tree->SetBranchAddress("TPCCluster", &mClusterBuffer);
//    mClusTree = tree;
//}

//void Data::loadClusters(int entry)
//{
//    static int lastLoaded = -1;
//
//    if (mClusTree == nullptr)
//        return;
//
//    auto event = entry;
//    if ((event < 0) || (event >= mClusTree->GetEntries())) {
//        std::cerr << "Clusters: Out of event range ! " << event << '\n';
//        return;
//    }
//    if (event != lastLoaded) {
//        mClusterBuffer->clear();
//        mClusTree->GetEntry(event);
//        lastLoaded = event;
//    }
//
//    int first = 0, last = mClusterBuffer->size();
//    mClusters = gsl::make_span(&(*mClusterBuffer)[first], last - first);
//
//    std::cout << "Number of TPC Clusters: " << mClusters.size() << '\n';
//}

void Data::setTracTree(TTree* tree)
{
    if (tree == nullptr) {
        std::cerr << "No tree for tracks !\n";
        return;
    }
    tree->SetBranchAddress("TPCTracks", &mTrackBuffer);
    mTracTree = tree;
}

void Data::loadTracks(int entry)
{
    static int lastLoaded = -1;

    if (mTracTree == nullptr)
        return;

    std::cout << "Tracks: Number of events: " << mTracTree->GetEntries() << std::endl;
    if ((entry < 0) || (entry >= mTracTree->GetEntries())) {
        std::cerr << "Tracks: Out of event range ! " << entry << '\n';
        return;
    }
    if (entry != lastLoaded) {
        mTrackBuffer->clear();
        mTracTree->GetEntry(entry);
        lastLoaded = entry;
    }

    int first = 0, last = mTrackBuffer->size();
    mTracks = gsl::make_span(&(*mTrackBuffer)[first], last - first);

    std::cout << "Number of TPC Tracks: " << mTracks.size() << '\n';
}

void Data::loadData(int entry)
{
    loadHits(entry);
    if(evdata.getRawData()) {
      loadRawDigits(entry);
    } else {
      loadDigits(entry);
    }
    //loadClusters(entry);
    loadTracks(entry);
}

TEveElement* Data::getEveHits()
{
  const auto& mapper = Mapper::instance();
  TEvePointSet* hits = new TEvePointSet("hits");
  hits->SetMarkerColor(kYellow);

  for(int i = 0; i < mHits.size(); i++) {
    for (const auto& hv : mHits[i]) {
      for(int j = 0; j < hv.getSize(); j++) {
        const auto& h = hv.getHit(j);
        hits->SetNextPoint(h.GetX(), h.GetY(), h.GetZ());
      }
    }
  }
  return hits;
}

TEveElement* Data::getEveDigits()
{
    const auto& mapper = Mapper::instance();
    TEvePointSet* digits = new TEvePointSet("digits");
    digits->SetMarkerColor(kBlue);

    for(int i = 0; i < mDigits.size(); i++) {
      for (const auto& d : mDigits[i]) {
        const auto pad = mapper.globalPadNumber(PadPos(d.getRow(),
                                                       d.getPad()));
        const auto& localXYZ = mapper.padCentre(pad);
        const auto globalXYZ = mapper.LocalToGlobal(localXYZ,
                                                    CRU(d.getCRU()).sector());
        // TODO: One needs event time0 to get z-coordinate
        digits->SetNextPoint(globalXYZ.X(), globalXYZ.Y(), 0.0f);
      }
    }
    return digits;
}

//TEveElement* Data::getEveClusters()
//{
//    const auto& mapper = Mapper::instance();
//    TEvePointSet* clusters = new TEvePointSet("clusters");
//    clusters->SetMarkerColor(kBlue);
//    for (const auto& c : mClusters) {
//        const auto pad = mapper.globalPadNumber(PadPos(c.getRow(),
//                                                       c.getPad()));
//        const auto& localXYZ = mapper.padCentre(pad);
//        const auto globalXYZ = mapper.LocalToGlobal(localXYZ,
//                                                    CRU(c.getCRU()).sector());
//        clusters->SetNextPoint(globalXYZ.X(), globalXYZ.Y(), globalXYZ.Z());
//    }
//    return clusters;
//}

TEveElement* Data::getEveTracks()
{
    TEveTrackList* tracks = new TEveTrackList("tracks");
    auto prop = tracks->GetPropagator();
    prop->SetMagField(0.5);
    //prop->SetMaxR(50.);
    for (const auto& rec : mTracks) {
        std::array<float, 3> p;
        rec.getPxPyPzGlo(p);
        TEveRecTrackD t;
        t.fP = { p[0], p[1], p[2] };
        t.fSign = (rec.getSign() < 0) ? -1 : 1;
        TEveTrack* track = new TEveTrack(&t, prop);
        track->SetLineColor(kMagenta);
        tracks->AddElement(track);

//        TEvePointSet* tpoints = new TEvePointSet("tclusters");
//        tpoints->SetMarkerColor(kGreen);
//        int nc = rec.getNClusterReferences();
//        while (nc--) {
//            uint8_t sector, row;
//            uint32_t clusterIndexInRow;
//            rec.getClusterReference(nc, sector, row, clusterIndexInRow);
//            const ClusterNative& cl = rec.getCluster(nc, *clusters, sector, row);
//            const ClusterNative& clLast = rec.getCluster(0, *clusters);
//            const auto& gloC = c.getXYZGloRot(*gman);
//            tpoints->SetNextPoint(gloC.X(), gloC.Y(), gloC.Z());
//        }
//        track->AddElement(tpoints);
    }
    tracks->MakeTracks();

    return tracks;
}

void Data::displayData(int entry)
{
    std::string ename("Event #");
    ename += std::to_string(entry);

    // Event display
    auto hits = getEveHits();
    auto digits = getEveDigits();
    //auto clusters = getEveClusters();
    auto tracks = getEveTracks();
    delete mEvent;
    mEvent = new TEveElementList(ename.c_str());
    mEvent->AddElement(hits);
    mEvent->AddElement(digits);
    //mEvent->AddElement(clusters);
    mEvent->AddElement(tracks);
    auto multi = o2::event_visualisation::MultiView::getInstance();
    multi->registerElement(mEvent);

    gEve->Redraw3D(kFALSE);
}

void load(int entry)
{
    int lastEvent = evdata.getLastEvent();
    if (lastEvent > entry) {
        std::cerr << "\nERROR: Cannot stay or go back over events. Please increase the event number !\n\n";
        gEntry->SetIntNumber(lastEvent - 1);
        return;
    }

    gEntry->SetIntNumber(entry);

    std::cout << "\n*** Event #" << entry << " ***\n";
    evdata.loadData(entry);
    evdata.displayData(entry);
}


void load()
{
    auto event = gEntry->GetNumberEntry()->GetIntNumber();
    load(event);
}

void next()
{
    auto event = gEntry->GetNumberEntry()->GetIntNumber();
    event++;
    load(event);
}

void prev()
{
    auto event = gEntry->GetNumberEntry()->GetIntNumber();
    event--;
    load(event);
}

void setupGeometry()
{
  // read path to geometry files from config file
  TEnv settings;
  ConfigurationManager::getInstance().getConfig(settings);

  // get geometry from Geometry Manager and register in multiview
  auto multiView = MultiView::getInstance();

  for(int iDet=0;iDet<NvisualisationGroups;++iDet){
    EVisualisationGroup det = static_cast<EVisualisationGroup>(iDet);
    std::string detName = gVisualisationGroupName[det];
    if(settings.GetValue((detName+".draw").c_str(), false))
    {
      if(   detName=="TPC" || detName=="MCH" || detName=="MTR"
            || detName=="MID" || detName=="MFT" || detName=="AD0"
            || detName=="FMD"){// don't load MUON+MFT and AD and standard TPC to R-Phi view

        multiView->drawGeometryForDetector(detName, false, true, false);
      }
      else if(detName=="RPH"){// special TPC geom from R-Phi view

        multiView->drawGeometryForDetector(detName, false, false, true, false);
      }
      else{// default
        multiView->drawGeometryForDetector(detName);
      }
    }
  }
}

int main(int argc, char **argv)
{
    // create ROOT application environment
    TApplication *app = new TApplication("o2-tpc-eve", &argc, argv);
    app->Connect("TEveBrowser", "CloseWindow()", "TApplication", app, "Terminate()");

    cout<<"Initializing TEveManager"<<endl;
    if(!(gEve=TEveManager::Create())){
        cout<<"FATAL -- Could not create TEveManager!!"<<endl;
        exit(0);
    }

    int entry = 0;
    std::string rawfile = "o2sim.root"; // should be e.g. GBTx0_Run005 - not a file per se?
    std::string hitsfile = "o2sim.root";
    std::string digifile = "tpcdigits.root";
    evdata.setRawData(false);
    std::string tracfile = "tpctracks.root";
    std::string inputGeom = "O2geometry.root";

    // Geometry
    o2::base::GeometryManager::loadGeometry(inputGeom, "FAIRGeom");
    TEveBrowser* browser = gEve->GetBrowser();

    // Event View
    std::cout << "Going to setup the geometry..." << std::endl;
    setupGeometry();

    // Event navigation
    browser->StartEmbedding(TRootBrowser::kBottom);
    auto frame = new TGMainFrame(gClient->GetRoot(), 1000, 600, kVerticalFrame);

    auto h = new TGHorizontalFrame(frame);
    auto b = new TGTextButton(h, "PrevEvnt", "prev()");
    h->AddFrame(b);
    gEntry = new TGNumberEntry(h, 0, 5, -1, TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 10000);
    gEntry->Connect("ValueSet(Long_t)", 0, 0, "load()");
    h->AddFrame(gEntry);
    b = new TGTextButton(h, "NextEvnt", "next()");
    h->AddFrame(b);

    frame->AddFrame(h);

    frame->MapSubwindows();
    frame->MapWindow();
    browser->StopEmbedding("Navigator");

    TFile* file;

    // Data sources
    file = TFile::Open(hitsfile.data());
    if (file && gFile->IsOpen()) {
      evdata.setHitsTree((TTree*)gFile->Get("o2sim"));
    } else
      std::cerr << "\nERROR: Cannot open file: " << hitsfile << "\n\n";

    if (evdata.getRawData()) {
      // TODO: Another raw file for TPC? If at all?
        std::ifstream* rawfileStream = new std::ifstream(rawfile.data(), std::ifstream::binary);
        if (rawfileStream->good()) {
            delete rawfileStream;
            std::cout << "Running with raw digits...\n";
            evdata.setRawReader(rawfile.data());
        } else
            std::cerr << "\nERROR: Cannot open file: " << rawfile << "\n\n";
    } else {
        file = TFile::Open(digifile.data());
        if (file && gFile->IsOpen()) {
            std::cout << "Running with MC digits...\n";
            evdata.setDigiTree((TTree*)gFile->Get("o2sim"));
        } else
            std::cerr << "\nERROR: Cannot open file: " << digifile << "\n\n";
    }

    file = TFile::Open(tracfile.data());
    if (file && gFile->IsOpen()) {
      evdata.setTracTree((TTree*)gFile->Get("tpcrec"));
    }
    else
        std::cerr << "\nERROR: Cannot open file: " << tracfile << "\n\n";

    // Manually adding an event
    gEve->AddEvent(new TEveEventManager("Event", "ALICE TPC Event"));

    load(entry);

    // Start the application
    app->Run(kTRUE);

    // Terminate application
    TEveManager::Terminate();
    app->Terminate();

    return 0;
}

void tpc_main()
{
  // A dummy function with the same name as this macro
}