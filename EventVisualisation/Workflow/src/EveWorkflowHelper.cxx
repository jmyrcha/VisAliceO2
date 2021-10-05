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

/// \file EveWorkflowHelper.cxx
/// \author julian.myrcha@cern.ch

#include <EveWorkflow/EveWorkflowHelper.h>
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "EveWorkflow/FileProducer.h"
#include "DataFormatsTRD/TrackTRD.h"

#include "ITStracking/IOUtils.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "SpacePoints/SpacePointsCalibParam.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#include "DetectorsBase/Propagator.h"
#include "TPCBase/ParameterElectronics.h"
#include "DataFormatsTPC/Defs.h"
#include "TPCFastTransform.h"
#include "TPCReconstruction/TPCFastTransformHelperO2.h"



using namespace o2::event_visualisation;


void EveWorkflowHelper::selectTracks(const CalibObjectsConst* calib,
                                     GID::mask_t maskCl, GID::mask_t maskTrk, GID::mask_t maskMatch)
{
  auto creator = [maskTrk, this](auto& trk, GID gid, float time, float) {
    if (!maskTrk[gid.getSource()]) {
      return true;
    }
    mTrackSet.trackGID.push_back(gid);
    mTrackSet.trackTime.push_back(time);
    return true;
  };
  this->mRecoCont.createTracksVariadic(creator);
}

void EveWorkflowHelper::draw(std::string jsonPath, int numberOfFiles, int numberOfTracks)
{
  prepareITSClusters();
  prepareMFTClusters();

  size_t nTracks = mTrackSet.trackGID.size();
  if (numberOfTracks != -1 && numberOfTracks < nTracks) {
    nTracks = numberOfTracks; // less than available
  }
  for (size_t it = 0; it < nTracks; it++) {
    const auto& gid = mTrackSet.trackGID[it];
    auto tim = mTrackSet.trackTime[it];

    switch(gid.getSource()) {
      case GID::TPC:
        drawTPC(gid, tim);
        break;
      case GID::MFT:
        drawMFT(gid, tim);
        break;
      case GID::MCH:
        drawMCH(gid, tim);
        break;
      case GID::ITS:
        drawITS(gid, tim);
        break;
      case GID::ITSTPCTOF:
        drawITSTPCTOF(gid, tim);
        break;
      case GID::ITSTPC:
        drawITSTPC(gid, tim);
        break;
      default:
        LOG(INFO) << "Track type " << gid.getSource() << " not handled";
    }

  }
  FileProducer producer(jsonPath, numberOfFiles);
  mEvent.toFile(producer.newFileName());
}

std::vector<PNT> EveWorkflowHelper::getTrackPoints(const o2::track::TrackPar& trc, float minR, float maxR, float maxStep)
{
  // adjust minR according to real track start fro track starting point
  float rMin = std::sqrt(trc.getX() * trc.getX() + trc.getY() * trc.getY());
  if (rMin > minR) {
    minR = rMin;
  }
  // prepare space points from the track param
  std::vector<PNT> pnts;
  int nSteps = std::max(2, int((maxR - minR) / maxStep));
  const auto prop = o2::base::Propagator::Instance();
  float xMin = trc.getX(), xMax = maxR * maxR - trc.getY() * trc.getY();
  if (xMax > 0) {
    xMax = std::sqrt(xMax);
  }
  //LOG(INFO) << "R: " << minR << " " << maxR << " || X: " << xMin << " " << xMax;
  float dx = (xMax - xMin) / nSteps;
  auto tp = trc;
  float dxmin = std::abs(xMin - tp.getX()), dxmax = std::abs(xMax - tp.getX());
  bool res = false;
  if (dxmin > dxmax) { //start from closest end
    std::swap(xMin, xMax);
    dx = -dx;
  }
  if (!prop->propagateTo(tp, xMin, false, 0.99, maxStep, o2::base::PropagatorF::MatCorrType::USEMatCorrNONE)) {
    return pnts;
  }
  auto xyz = tp.getXYZGlo();
  pnts.emplace_back(PNT{xyz.X(), xyz.Y(), xyz.Z()});
  for (int is = 0; is < nSteps; is++) {
    if (!prop->propagateTo(tp, tp.getX() + dx, false, 0.99, 999., o2::base::PropagatorF::MatCorrType::USEMatCorrNONE)) {
      return pnts;
    }
    xyz = tp.getXYZGlo();
    pnts.emplace_back(PNT{xyz.X(), xyz.Y(), xyz.Z()});
  }
  return pnts;
}




template <typename Functor>
void EveWorkflowHelper::addTrackToEvent(Functor source, GID gid, float trackTime, float dz)
{
  const auto& tr = source(gid);

  auto vTrack = mEvent.addTrack({.time = trackTime,
                                  .charge = tr.getCharge(),
                                  .PID = tr.getPID(),
                                  .startXYZ = {tr.getX(), tr.getY(), tr.getZ()},
                                  .phi = tr.getPhi(),
                                  .theta = tr.getTheta(),
                                  .source = (o2::dataformats::GlobalTrackID::Source)gid.getSource()});
  auto pnts = getTrackPoints(tr, minmaxR[gid.getSource()].first, minmaxR[gid.getSource()].second, 4);

  for (size_t ip = 0; ip < pnts.size(); ip++) {
    vTrack->addPolyPoint(pnts[ip][0], pnts[ip][1], pnts[ip][2] + dz);
  }
}

void EveWorkflowHelper::prepareITSClusters(std::string dictfile)
{
  o2::itsmft::TopologyDictionary dict;
  if (dictfile.empty()) {
    dictfile = o2::base::NameConf::getAlpideClusterDictionaryFileName(o2::detectors::DetID::ITS, "", "bin");
    dict.readBinaryFile(dictfile);
  }
  const auto& ITSClusterROFRec = mRecoCont.getITSClustersROFRecords();
  const auto& clusITS = mRecoCont.getITSClusters();
  if (clusITS.size() && ITSClusterROFRec.size()) {
    const auto& patterns = mRecoCont.getITSClustersPatterns();
    auto pattIt = patterns.begin();
    mITSClustersArray.reserve(clusITS.size());
    o2::its::ioutils::convertCompactClusters(clusITS, pattIt, mITSClustersArray, dict);
  }
}

void EveWorkflowHelper::prepareMFTClusters(std::string dictionaryFile)  // do we also have something as ITS...dict?
{
  o2::itsmft::TopologyDictionary dict;
  if (dictionaryFile.empty()) {
    dictionaryFile = o2::base::NameConf::getAlpideClusterDictionaryFileName(o2::detectors::DetID::MFT, "", "bin");
    dict.readBinaryFile(dictionaryFile);
  }
  const auto& MFTClusterROFRec = this->mRecoCont.getMFTClustersROFRecords();
  const auto& clusMFT = this->mRecoCont.getMFTClusters();
  if (clusMFT.size() && MFTClusterROFRec.size()) {
    const auto& patterns = this->mRecoCont.getMFTClustersPatterns();
    auto pattIt = patterns.begin();
    this->mMFTClustersArray.reserve(clusMFT.size());
    o2::its::ioutils::convertCompactClusters(clusMFT, pattIt, this->mMFTClustersArray, dict);  // it is its not mft
  }
}


void EveWorkflowHelper::drawITSTPC(GID gid, float trackTime)
{
  const auto& track = mRecoCont.getTPCITSTrack(gid);
  auto pnts = getTrackPoints(track, minmaxR[gid.getSource()].first, minmaxR[gid.getSource()].second, 4);
  addTrackToEvent([this, trackTime](GID gid) { return mRecoCont.getTPCITSTrack(gid); }, trackTime, 0.);
  GID gidTPC = track.getRefTPC();
  GID gidITS = track.getRefITS();
  drawITSClusters(gidITS, trackTime);
  drawTPCClusters(gidTPC, trackTime);
}

void EveWorkflowHelper::drawITSTPCTOF(GID gid, float trackTime)
{
  const auto& track = mRecoCont.getITSTPCTOFTrack(gid);
  addTrackToEvent([this, trackTime](GID gid) { return mRecoCont.getITSTPCTOFTrack(gid); }, trackTime, 0.);
  GID gidTPC = track.getRefTPC();
  GID gidITS = track.getRefITS();
  drawITSClusters(gidITS, trackTime);
  drawTPCClusters(gidTPC, trackTime);
}



void EveWorkflowHelper::drawITSClusters(GID gid, float trackTime)
{
    LOG(INFO) << "+++++++++++++ drawITSClusters" ;
  const auto& trc = mRecoCont.getITSTrack(gid);
  auto refs = mRecoCont.getITSTracksClusterRefs();
  //int entry0 = trc.getClusterEntry(gid.getIndex()); // TODO correct?
  int ncl = trc.getNumberOfClusters();
  for (int icl = 0; icl < ncl; icl++) {
    const auto& pnt = mITSClustersArray[refs[icl]];
    drawPoint(pnt, trackTime);
  }
    LOG(INFO) << "------------- drawITSClusters" ;
}



std::vector<PNT> EveWorkflowHelper::getMFTTrackPoints(o2::mft::TrackMFT &mftTrack, float maxStep)
{
    LOG(INFO) << "+++++++++++++ EveWorkflowHelper::getMFTTrackPoints" ;
  std::vector<PNT> pnts;                                  // list of created points
  auto noOfClusters = mftTrack.getNumberOfPoints();       // number of clusters in MFT Track
  auto offset = mftTrack.getExternalClusterIndexOffset(); // TODO first cluster of that track on list of mft cluster refs
  auto refs = mRecoCont.getMFTTracksClusterRefs();        // list of references to clusters
  auto MFTClusters = mRecoCont.getMFTClusters();          // list of all clusters

  std::vector<float> zCoordinates ;
    LOG(INFO) << "+++++++++++++ EveWorkflowHelper::getMFTTrackPoints-1" ;

  for (int icl = noOfClusters - 1; icl > -1; --icl) {
    const auto& pnt = mMFTClustersArray[refs[icl]];
    zCoordinates.push_back(pnt.getZ());
  }

  if(zCoordinates.size())                       // should be, just safe
  {
    auto minZ = *(std::min_element(zCoordinates.begin(), zCoordinates.end()));  // TODO can we assume that order is preserved?
    auto maxZ = *(std::max_element(zCoordinates.begin(), zCoordinates.end()));
      LOG(INFO) << "+++++++++++++++++++++++++++++++++++++++++++++++++++=minZ " << minZ;
      LOG(INFO) << zCoordinates[0];
      LOG(INFO) << "+++++++++++++++++++++++++++++++++++++++++++++++++++=maxZ "<< maxZ;


    for(float z = minZ; z < maxZ; z+= maxStep)
      mftTrack.propagateToZlinear(z);         // TODO max_z should be global or cluster coordinate
      //auto gloXYZ = mMFTGeom->getMatrixL2G(pnt.getSensorID()) * pnt.getXYZ();
    }

    LOG(INFO) << "------------- EveWorkflowHelper::getMFTTrackPoints" ;
  return pnts;
}



// TPC cluseters for given TPC track (gid)
//  time should be in time bins
void EveWorkflowHelper::drawTPCClusters(GID gid, float trackTime)
{
    //LOG(INFO) << "++++++++++++++++ EveWorkflowHelper::drawTPCClusters" ;
    LOG(INFO) << "++++++++++++++++ EveWorkflowHelper::drawTPCClusters" << gid;
  const auto& trc = mRecoCont.getTPCTrack(gid);

  auto mTPCTracksClusIdx = mRecoCont.getTPCTracksClusterRefs();
  auto mTPCClusterIdxStruct = &mRecoCont.getTPCClusters();
  const auto& elParam = o2::tpc::ParameterElectronics::Instance();

  auto mTPCTimeBinMUS = elParam.ZbinWidth;
  float clusterTimeBinOffset = trackTime / mTPCTimeBinMUS;
  std::unique_ptr<gpu::TPCFastTransform> fastTransform = (o2::tpc::TPCFastTransformHelperO2::instance()->create(0));
  auto mFastTransform = std::move(fastTransform);

  // store the TPC cluster positions
  for (int iCl = trc.getNClusterReferences(); iCl--;) {
      uint8_t sector, row;     // TODO - how to set sector ???

      const auto& clTPC = trc.getCluster(mTPCTracksClusIdx, iCl, *mTPCClusterIdxStruct, sector, row);
      //clTPC.getTime() it is in time beans should be converted to ms by multiply   track time may be different
      //
      const float TB2MUSEC = o2::constants::lhc::LHCOrbitMUS / o2::constants::lhc::LHCMaxBunches * 8;

      //sector = clTPC.getPad();    // TODO - how to get sector ???
      float clTPCX;
      std::array<float, 2> clTPCYZ;
      mFastTransform->TransformIdeal(sector, row, clTPC.getPad(), clTPC.getTime(), clTPCX, clTPCYZ[0], clTPCYZ[1], clusterTimeBinOffset);  // which time should used here? this olso time means
      sector %= o2::tpc::SECTORSPERSIDE;

      double xyz[] = {o2::tpc::param::RowX[row],clTPCYZ[0],clTPCYZ[1]};
      o2::math_utils::rotateZd(o2::tpc::param::RowX[row], xyz[0], xyz[1], o2::math_utils::sector2Angle(sector) );
      mEvent.addCluster(xyz[0], xyz[1], xyz[2], trackTime);
  }
    //LOG(INFO) << "------------- EveWorkflowHelper::drawTPCClusters" ;
}

void EveWorkflowHelper::drawMFTClusters(GID gid, float trackTime)
{

    LOG(INFO) << "+++++++++++++++ drawMFTClusters";
  const auto& mftTrack = mRecoCont.getMFTTrack(gid);
  auto noOfClusters = mftTrack.getNumberOfPoints();       // number of clusters in MFT Track
  auto refs = mRecoCont.getMFTTracksClusterRefs();        // list of references to clusters
  for (int icl = noOfClusters - 1; icl > -1; --icl) {
    const auto& pnt = mMFTClustersArray[refs[icl]];
    auto gloXYZ = mMFTGeom->getMatrixL2G(pnt.getSensorID()) * pnt.getXYZ();
    float xyz[] = {gloXYZ.X(), gloXYZ.Y(), gloXYZ.Z()};
    drawPoint(xyz, trackTime);
  }
    LOG(INFO) << "-------------- drawMFTClusters";

}

void EveWorkflowHelper::drawTPC(GID gid, float trackTime)
{

    LOG(INFO) << "+++++++++++++++ drawTPC";
  const auto& tr = mRecoCont.getTPCTrack(gid);
  auto vTrack = mEvent.addTrack({.time = static_cast<float>(trackTime * 8 * o2::constants::lhc::LHCBunchSpacingMUS),
                                 .charge = tr.getCharge(),
                                 .PID = tr.getPID(),
                                 .startXYZ = {tr.getX(), tr.getY(), tr.getZ()},
                                 .phi = tr.getPhi(),
                                 .theta = tr.getTheta(),
                                 .source = GID::TPC});
  auto pnts = getTrackPoints(tr, minmaxR[gid.getSource()].first, minmaxR[gid.getSource()].second, 4);
  float dz = 0.0;
  for (size_t ip = 0; ip < pnts.size(); ip++) {
    vTrack->addPolyPoint(pnts[ip][0], pnts[ip][1], pnts[ip][2] + dz);
  }
  drawTPCClusters(gid, trackTime);
    LOG(INFO) << "----------------- drawTPC";
}

void EveWorkflowHelper::drawITS(GID gid, float trackTime)
{
    LOG(INFO) << "+++++++++++++++ drawITS";
  addTrackToEvent([this, trackTime](GID gid) { return mRecoCont.getITSTrack(gid); }, trackTime, 0.);
  drawITSClusters(gid, trackTime);
    LOG(INFO) << "---------------- drawITS";
}



void EveWorkflowHelper::drawMFT(GID gid, float trackTime) {
    LOG(INFO) << "++++++++++++++++++++++++++drawMFT ";
    auto tr = mRecoCont.getMFTTrack(gid);
    auto vTrack = mEvent.addTrack({.time = static_cast<float>(trackTime * 8 * o2::constants::lhc::LHCBunchSpacingMUS),
                                          .charge = (int)tr.getCharge(),
                                          .PID = o2::track::PID::Muon,
                                          .startXYZ = {(float)tr.getX(), (float)tr.getY(), (float)tr.getZ()},
                                          .phi = (float)tr.getPhi(),
                                          .theta = (float)tr.getTanl(),
                                          .source = GID::MFT});
    auto pnts = getMFTTrackPoints(tr,  4);
    float dz = 0.0;
    for (size_t ip = 0; ip < pnts.size(); ip++) {
       vTrack->addPolyPoint(pnts[ip][0], pnts[ip][1], pnts[ip][2] + dz);
    }
    drawMFTClusters(gid, trackTime);
    LOG(INFO) << "-----------------------------drawMFT ";
}


void EveWorkflowHelper::drawMCH(GID gid, float trackTime) {
    LOG(INFO) << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++MCH ";
}

EveWorkflowHelper::EveWorkflowHelper()
{
  this->mMFTGeom = o2::mft::GeometryTGeo::Instance();
  this->mMFTGeom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));
}







