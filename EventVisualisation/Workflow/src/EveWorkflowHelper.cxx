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
#include "MFTTracking/IOUtils.h"
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
  // adjust minR according to real track start from track starting point
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

  float dx = (xMax - xMin) / nSteps;
  auto tp = trc;
  float dxmin = std::abs(xMin - tp.getX()), dxmax = std::abs(xMax - tp.getX());

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
    o2::mft::ioutils::convertCompactClusters(clusMFT, pattIt, this->mMFTClustersArray, dict);
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
  const auto& trc = mRecoCont.getITSTrack(gid);
  auto refs = mRecoCont.getITSTracksClusterRefs();
  int ncl = trc.getNumberOfClusters();
  for (int icl = 0; icl < ncl; icl++) {
    const auto& pnt = mITSClustersArray[refs[icl]];
    const auto glo = mITSGeom->getMatrixT2G(pnt.getSensorID()) *pnt.getXYZ() ;
    drawPoint(glo.X(),glo.Y(),glo.Z(), trackTime);
  }
}



std::vector<PNT> EveWorkflowHelper::getMFTTrackPoints(o2::mft::TrackMFT &mftTrack, float maxStep)
{
  std::vector<PNT> pnts;                                  // list of created points
  auto noOfClusters = mftTrack.getNumberOfPoints();       // number of clusters in MFT Track
  auto offset = mftTrack.getExternalClusterIndexOffset(); // first external cluster index offset:
  auto refs = mRecoCont.getMFTTracksClusterRefs();        // list of references to clusters, offset:offset+no
  for (int mcl = noOfClusters - 1; mcl > -1; --mcl) {
    const auto& pnt = mMFTClustersArray[refs[mcl+offset]];
    auto gloXYZ = mMFTGeom->getMatrixL2G(pnt.getSensorID()) * pnt.getXYZ();
    pnts.emplace_back(PNT{gloXYZ.X(), gloXYZ.Y(), gloXYZ.Z()});
  }
  return pnts;
}



// TPC cluseters for given TPC track (gid)
void EveWorkflowHelper::drawTPCClusters(GID gid, float trackTime)
{
  const auto& trc = mRecoCont.getTPCTrack(gid);
  auto mTPCTracksClusIdx = mRecoCont.getTPCTracksClusterRefs();
  auto mTPCClusterIdxStruct = &mRecoCont.getTPCClusters();
  const auto& elParam = o2::tpc::ParameterElectronics::Instance();

  auto mTPCTimeBinMUS = elParam.ZbinWidth;
  float clusterTimeBinOffset = trackTime / mTPCTimeBinMUS;
  std::unique_ptr<gpu::TPCFastTransform> fastTransform = (o2::tpc::TPCFastTransformHelperO2::instance()->create(0));

  // store the TPC cluster positions
  for (int iCl = trc.getNClusterReferences(); iCl--;) {
      uint8_t sector, row;
      const auto& clTPC = trc.getCluster(mTPCTracksClusIdx, iCl, *mTPCClusterIdxStruct, sector, row);
      sector %= o2::tpc::SECTORSPERSIDE;

    std::array<float, 3> xyz;
    fastTransform->TransformIdeal(sector, row, clTPC.getPad(), clTPC.getTime(), xyz[0], xyz[1], xyz[2], clusterTimeBinOffset);
    o2::math_utils::rotateZ(xyz, o2::math_utils::sector2Angle(sector%o2::tpc::SECTORSPERSIDE) );
    mEvent.addCluster(xyz[0], xyz[1], xyz[2], trackTime);
  }
}

void EveWorkflowHelper::drawMFTClusters(GID gid, float trackTime)
{
  const auto& mftTrack = mRecoCont.getMFTTrack(gid);
  auto noOfClusters = mftTrack.getNumberOfPoints();       // number of clusters in MFT Track
  auto offset = mftTrack.getExternalClusterIndexOffset(); // first external cluster index offset:
  auto refs = mRecoCont.getMFTTracksClusterRefs();        // list of references to clusters, offset:offset+no
  for (int icl = noOfClusters - 1; icl > -1; --icl) {
    const auto& pnt = mMFTClustersArray[refs[offset+icl]];
    auto gloXYZ = mMFTGeom->getMatrixL2G(pnt.getSensorID()) * pnt.getXYZ();
    drawPoint(gloXYZ.X(), gloXYZ.Y(), gloXYZ.Z(), trackTime);
  }
}

void EveWorkflowHelper::drawTPC(GID gid, float trackTime)
{
  //LOG(INFO) << "+++++++++++++++ drawTPC";
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
  //LOG(INFO) << "----------------- drawTPC";
}

void EveWorkflowHelper::drawITS(GID gid, float trackTime)
{
  //LOG(INFO) << "+++++++++++++++ drawITS";
  addTrackToEvent([this, trackTime](GID gid) { return mRecoCont.getITSTrack(gid); }, trackTime, 0.);
  drawITSClusters(gid, trackTime);
  //LOG(INFO) << "---------------- drawITS";
}



void EveWorkflowHelper::drawMFT(GID gid, float trackTime) {
    //LOG(INFO) << "++++++++++++++++++++++++++drawMFT ";
    auto tr = mRecoCont.getMFTTrack(gid);
    auto vTrack = mEvent.addTrack({.time = static_cast<float>(trackTime),
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
    //LOG(INFO) << "-----------------------------drawMFT ";
}


void EveWorkflowHelper::drawMCH(GID gid, float trackTime) {
    LOG(INFO) << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++MCH ";
}

EveWorkflowHelper::EveWorkflowHelper()
{
  this->mMFTGeom = o2::mft::GeometryTGeo::Instance();
  this->mMFTGeom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::L2G));

  this->mITSGeom = o2::its::GeometryTGeo::Instance();
  this->mITSGeom->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot,o2::math_utils::TransformType::L2G));
}







