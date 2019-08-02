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
/// \file    EventManager.h
/// \author  Jeremi Niedziela
/// \author julian.myrcha@cern.ch
/// \author p.nowakowski@cern.ch

#ifndef ALICE_O2_EVENTVISUALISATION_VIEW_EVENTMANAGER_H
#define ALICE_O2_EVENTVISUALISATION_VIEW_EVENTMANAGER_H

#include "CCDB/Manager.h"

#include <TEveElement.h>
#include <TEveEventManager.h>

#include <string>
#include <EventVisualisationDataConverter/VisualisationEvent.h>
#include <EventVisualisationBase/DataInterpreter.h>
#include <EventVisualisationBase/DataReader.h>


namespace o2
{
namespace event_visualisation
{


/// EventManager is a singleton class managing event loading.
///
/// This class is a hub for data visualisation classes, providing them with objects of requested type
/// (Raw data, hits, digits, clusters, ESDs, AODs...). It is a role of detector-specific data macros to
/// interpret data from different formats as visualisation objects (points, lines...) and register them
/// for drawing in the MultiView.

class DataSource;

class EventManager : public TEveEventManager
{
public:
    enum EDataSource{
      SourceOnline,   ///< Online reconstruction is a source of events
      SourceOffline,  ///< Local files are the source of events
      SourceHLT       ///< HLT reconstruction is a source of events
    };

    /// Returns an instance of EventManager
    static EventManager& getInstance();

    /// Setter of the current data source type
    inline void setDataSourceType(EDataSource source)   { mCurrentDataSourceType = source; }
    /// Setter of the current data source path
    inline void setDataSourcePath(const TString& path)  { dataPath = path;}
    /// Sets the CDB path in CCDB Manager
    inline void setCdbPath(const TString& path)  {
      o2::ccdb::Manager::Instance()->setDefaultStorage(path.Data());
    }

    Int_t getCurrentEvent() const {return currentEvent;}
    DataSource *getDataSource() const {return dataSource;}

    void setDataSource(DataSource *dataSource)          { this->dataSource = dataSource; }

    void Open() override ;
    void GotoEvent(Int_t /*event*/) override ;
    void NextEvent() override ;
    void PrevEvent() override ;
    void Close() override ;

    void AfterNewEventLoaded() override;

    void AddNewEventCommand(const TString& cmd) override ;
    void RemoveNewEventCommand(const TString& cmd) override ;
    void ClearNewEventCommands() override ;

    void registerDetector(DataReader *reader, DataInterpreter *interpreter, EVisualisationGroup type);

private:
    static EventManager *instance;
    DataInterpreter* dataInterpreters[EVisualisationGroup::NvisualisationGroups];
    DataReader* dataReaders[EVisualisationGroup::NvisualisationGroups];


    /// store lists of visualisation element in current event (row, clusters ...)
    TEveElementList* dataTypeLists[EVisualisationDataType::NdataTypes];

    /// store user request to see data for given type  (row, clusters ...)
    TEveElementList* mRequestedDataTypes[EVisualisationDataType::NdataTypes];

    EDataSource mCurrentDataSourceType = EDataSource::SourceOffline;
    DataSource *dataSource = nullptr;
    TString dataPath = "";
    Int_t currentEvent = 0;

    /// Default constructor
    EventManager();
    /// Default destructor
    ~EventManager() final;
    /// Deleted copy constructor
    EventManager(EventManager const&) = delete;
    /// Deleted assignemt operator
    void operator=(EventManager const&) = delete;

    //Display visualisation event
    void displayVisualisationEvent(VisualisationEvent &event, const std::string &detectorName);
};

}
}

#endif


