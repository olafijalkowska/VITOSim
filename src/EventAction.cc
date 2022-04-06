
// $Id: EventAction.cc 04.04.2022 A. Fijalkowska $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
//
#include "EventAction.hh"

EventAction::EventAction()
{
    outputFile=OutputRoot::GetInstance();
    outputFile->CreateFile("ScintDetOutput.root");
}
 
EventAction::~EventAction()
{
	outputFile->SaveOutput();
}


void EventAction::BeginOfEventAction(const G4Event* anEvent)
{

	
}

double  EventAction::GetEnergy(G4HCofThisEvent* HCE, G4String detectorName)
{
    G4SDManager* SDmanager = G4SDManager::GetSDMpointer();
    G4int collectionId = SDmanager->GetCollectionID(detectorName);
    
	G4THitsMap<G4double>* eventMap = (G4THitsMap<G4double>*)(HCE->GetHC(collectionId));
	double energyDeposit = 0;
	if((*eventMap)[0] != 0)
	{
		energyDeposit = *((*eventMap)[0]);	    
	}
	//if(eventMap->entries() !=0)
	//	std::cout << detectorName << " " << collectionId << " " << energyDeposit << " " << eventMap->entries() << std::endl;
	return energyDeposit;
}

void EventAction::EndOfEventAction(const G4Event* anEvent)
{
	
  G4int eventID = anEvent->GetEventID();
  if( eventID % 1000 == 0 )
  {
    G4cout << "Finished Running Event # " << eventID << G4endl;
  }
  
    G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
    if(!HCE)
        return;
    
    int eventId = anEvent->GetEventID();    
	double enFront = GetEnergy(HCE, "scintFrontSD/eDep");
	double enRear2h = GetEnergy(HCE, "scintRear2hSD/eDep");
	double enRear4h = GetEnergy(HCE, "scintRear4hSD/eDep");
	double enRear6h = GetEnergy(HCE, "scintRear6hSD/eDep");
	double enRear8h = GetEnergy(HCE, "scintRear8hSD/eDep");
	double enRear10h = GetEnergy(HCE, "scintRear10hSD/eDep");
	double enRear12h = GetEnergy(HCE, "scintRear12hSD/eDep");
	
	outputFile->AddHit(eventId, enFront, enRear2h, enRear4h, enRear6h, enRear8h, enRear10h, enRear12h);
	
	
}

 
 

