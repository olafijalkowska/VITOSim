//
// $Id: DetectorConstruction.cc 12.16.2016 A. Fijalkowska $
//
/// \file DetectorConstruction.cc
/// \brief DetectorConstruction class
//
//
#include "DetectorConstruction.hh"
#include "G4PVPlacement.hh" 
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh" 
#include "G4NistManager.hh" 
#include "G4Element.hh" 
#include "G4Box.hh" 
#include "G4ThreeVector.hh" 
#include "globals.hh"
#include "Scintillator.hh"
#include "SampleHolder2.hh"
#include "CADMesh.hh"
#include "VITOMagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

DetectorConstruction::DetectorConstruction()
{
  worldLogic = 0L;
}



DetectorConstruction::~DetectorConstruction() {}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4VPhysicalVolume* worldPhys = ConstructWorld();
    ConstructScintillator();
    ConstructAluFrame();
    ConstructHolder();
    return worldPhys;
}


G4VPhysicalVolume* DetectorConstruction::ConstructWorld()
{

    G4double worldX = 1*m;
    G4double worldY = 1*m;
    G4double worldZ = 1*m;
    G4Material* vaccum = new G4Material("GalacticVacuum", 1., 1.01*g/mole,
                           CLHEP::universe_mean_density, 
                           kStateGas, 3.e-18*pascal, 2.73*kelvin);
  
    G4Box* worldSolid = new G4Box("worldSolid",worldX,worldY,worldZ);
    worldLogic = new G4LogicalVolume(worldSolid, vaccum,"worldLogic", 0,0,0);
                                             
    //worldLogic->SetVisAttributes(G4ktVisAttributes::Invisible);
    G4VPhysicalVolume* worldPhys = new G4PVPlacement(0, G4ThreeVector(), worldLogic, "world", 0, false, 0);
    return worldPhys;

}

void DetectorConstruction::ConstructScintillator()
{
    Scintillator* scint =  new Scintillator();
    G4ThreeVector pos(0,0,0);  
    scint->Place(0, pos, "scintDet", worldLogic, 0);
    scint->ConstructSDandField();
}

void DetectorConstruction::ConstructAluFrame()
{
   CADMesh* frameMesh = new CADMesh("../stl/alumin.stl",
                                    1*mm,
                                    G4ThreeVector(0*cm, 0*cm, 0*cm), 
                                    false);                                    
    G4VSolid* frameSolid =  frameMesh->TessellatedMesh(); 
    
    double z, atomicMass, density;
    int numberElements;
   	G4Element* Al = new G4Element("Aluminium","Al", z=13., atomicMass = 26.98*g/mole );
    G4Material* aluminum = new G4Material( "Aluminium", density= 2.7*g/cm3, numberElements=1 );
	aluminum->AddElement( Al, 1 );
    
    G4LogicalVolume* frameLog = new G4LogicalVolume(frameSolid, aluminum, "frameSolid");
    G4VisAttributes* frameVisAtt = new G4VisAttributes(G4Colour(1, 0, 0, 0.7));
	frameVisAtt->SetForceAuxEdgeVisible(true);// Can see outline when drawn with lines
	frameVisAtt->SetForceSolid(true);
	frameLog->SetVisAttributes(frameVisAtt);
	new G4PVPlacement(0, G4ThreeVector(0,0,0), frameLog, "framePhys",  worldLogic, 1,0);
}

void DetectorConstruction::ConstructHolder()
{
    SampleHolder2* holder =  new SampleHolder2();
    G4ThreeVector pos(0,0,0);  
    holder->Place(0, pos, "holderDet", worldLogic, 0);

}

void DetectorConstruction::ConstructSDandField() 
{
/*
	if (fField.Get() == 0)
    {
		G4MagneticField* magField= new VITOMagneticField("../field1Axial.txt", 
														 "../field1Radial.txt", 
														 "../field2Axial.txt", 
														 "../field2Radial.txt");
		fField.Put(magField);      
		//This is thread-local
		G4FieldManager* pFieldMgr = 
		G4TransportationManager::GetTransportationManager()->GetFieldManager();
		pFieldMgr->SetDetectorField(fField.Get());
		pFieldMgr->CreateChordFinder(fField.Get());  
	}  */
}






