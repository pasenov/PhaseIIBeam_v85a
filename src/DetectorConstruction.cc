//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DetectorConstruction.cc 68348 2013-03-22 10:00:19Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ElectricFieldSetup.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <iomanip>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4AutoDelete.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4NistMaterialBuilder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),fPBoxW(0), fPBox1(0), fPBox2(0), fPBoxBPIX12(0), fPBoxBPIX34(0), fLBoxW(0), fLBox1(0), fLBox2(0), fLBoxBPIX12(0), fLBoxBPIX34(0), fMaterialW(0), fMaterialD(0), fMaterialBPIX(0), fMaterialScint(0), fMaterialAl(0), fDetectorMessenger(0)
{
  zOff = 13.5*cm;
  fWorldSizeX   = 47.0*cm; 
  fWorldSizeY   = 47.0*cm; 
  fWorldSizeZ   = 47.0*cm;
  fDet1SizeX = 102700.0*um;
  fDet1SizeY = 94108.0*um;
  fDet1SizeZ = 320.0*um;
  fDet2SizeX = 102700.0*um;
  fDet2SizeY = 94108.0*um;
  fDet2SizeZ = 320.0*um;
  fBPIXSizeX = 66.6*mm;
  fBPIXSizeY = 25.0*mm;
  fBPIXSizeZ = 460.0*um;
  fBPIXSenSizeX = 66.6*mm;
  fBPIXSenSizeY = 25.0*mm;
  fBPIXSenSizeZ = 285.0*um;
  fBPIXROCSizeX = 66.6*mm;
  fBPIXROCSizeY = 25.0*mm;
  fBPIXROCSizeZ = 460.0*um - 285.0*um;
  fScintSizeX = 67.0*mm;
  fScintSizeY = 40.0*mm;
  fScintSizeZ = 2.0*mm;
  fAlSizeX = 66.6*mm;
  fAlSizeY = 25.0*mm;
  fAlSizeZ = 50.0*um;

  fDist = 2.0*mm;
  fDist2 = 1.0*cm;
  fInterArmDeltaZ = fDist2 + fDet1SizeZ + fDist + fDet2SizeZ + fDist2;

  fStrip1Depth = fStrip2Depth = 290.0*um;
  fStrip1Length = fStrip2Length = (102700.0 - 2*1368.0)*um;
  fStripDist = 68.0*um;
  fStripWidth = 22.0*um;
  fStripPitch = fStripDist + fStripWidth;

  fPixSensorSizeX = 64.8*mm;
  fPixSensorSizeY = 18.6*mm;

  fPixelPitchX = 150.0*um;
  fPixelDoublePitchX = 2*150.0*um;
  fPixelPitchY = 100.0*um;
  fPixelDoublePitchY = 2*100.0*um;
  fPixelDepth = 285.0*um;

  fRestZ = fDet1SizeZ - fStrip1Depth;

  fDUTSizeX = fDet1SizeX;
  fDUTSizeY = fDet1SizeY;
  fDUTSizeZ = fDet1SizeZ + fDist + fDet2SizeZ;

  pos_EndArm1Abs = 1.8*cm - fBPIXSizeZ/2; //pos_EndArm1Abs = fInterArmDeltaZ/2;
  pos_BeginningArm2Abs = 1.8*cm - fBPIXSizeZ/2; //pos_BeginningArm2Abs = fInterArmDeltaZ/2;

  fPotStrip1 = fPotStrip2 = 0*kilovolt;
  fPotBackplane1 = -0.4*kilovolt;
  fPotBackplane2 = -0.15*kilovolt;

  DefineMaterials();
  SetMaterialW("Air");
  SetMaterialD("Silicon"); 
  SetMaterialBPIX("Silicon"); 
  SetMaterialScint("Scintillator");
  SetMaterialAl("Aluminium");
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  //
  // define Elements
  //
  G4double z,a;
  
  G4Element* H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"   ,"C" , z= 6., a=  12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"   ,"O" , z= 8., a=  16.00*g/mole);
  
  //
  // define materials
  //
  G4double density;
  G4int ncomponents, natoms;
  G4double fractionmass;  

  G4Material* H2O = 
  new G4Material("Water", density= 1.000*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  G4Material* CO2 = 
  new G4Material("CarbonicGas", density= 1.87*mg/cm3, ncomponents=2, kStateGas, 273.15*kelvin, 1*atmosphere);
  CO2->AddElement(C, natoms=1);
  CO2->AddElement(O, natoms=2);
  CO2->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  G4Material* vapor = 
  new G4Material("Water_vapor", density= 1.000*mg/cm3, ncomponents=2);
  vapor->AddElement(H, natoms=2);
  vapor->AddElement(O, natoms=1);
  vapor->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  new G4Material("Carbon"     , z=6.,  a= 12.01*g/mole, density= 2.267*g/cm3);
  new G4Material("Aluminium"  , z=13., a= 26.98*g/mole, density= 2.700*g/cm3);
  new G4Material("Silicon"    , z=14., a= 28.09*g/mole, density= 2.330*g/cm3);
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  new G4Material("Iron"       , z=26., a= 55.85*g/mole, density= 7.870*g/cm3);  
  new G4Material("Germanium"  , z=32., a= 72.61*g/mole, density= 5.323*g/cm3);
  new G4Material("Tungsten"   , z=74., a=183.85*g/mole, density= 19.30*g/cm3);
  new G4Material("Lead"       , z=82., a=207.19*g/mole, density= 11.35*g/cm3);
  new G4Material("Nitrogen"   , z=7.,  a= 14.01*g/mole, density= 1.145*mg/cm3);
  
  G4Material* ArgonGas =   
  new G4Material("ArgonGas"   , z=18., a=39.948*g/mole, density= 1.782*mg/cm3,
                 kStateGas, 273.15*kelvin, 1*atmosphere);

  G4Material* Air = 
  new G4Material("Air", density= 1.290*mg/cm3, ncomponents=4);
  Air->AddElement(N, fractionmass=78.08*perCent);
  Air->AddElement(O, fractionmass=20.95*perCent);
  Air->AddMaterial(ArgonGas, fractionmass=0.93*perCent);
  Air->AddMaterial(CO2, fractionmass=0.04*perCent);

  G4Material* Scintillator = 
  new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
  Scintillator->AddElement(C, natoms=9);
  Scintillator->AddElement(H, natoms=10);
  Scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  fDUTangleX = 0.*deg;
  fBPIXangleX = 0.*deg;
  fBPIXangleY = 0.*deg;
  fScintangleX = 0.*deg;
  fScintangleY = 0.*deg;

  G4RotationMatrix* myRotationDUT = new G4RotationMatrix();
  myRotationDUT->rotateX(0.*deg);
  myRotationDUT->rotateY(0.*deg);
  myRotationDUT->rotateZ(0.*deg);

  G4RotationMatrix* myRotationBPIX = new G4RotationMatrix();
  myRotationBPIX->rotateX(fBPIXangleX);
  myRotationBPIX->rotateY(fBPIXangleY);
  myRotationBPIX->rotateZ(0.*deg);

  G4RotationMatrix* myRotationScint = new G4RotationMatrix();
  myRotationScint->rotateX(fScintangleX);
  myRotationScint->rotateY(fScintangleY);
  myRotationScint->rotateZ(0.*deg);

  G4Box*
  sBoxW = new G4Box("World",                              //its name
                   fWorldSizeX/2,fWorldSizeY/2,fWorldSizeZ/2);       //its dimensions

  fLBoxW = new G4LogicalVolume(sBoxW,                        //its shape
                             fMaterialW,                    //its material
                             fMaterialW->GetName());        //its name

  fPBoxW = new G4PVPlacement(0,                          //no rotation
                           G4ThreeVector(),           //at (0,0,0)
                           fLBoxW,                       //its logical volume
                           fMaterialW->GetName(),        //its name
                           0,                           //its mother  volume
                           false,                       //no boolean operation
                           0);                          //copy number

  // DUT 2S

  G4Box*
  sBoxDUT = new G4Box("DUT",				     //its name
				 fDUTSizeX/2,fDUTSizeY/2,fDUTSizeZ/2);   //its dimensions
  

  fLBoxDUT = new G4LogicalVolume(sBoxDUT,                   //its shape
			     fMaterialW,                   //its material
			     fMaterialW->GetName());       //its name
  
  pos_xDUT =  0.0*cm;
  pos_yDUT =  0.0*cm;
  pos_zDUT =  -fWorldSizeZ/2 + 13.5*cm;
  
  fPBoxDUT =  new G4PVPlacement(myRotationDUT,			       //rotation
		    G4ThreeVector(pos_xDUT, pos_yDUT, pos_zDUT),	//at (pos_xDUT, pos_yDUT, pos_zDUT)
                    "DUT",                  //its name
                    fLBoxDUT,                       //its logical volume
                    fPBoxW,                           //its mother  volume
                    false,                       //no boolean operation
                    0);                          //copy number

  // 2S sensor 1

  G4Box*
  sBox1 = new G4Box("Detector1",				     //its name
				 fDet1SizeX/2,fDet1SizeY/2,fDet1SizeZ/2);   //its dimensions
  

  fLBox1 = new G4LogicalVolume(sBox1,                      //its shape
			     fMaterialD,                   //its material
			     fMaterialD->GetName());        //its name

  pos_x1 =  0.0*cm;
  pos_y1 =  0.0*cm;
  pos_z1 =  -fDist/2 - fDet1SizeZ/2;
  
  fPBox1 =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_x1, pos_y1, pos_z1),	       //at (pos_x1, pos_y1, pos_z1)
                    "Detector1",                  //its name
                    fLBox1,                       //its logical volume
                    fPBoxDUT,                           //its mother  volume
                    false,                       //no boolean operation
                    0);                          //copy number

  // 2S sensor 2

  G4Box*
  sBox2 = new G4Box("Detector2",				     //its name
				 fDet2SizeX/2,fDet2SizeY/2,fDet2SizeZ/2);   //its dimensions
  

  fLBox2 = new G4LogicalVolume(sBox2,                      //its shape
			     fMaterialD,                   //its material
			     fMaterialD->GetName());        //its name
  
  pos_x2 =  0.0*cm;
  pos_y2 =  0.0*cm;
  pos_z2 =  fDist/2 + fDet2SizeZ/2;
  
  fPBox2 =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_x2, pos_y2, pos_z2),	       //at (pos_x2, pos_y2, pos_z2)
                    "Detector2",                  //its name
                    fLBox2,                       //its logical volume
                    fPBoxDUT,                           //its mother  volume
                    false,                       //no boolean operation
                    0);                          //copy number

  // Pixel modules 1, 2

  G4Box*
  sBoxBPIX12 = new G4Box("Pixel 1",				     //its name
				 fBPIXSizeX/2,2*fBPIXSizeY/2,fBPIXSizeZ/2);   //its dimensions
  

  fLBoxBPIX12 = new G4LogicalVolume(sBoxBPIX12,                      //its shape
			     fMaterialBPIX,                   //its material
			     fMaterialBPIX->GetName());        //its name
  
  pos_xBPIX12 =  0.0*cm;
  pos_yBPIX12 =  0.0*cm;
  pos_zBPIX12 = -fWorldSizeZ/2 + 11.7*cm; //G4double pos_zBPIX12 =  -pos_EndArm1Abs - fBPIXSizeZ/2;
  
  fPBoxBPIX12 =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xBPIX12, pos_yBPIX12, pos_zBPIX12),	       //at its position
                    "Pixel 1",                  //its name
                    fLBoxBPIX12,                       //its logical volume
                    fPBoxW,                           //its mother  volume
                    false,                       //no boolean operation
                    0);  

  // Sensors of pixel modules 1, 2

  G4Box*
  sBoxBPIX12Sen = new G4Box("Pixel 1 Sensor",				     //its name
				 fBPIXSenSizeX/2,2*fBPIXSenSizeY/2,fBPIXSenSizeZ/2);   //its dimensions
  

  fLBoxBPIX12Sen = new G4LogicalVolume(sBoxBPIX12Sen,                      //its shape
			     fMaterialBPIX,                   //its material
			     fMaterialBPIX->GetName());        //its name
  
  pos_xBPIX12Sen =  0.0*cm;
  pos_yBPIX12Sen =  0.0*cm;
  pos_zBPIX12Sen = -fBPIXSizeZ/2 + fBPIXSenSizeZ/2;
  
  fPBoxBPIX12Sen =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xBPIX12Sen, pos_yBPIX12Sen, pos_zBPIX12Sen),	       //at its position
                    "Pixel 1 Sensor",                  //its name
                    fLBoxBPIX12Sen,                       //its logical volume
                    fPBoxBPIX12,                           //its mother  volume
                    false,                       //no boolean operation
                    0); 

  // ROC of pixel modules 1, 2

  G4Box*
  sBoxBPIX12ROC = new G4Box("Pixel 1 ROC",				     //its name
				 fBPIXROCSizeX/2,2*fBPIXROCSizeY/2,fBPIXROCSizeZ/2);   //its dimensions
  

  fLBoxBPIX12ROC = new G4LogicalVolume(sBoxBPIX12ROC,                      //its shape
			     fMaterialBPIX,                   //its material
			     fMaterialBPIX->GetName());        //its name
  
  pos_xBPIX12ROC =  0.0*cm;
  pos_yBPIX12ROC =  0.0*cm;
  pos_zBPIX12ROC = fBPIXSizeZ/2 - fBPIXROCSizeZ/2;
  
  fPBoxBPIX12ROC =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xBPIX12ROC, pos_yBPIX12ROC, pos_zBPIX12ROC),	       //at its position
                    "Pixel 1 ROC",                  //its name
                    fLBoxBPIX12ROC,                       //its logical volume
                    fPBoxBPIX12ROC,                           //its mother  volume
                    false,                       //no boolean operation
                    0); 

  // Pixel modules 3, 4

  G4Box*
  sBoxBPIX34 = new G4Box("Pixel 2",				     //its name
				 fBPIXSizeX/2,2*fBPIXSizeY/2,fBPIXSizeZ/2);   //its dimensions
  

  fLBoxBPIX34 = new G4LogicalVolume(sBoxBPIX34,                      //its shape
			     fMaterialBPIX,                   //its material
			     fMaterialBPIX->GetName());        //its name
  
  pos_xBPIX34 =  0.0*cm;
  pos_yBPIX34 =  0.0*cm;
  pos_zBPIX34 = -fWorldSizeZ/2 + 15.3*cm; //G4double pos_zBPIX34 =  pos_BeginningArm2Abs + fBPIXSizeZ/2;
  
  fPBoxBPIX34 =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xBPIX34, pos_yBPIX34, pos_zBPIX34),	       //at its position
                    "Pixel 2",                  //its name
                    fLBoxBPIX34,                       //its logical volume
                    fPBoxW,                           //its mother  volume
                    false,                       //no boolean operation
                    0);

  // Sensors of pixel modules 3, 4

  G4Box*
  sBoxBPIX34Sen = new G4Box("Pixel 2 Sensor",				     //its name
				 fBPIXSenSizeX/2,2*fBPIXSenSizeY/2,fBPIXSenSizeZ/2);   //its dimensions
  

  fLBoxBPIX34Sen = new G4LogicalVolume(sBoxBPIX34Sen,                      //its shape
			     fMaterialBPIX,                   //its material
			     fMaterialBPIX->GetName());        //its name
  
  pos_xBPIX34Sen =  0.0*cm;
  pos_yBPIX34Sen =  0.0*cm;
  pos_zBPIX34Sen = -fBPIXSizeZ/2 + fBPIXSenSizeZ/2;
  
  fPBoxBPIX34Sen =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xBPIX34Sen, pos_yBPIX34Sen, pos_zBPIX34Sen),	       //at its position
                    "Pixel 2 Sensor",                  //its name
                    fLBoxBPIX34Sen,                       //its logical volume
                    fPBoxBPIX34,                           //its mother  volume
                    false,                       //no boolean operation
                    0); 

  // ROC of pixel modules 3, 4

  G4Box*
  sBoxBPIX34ROC = new G4Box("Pixel 2 ROC",				     //its name
				 fBPIXROCSizeX/2,2*fBPIXROCSizeY/2,fBPIXROCSizeZ/2);   //its dimensions
  

  fLBoxBPIX34ROC = new G4LogicalVolume(sBoxBPIX34ROC,                      //its shape
			     fMaterialBPIX,                   //its material
			     fMaterialBPIX->GetName());        //its name
  
  pos_xBPIX34ROC =  0.0*cm;
  pos_yBPIX34ROC =  0.0*cm;
  pos_zBPIX34ROC = fBPIXSizeZ/2 - fBPIXROCSizeZ/2;
  
  fPBoxBPIX34ROC =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xBPIX34ROC, pos_yBPIX34ROC, pos_zBPIX34ROC),	       //at its position
                    "Pixel 2 ROC",                  //its name
                    fLBoxBPIX34ROC,                       //its logical volume
                    fPBoxBPIX34ROC,                           //its mother  volume
                    false,                       //no boolean operation
                    0); 

  // Scintillator 1

  G4Box*
  sBoxScint1 = new G4Box("Scintillator1",				     //its name
				 fScintSizeX/2,fScintSizeY/2,fScintSizeZ/2);   //its dimensions
  

  fLBoxScint1 = new G4LogicalVolume(sBoxScint1,                      //its shape
			     fMaterialScint,                   //its material
			     fMaterialScint->GetName());        //its name
  
  pos_xScint1 =  0.0*cm;
  pos_yScint1 =  0.0*cm;
  pos_zScint1 = -fWorldSizeZ/2 + 7.6*cm; // G4double pos_zScint = -10.0*cm;
  
  fPBoxScint1 =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xScint1, pos_yScint1, pos_zScint1),	       //at its position
                    "Scintillator1",                  //its name
                    fLBoxScint1,                       //its logical volume
                    fPBoxW,                           //its mother  volume
                    false,                       //no boolean operation
                    0);

  // Scintillator 2

  G4Box*
  sBoxScint2 = new G4Box("Scintillator2",				     //its name
				 fScintSizeX/2,fScintSizeY/2,fScintSizeZ/2);   //its dimensions
  

  fLBoxScint2 = new G4LogicalVolume(sBoxScint2,                      //its shape
			     fMaterialScint,                   //its material
			     fMaterialScint->GetName());        //its name
  
  pos_xScint2 =  0.0*cm;
  pos_yScint2 =  0.0*cm;
  pos_zScint2 = -fWorldSizeZ/2 + 16.7*cm; // G4double pos_zScint = -10.0*cm;

  pos_zEndScint2 = pos_zScint2 + fScintSizeZ/2;
  
  fPBoxScint2 =  new G4PVPlacement(0,			       //rotation
		    G4ThreeVector(pos_xScint2, pos_yScint2, pos_zScint2),	       //at its position
                    "Scintillator2",                  //its name
                    fLBoxScint2,                       //its logical volume
                    fPBoxW,                           //its mother  volume
                    false,                       //no boolean operation
                    0);

  // Aluminum

  G4Box*
  sBoxAl = new G4Box("Aluminum",				     //its name
				 fAlSizeX/2,fAlSizeY/2,fAlSizeZ/2);   //its dimensions
  

  fLBoxAl = new G4LogicalVolume(sBoxAl,                      //its shape
			     fMaterialAl,                   //its material
			     fMaterialAl->GetName());        //its name
  
  pos_xAl =  0.0*cm;
  pos_yAl =  0.0*cm;
  pos_zAl = -fWorldSizeZ/2 + fAlSizeZ/2;
  
  fPBoxAl =  new G4PVPlacement(0,			       //no rotation
		    G4ThreeVector(pos_xAl, pos_yAl, pos_zAl),	       //at its position
                    "Aluminum",                  //its name
                    fLBoxAl,                       //its logical volume
                    fPBoxW,                           //its mother  volume
                    false,                       //no boolean operation
                    0);

  // Visualization attributes
  G4VisAttributes* worldVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0)); //White
  worldVisAtt->SetVisibility(true);
  fLBoxW->SetVisAttributes(worldVisAtt);

  G4VisAttributes* det1VisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //Blue
  det1VisAtt->SetVisibility(true);
  fLBox1->SetVisAttributes(det1VisAtt);

  G4VisAttributes* det2VisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //Blue
  det2VisAtt->SetVisibility(true);
  fLBox2->SetVisAttributes(det2VisAtt);

  G4VisAttributes* det12VisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //White, to be set to blue
  det12VisAtt->SetVisibility(true);
  fLBoxBPIX12->SetVisAttributes(det12VisAtt);

  G4VisAttributes* det34VisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); //Blue
  det34VisAtt->SetVisibility(true);
  fLBoxBPIX34->SetVisAttributes(det34VisAtt);

  G4VisAttributes* scint1VisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); //Green
  scint1VisAtt->SetVisibility(true);
  fLBoxScint1->SetVisAttributes(scint1VisAtt);

  G4VisAttributes* scint2VisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); //Green
  scint2VisAtt->SetVisibility(true);
  fLBoxScint2->SetVisAttributes(scint2VisAtt);

  G4VisAttributes* AlVisAtt = new G4VisAttributes(G4Colour(0.75,0.75,0.75)); //Silver
  AlVisAtt->SetVisibility(true);
  fLBoxAl->SetVisAttributes(AlVisAtt);
                           
  PrintParameters();
  
  //always return the root volume
  //
  return fPBoxW;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The World is " << G4BestUnit(fWorldSizeZ,"Length")
         << " of " << fMaterialW->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterialW(const G4String& nameW)
{
  // search the material by its name
  G4Material* matW = G4Material::GetMaterial(nameW, false);

  // create the material by its name
  if(!matW) { matW = G4NistManager::Instance()->FindOrBuildMaterial(nameW); }

  if(matW != fMaterialW) {
    G4cout << "### New material " << matW->GetName() << G4endl;
    fMaterialW = matW;
    UpdateGeometry();
  }

  if(!matW) {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterialW : "
           << nameW << " not found" << G4endl;  
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterialD(const G4String& nameD)
{
  // search the material by its name
  G4Material* matD = G4Material::GetMaterial(nameD, false);

  // create the material by its name
  if(!matD) { matD = G4NistManager::Instance()->FindOrBuildMaterial(nameD); }

  if(matD != fMaterialD) {
    G4cout << "### New material " << matD->GetName() << G4endl;
    fMaterialD = matD;
    UpdateGeometry();
  }

  if(!matD) {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterialD : "
           << nameD << " not found" << G4endl;  
  } 
}

void DetectorConstruction::SetMaterialBPIX(const G4String& nameBPIX)
{
  // search the material by its name
  G4Material* matBPIX = G4Material::GetMaterial(nameBPIX, false);

  // create the material by its name
  if(!matBPIX) { matBPIX = G4NistManager::Instance()->FindOrBuildMaterial(nameBPIX); }

  if(matBPIX != fMaterialBPIX) {
    G4cout << "### New material " << matBPIX->GetName() << G4endl;
    fMaterialBPIX = matBPIX;
    UpdateGeometry();
  }

  if(!matBPIX) {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterialBPIX : "
           << nameBPIX << " not found" << G4endl;  
  } 
}

void DetectorConstruction::SetMaterialScint(const G4String& nameScint)
{
  // search the material by its name
  G4Material* matScint = G4Material::GetMaterial(nameScint, false);

  // create the material by its name
  if(!matScint) { matScint = G4NistManager::Instance()->FindOrBuildMaterial(nameScint); }

  if(matScint != fMaterialScint) {
    G4cout << "### New material " << matScint->GetName() << G4endl;
    fMaterialScint = matScint;
    UpdateGeometry();
  }

  if(!matScint) {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterialScint : "
           << nameScint << " not found" << G4endl;  
  } 
}

void DetectorConstruction::SetMaterialAl(const G4String& nameAl)
{
  // search the material by its name
  G4Material* matAl = G4Material::GetMaterial(nameAl, false);

  // create the material by its name
  if(!matAl) { matAl = G4NistManager::Instance()->FindOrBuildMaterial(nameAl); }

  if(matAl != fMaterialAl) {
    G4cout << "### New material " << matAl->GetName() << G4endl;
    fMaterialAl = matAl;
    UpdateGeometry();
  }

  if(!matAl) {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterialAl : "
           << nameAl << " not found" << G4endl;  
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeW(G4double valueW)
{
  fWorldSizeZ = valueW;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize1(G4double value1)
{
  fDet1SizeZ = value1;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize2(G4double value2)
{
  fDet2SizeZ = value2;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void DetectorConstruction::SetDist2(G4double value2)
{
  fDist2 = value2;
  UpdateGeometry();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive Detectors Absorber

  //if (!fCalorimeterSD.Get()) {
    //CalorimeterSD* calorimeterSD = new CalorimeterSD("CalorSD",this);
    //fCalorimeterSD.Put(calorimeterSD);
  //}  
  //G4SDManager::GetSDMpointer()->AddNewDetector(fCalorimeterSD.Get());
  //SetSensitiveDetector(fLogicAbsorber, fCalorimeterSD.Get());

  // Construct the field creator - this will register the field it creates

  if (!fEmFieldSetup.Get()) { 
    ElectricFieldSetup* fieldSetup = new ElectricFieldSetup();
    G4AutoDelete::Register(fieldSetup); //Kernel will delete the messenger
    fEmFieldSetup.Put(fieldSetup);
  } 
 

  fElField1z = -(fPotBackplane1-fPotStrip1)/fDet1SizeZ;
  //G4cout 
      //<< "\n Electric field inside Detector 1: " 
      //<< G4BestUnit(fElField1z, "Electric field")
      //<< G4endl;

  //2S sensor 1
  fElField1F = new G4UniformElectricField(G4ThreeVector(0.0, 0.0, fElField1z));

  fLocalEquation1F = new G4EqMagElectricField(fElField1F);

  G4int nvar1F = 8;
  fLocalStepper1F = new G4ClassicalRK4(fLocalEquation1F, nvar1F);

  G4double fMinStep1F = 0.010*mm;

  fIntgrDriver1F = new G4MagInt_Driver(fMinStep1F, fLocalStepper1F, fLocalStepper1F->GetNumberOfVariables());
  fLocalChordFinder1F = new G4ChordFinder(fIntgrDriver1F);

  fLocalFieldManager1F = new G4FieldManager();
  fLocalFieldManager1F->SetDetectorField(fElField1F);
  fLocalFieldManager1F->SetChordFinder(fLocalChordFinder1F);

  G4bool allLocal1F = true ;
  fLBox1->SetFieldManager(fLocalFieldManager1F, allLocal1F);

  //2S sensor 2
  fElField1B = new G4UniformElectricField(G4ThreeVector(0.0, 0.0, -fElField1z));

  fLocalEquation1B = new G4EqMagElectricField(fElField1B);

  G4int nvar1B = 8;
  fLocalStepper1B = new G4ClassicalRK4(fLocalEquation1B, nvar1B);

  G4double fMinStep1B = 0.010*mm;

  fIntgrDriver1B = new G4MagInt_Driver(fMinStep1B, fLocalStepper1B, fLocalStepper1B->GetNumberOfVariables());
  fLocalChordFinder1B = new G4ChordFinder(fIntgrDriver1B);

  fLocalFieldManager1B = new G4FieldManager();
  fLocalFieldManager1B->SetDetectorField(fElField1B);
  fLocalFieldManager1B->SetChordFinder(fLocalChordFinder1B);

  G4bool allLocal1B = true ;
  fLBox2->SetFieldManager(fLocalFieldManager1B, allLocal1B);


  //Pixel sensors
  fElField2z = -(fPotBackplane2-fPotStrip2)/fDet2SizeZ;
  fElField2 = new G4UniformElectricField(G4ThreeVector(0.0, 0.0, fElField2z));

  fLocalEquation2 = new G4EqMagElectricField(fElField2);

  G4int nvar2 = 8;
  fLocalStepper2 = new G4ClassicalRK4(fLocalEquation2, nvar2);

  G4double fMinStep2 = 0.010*mm;

  fIntgrDriver2 = new G4MagInt_Driver(fMinStep2, fLocalStepper2, fLocalStepper2->GetNumberOfVariables());
  fLocalChordFinder2 = new G4ChordFinder(fIntgrDriver2);

  fLocalFieldManager2 = new G4FieldManager();
  fLocalFieldManager2->SetDetectorField(fElField2);
  fLocalFieldManager2->SetChordFinder(fLocalChordFinder2);

  G4bool allLocal2 = true ;
  fLBoxBPIX12Sen->SetFieldManager(fLocalFieldManager2, allLocal2);
  fLBoxBPIX34Sen->SetFieldManager(fLocalFieldManager2, allLocal2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
