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
// $Id: DetectorConstruction.hh 66241 2012-12-13 18:34:42Z gunter $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include "G4Cache.hh"
#include "G4ElectricField.hh"
#include "G4EqMagElectricField.hh"


class MagneticField;

class G4LogicalVolume;
class G4Material;
class G4UniformMagField;
class G4UniformElectricField;
class DetectorMessenger;
class ElectricFieldSetup;
class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class G4MagInt_Driver;
class FieldMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
     virtual G4VPhysicalVolume* Construct();
     virtual void ConstructSDandField();
     
     void SetSizeW     (G4double);
     void SetSize1     (G4double);
     void SetSize2     (G4double);
     void SetDist2     (G4double);
     void SetMaterialW (const G4String&);
     void SetMaterialD (const G4String&);
     void SetMaterialBPIX (const G4String&);
     void SetMaterialScint (const G4String&);
     void SetMaterialAl (const G4String&);

     void GetFieldValue (G4double Point[4], G4double* Efield);


     void UpdateGeometry();
     
  public:
  
     const
     G4VPhysicalVolume* GetWorld()      {return fPBoxW;};
     G4VPhysicalVolume* GetDet1()       {return fPBox1;};
     G4VPhysicalVolume* GetDet2()       {return fPBox2;}; 
     G4VPhysicalVolume* GetDUT()        {return fPBoxDUT;};
     G4VPhysicalVolume* GetBPIX12()       {return fPBoxBPIX12;};
     G4VPhysicalVolume* GetBPIX12Sen()       {return fPBoxBPIX12Sen;};
     G4VPhysicalVolume* GetBPIX12ROC()       {return fPBoxBPIX12ROC;};
     G4VPhysicalVolume* GetBPIX34()       {return fPBoxBPIX34;}; 
     G4VPhysicalVolume* GetBPIX34Sen()       {return fPBoxBPIX34Sen;}; 
     G4VPhysicalVolume* GetBPIX34ROC()       {return fPBoxBPIX34ROC;}; 
     G4VPhysicalVolume* GetScint1()       {return fPBoxScint1;}; 
     G4VPhysicalVolume* GetScint2()       {return fPBoxScint2;}; 
     G4VPhysicalVolume* GetAl()       {return fPBoxAl;}; 
 
     G4LogicalVolume* GetLWorld()      {return fLBoxW;};
     G4LogicalVolume* GetLDet1()       {return fLBox1;};
     G4LogicalVolume* GetLDet2()       {return fLBox2;}; 
     G4LogicalVolume* GetLDUT()        {return fLBoxDUT;}; 
     G4LogicalVolume* GetLBPIX12()       {return fLBoxBPIX12;};
     G4LogicalVolume* GetLBPIX12Sen()       {return fLBoxBPIX12Sen;};
     G4LogicalVolume* GetLBPIX12ROC()       {return fLBoxBPIX12ROC;};
     G4LogicalVolume* GetLBPIX34()       {return fLBoxBPIX34;};  
     G4LogicalVolume* GetLBPIX34Sen()       {return fLBoxBPIX34Sen;}; 
     G4LogicalVolume* GetLBPIX34ROC()       {return fLBoxBPIX34ROC;};       
     G4LogicalVolume* GetLScint1()       {return fLBoxScint1;};  
     G4LogicalVolume* GetLScint2()       {return fLBoxScint2;};  
     G4LogicalVolume* GetLAl()       {return fLBoxAl;};                        
                    
     G4double           GetSizeW()       {return fWorldSizeZ;};
     G4double           GetSize1()       {return fDet1SizeZ;};       
     G4double           GetSize2()       {return fDet2SizeZ;};
     G4double           GetSizeBPIX()       {return fBPIXSizeZ;};   
     G4double           GetSizeBPIXSen()       {return fBPIXSenSizeZ;};           
     G4double           GetDist()        {return fDist;};

     G4double           GetStrip1Depth() {return fStrip1Depth;};
     G4double           GetStrip1Length() {return fStrip1Length;};       
     G4double           GetStrip2Depth() {return fStrip2Depth;};
     G4double           GetStrip2Length() {return fStrip2Length;};
     G4double           GetStripDist()   {return fStripDist;};
     G4double           GetStripWidth()  {return fStripWidth;};
     G4double           Getpos_EndArm1Abs()  {return pos_EndArm1Abs;};
     G4double           Getpos_BeginningArm2Abs()  {return pos_BeginningArm2Abs;};

     G4double           GetPosZAl() {return pos_zAl;};
     G4double           GetPosZScint1() {return pos_zScint1;};
     G4double           GetPosZScint2() {return pos_zScint2;};
     G4double           GetPosZBPIX12() {return pos_zBPIX12;};
     G4double           GetPosZBPIX12Sen() {return pos_zBPIX12Sen;};
     G4double           GetPosZBPIX12ROC() {return pos_zBPIX12ROC;};
     G4double           GetPosZBPIX34() {return pos_zBPIX34;};
     G4double           GetPosZBPIX34Sen() {return pos_zBPIX34Sen;};
     G4double           GetPosZBPIX34ROC() {return pos_zBPIX34ROC;};
     G4double           GetPosZ1() {return pos_z1;};
     G4double           GetPosZ2() {return pos_z2;};
     G4double           GetPosZDUT() {return pos_zDUT;};

     G4double           Getpos_zEndScint2()  {return pos_zEndScint2;};

     G4double           GetPixSensorSizeX() {return fPixSensorSizeX;};
     G4double           GetPixSensorSizeY() {return fPixSensorSizeY;};
     G4double           GetPixelPitchX() {return fPixelPitchX;};
     G4double           GetPixelDoublePitchX() {return fPixelDoublePitchX;};
     G4double           GetPixelPitchY() {return fPixelPitchY;};
     G4double           GetPixelDoublePitchY() {return fPixelDoublePitchY;};
     G4double           GetPixelDepth() {return fPixelDepth;};

     G4double           GetDUTangleX()   {return fDUTangleX;};
     G4double           GetBPIXangleX()   {return fBPIXangleX;};
     G4double           GetBPIXangleY()   {return fBPIXangleY;};
     G4double           GetScintangleX()   {return fScintangleX;};
     G4double           GetScintangleY()   {return fScintangleY;};

     G4double           GetElField1()   {return fElField1z;};
     G4double           GetElField2()   {return fElField2z;};

     G4Material*        GetMaterialW()   {return fMaterialW;};
     G4Material*        GetMaterialD()   {return fMaterialD;};
     G4Material*        GetMaterialBPIX()   {return fMaterialBPIX;};
     G4Material*        GetMaterialScint()   {return fMaterialScint;};
     G4Material*        GetMaterialAl()   {return fMaterialAl;};
     
     void               PrintParameters();
                       
  private:
  
     G4VPhysicalVolume* fPBoxW;
     G4VPhysicalVolume* fPBox1;
     G4VPhysicalVolume* fPBox2;
     G4VPhysicalVolume* fPBoxDUT;
     G4VPhysicalVolume* fPBoxBPIX12;
     G4VPhysicalVolume* fPBoxBPIX12Sen;
     G4VPhysicalVolume* fPBoxBPIX12ROC;
     G4VPhysicalVolume* fPBoxBPIX34;
     G4VPhysicalVolume* fPBoxBPIX34Sen;
     G4VPhysicalVolume* fPBoxBPIX34ROC;
     G4VPhysicalVolume* fPBoxScint1;
     G4VPhysicalVolume* fPBoxScint2;
     G4VPhysicalVolume* fPBoxAl;

     G4LogicalVolume*   fLBoxW;
     G4LogicalVolume*   fLBox1;
     G4LogicalVolume*   fLBox2;
     G4LogicalVolume*   fLBoxDUT;
     G4LogicalVolume*   fLBoxBPIX12;
     G4LogicalVolume*   fLBoxBPIX12Sen;
     G4LogicalVolume*   fLBoxBPIX12ROC;
     G4LogicalVolume*   fLBoxBPIX34;
     G4LogicalVolume*   fLBoxBPIX34Sen;
     G4LogicalVolume*   fLBoxBPIX34ROC;
     G4LogicalVolume*   fLBoxScint1;
     G4LogicalVolume*   fLBoxScint2;
     G4LogicalVolume*   fLBoxAl;
     
     G4double           zOff; //Offset in z-direction. z = 0 cm in simulation corresponds to z = 13.5 cm in real life.
     G4double           fWorldSizeX;
     G4double           fWorldSizeY;
     G4double           fWorldSizeZ;
     G4double           fDet1SizeX;
     G4double           fDet1SizeY;
     G4double           fDet1SizeZ;
     G4double           fDet2SizeX;
     G4double           fDet2SizeY;
     G4double           fDet2SizeZ;
     G4double           fBPIXSizeX;
     G4double           fBPIXSizeY;
     G4double           fBPIXSizeZ;
     G4double           fBPIXSenSizeX;
     G4double           fBPIXSenSizeY;
     G4double           fBPIXSenSizeZ;
     G4double           fBPIXROCSizeX;
     G4double           fBPIXROCSizeY;
     G4double           fBPIXROCSizeZ;
     G4double           fDUTSizeX;
     G4double           fDUTSizeY;
     G4double           fDUTSizeZ;
     G4double           fScintSizeX;
     G4double           fScintSizeY;
     G4double           fScintSizeZ;
     G4double           fAlSizeX;
     G4double           fAlSizeY;
     G4double           fAlSizeZ;
     G4double           fStrip1Depth;
     G4double           fStrip1Length;
     G4double           fStrip2Depth;
     G4double           fStrip2Length;
     G4double           fStripDist;
     G4double           fStripWidth;
     G4double           fStripPitch;
     G4double           fPixSensorSizeX;
     G4double           fPixSensorSizeY;
     G4double           fPixelPitchX;
     G4double           fPixelDoublePitchX;
     G4double           fPixelPitchY;
     G4double           fPixelDoublePitchY;
     G4double           fPixelDepth;
     G4double           fDist;
     G4double           fDist2;
     G4double           fInterArmDeltaZ;
     G4double           pos_EndArm1Abs;
     G4double           pos_BeginningArm2Abs;
     G4double           fDUTangleX;
     G4double           fBPIXangleX;
     G4double           fBPIXangleY;
     G4double           fScintangleX;
     G4double           fScintangleY;
     G4double           fElField1z;
     G4double           fElField2z;

     G4double		pos_xDUT, pos_yDUT, pos_zDUT;
     G4double		pos_x1, pos_y1, pos_z1;
     G4double		pos_x2, pos_y2, pos_z2;
     G4double		pos_xBPIX12, pos_yBPIX12, pos_zBPIX12;
     G4double		pos_xBPIX12Sen, pos_yBPIX12Sen, pos_zBPIX12Sen;
     G4double		pos_xBPIX12ROC, pos_yBPIX12ROC, pos_zBPIX12ROC;
     G4double		pos_xBPIX34, pos_yBPIX34, pos_zBPIX34;
     G4double		pos_xBPIX34Sen, pos_yBPIX34Sen, pos_zBPIX34Sen;
     G4double		pos_xBPIX34ROC, pos_yBPIX34ROC, pos_zBPIX34ROC;
     G4double		pos_xScint1, pos_yScint1, pos_zScint1;
     G4double		pos_xScint2, pos_yScint2, pos_zScint2;
     G4double		pos_xAl, pos_yAl, pos_zAl;

     G4double           fRestZ; //detZ - stripDepth

     G4double           fPotStrip1;
     G4double           fPotBackplane1;
     G4double           fPotStrip2;
     G4double           fPotBackplane2;
     G4double           pos_zEndScint2;

     G4Material*        fMaterialW;
     G4Material*        fMaterialD;      
     G4Material*        fMaterialBPIX;
     G4Material*        fMaterialScint; 
     G4Material*        fMaterialAl; 
     
     DetectorMessenger* fDetectorMessenger;
     G4Cache<ElectricFieldSetup*> fEmFieldSetup;

     G4FieldManager*         fLocalFieldManager1F;
     G4EqMagElectricField*   fLocalEquation1F;
     G4ChordFinder*          fLocalChordFinder1F;
     G4ElectricField*        fElField1F;
     G4MagIntegratorStepper* fLocalStepper1F;
     G4MagInt_Driver*        fIntgrDriver1F;

     G4FieldManager*         fLocalFieldManager1B;
     G4EqMagElectricField*   fLocalEquation1B;
     G4ChordFinder*          fLocalChordFinder1B;
     G4ElectricField*        fElField1B;
     G4MagIntegratorStepper* fLocalStepper1B;
     G4MagInt_Driver*        fIntgrDriver1B;

     G4FieldManager*         fLocalFieldManager2;
     G4EqMagElectricField*   fLocalEquation2;
     G4ChordFinder*          fLocalChordFinder2;
     G4ElectricField*        fElField2;
     G4MagIntegratorStepper* fLocalStepper2;
     G4MagInt_Driver*        fIntgrDriver2;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();
     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

