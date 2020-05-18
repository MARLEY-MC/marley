// Simple detector construction class
#pragma once

#include "G4VUserDetectorConstruction.hh"

class G4PhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();

    virtual G4VPhysicalVolume* Construct() override;

    const G4VPhysicalVolume* GetWorld() const;

  protected:

    G4VPhysicalVolume* world_ = nullptr;
};
