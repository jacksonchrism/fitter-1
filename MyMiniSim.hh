/* My MiniSim class
*
* hopefully this is simulate what i need it to for Anged Concentrator fitting
*
* Andy: This MiniSim throws some isotropic betas with a user-specified position and
* energy, and computes the mean NHIT for the events.
*/
#ifndef __RAT__MyMiniSim__
#define __RAT__MyMiniSim__

#include <RAT/MiniSim.hh>
#include <G4ThreeVector.hh>
#include <G4Event.hh>
#include <TH1D.h>
#include <string>
#include <TVector3.h>
#include <TMath.h>

namespace RAT {

class MyMiniSim : public MiniSim
{
public:
  MyMiniSim(double photonWavelength=386, double par0=0, double par1=0);
  virtual ~MyMiniSim();

  // some function from which we'll call MiniSim::BeamOn(int nevents)
  // it's not necessary to do things this way

  // getter to extract data after running simulation
  TH1D* GetSimAngResp();
  void ScaleProperly(TH1D* result);
  double GetMeanNhits();

  int fEventCtr;

  // must override GeneratePrimaries
  virtual void GeneratePrimaries(G4Event* event);

  // overriding other methods is optional
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);
  

protected:

  // store the simulation output in class members
  TH1D* fIncidentPhotons;
  TH1D* fSuccessfulPhotons;
  TH1D* fReflectedPhotons;


};

} // namespace RAT

#endif
