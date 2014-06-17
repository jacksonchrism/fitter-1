#include "MyMiniSim.hh"
#include <Randomize.hh>
#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <RAT/DB.hh>
#include <TCanvas.h>
#include <TFile.h>
#include <TRandom.h>
#include <RAT/Producer.hh>
#include <RAT/ProcBlock.hh>
#include <RAT/PMTCharge.hh>
#include <RAT/ChannelEfficiency.hh>
#include <RAT/PMTTransitTime.hh>
#include <RAT/PhotonTrajectory.hh>
#include <RAT/DS/MCPhoton.hh>
#include <RAT/TrackingAction.hh>

//class G4VPhysicalVolume;
//class G4VSolid;

double fCurrentParam;
double fPhotonWavelength;

namespace RAT {

MyMiniSim::MyMiniSim(double photonWavelength, double parameterIterationValue)
: MiniSim()
{


  //fix the random number seed
//  CLHEP::HepRandom::setTheSeed(1234567);
  std::cout << "making the minisim" << std::endl;
  fCurrentParam = parameterIterationValue;
  fPhotonWavelength = photonWavelength;
  fEventCtr = 0;
  const int nchannels = 90;
  const double maxangle = 90.;
  fIncidentPhotons = new TH1D("IncidentPhotons", "IncidentPhotons", nchannels, 0., maxangle);
  fSuccessfulPhotons = new TH1D("SuccessfulPhotons", "SuccessfulPhotons", nchannels, 0., maxangle);
  fReflectedPhotons = new TH1D("ReflectedPhotons", "ReflectedPhotons", nchannels, 0., maxangle);

  fSuccessfulPhotons->SetDirectory(0);
  fIncidentPhotons->SetDirectory(0);
  //these are here for debugging
  fReflectedPhotons->SetDirectory(0);
}

MyMiniSim::~MyMiniSim() {
   delete fIncidentPhotons;
   delete fSuccessfulPhotons;
   delete fReflectedPhotons;
}

void MyMiniSim::ScaleProperly(TH1D* result)
{      
      const double normFactor = result->GetBinContent( 1 );
      const double normError = result->GetBinError( 1 );
      for( int iLoop = 1; iLoop <= result->GetNbinsX(); iLoop++ )
        {
          if( result->GetBinContent( iLoop ) == 0.0 )
            continue;
          const double binVal = result->GetBinContent( iLoop ) / normFactor;
          const double errVal = sqrt( pow( result->GetBinError( iLoop ) / result->GetBinContent( iLoop ), 2 ) + pow( normError / normFactor, 2 ) ) * binVal;

          result->SetBinContent( iLoop, binVal );
          result->SetBinError( iLoop, errVal );
        }
}


TH1D* MyMiniSim::GetSimAngResp(){

   //call this from the fitter to get minisim angular response
   fSuccessfulPhotons->Sumw2();
   fIncidentPhotons->Sumw2();
   TH1D* SimAngResponse = new TH1D("SimAngResponse", "SimAngResponse", 90, 0., 90.);
   SimAngResponse->SetDirectory(0);
   SimAngResponse->Divide(fSuccessfulPhotons,fIncidentPhotons);
   return SimAngResponse; 
}

void MyMiniSim::GeneratePrimaries(G4Event* event)
{
  double dt = 0.0;

  const int fNumPhotons = 1; //number of photons in each event
  //const int fNumPhotons = 90000; //number of photons in each event
  const double fEnergy = hbarc * twopi / (fPhotonWavelength * nm);
  G4ParticleDefinition* fOpticalPhoton = G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton");

	for(int i=0; i < fNumPhotons; i++){
        double an = (i%180)/2.;
        an = an * TMath::DegToRad();
        G4ThreeVector fPosition = G4ThreeVector(0., 1000.*sin(an), 1000.*cos(an));
        G4ThreeVector fNormal = -fPosition.unit();
        G4ThreeVector fX = fNormal.orthogonal().unit();
        G4ThreeVector fY = fNormal.cross( fX ).unit();

        const double theta = G4UniformRand() * 2.0 * pi;
        const double radius = G4UniformRand() * 300. * 300.;
        G4ThreeVector dx = ( cos( theta ) * fX + sin( theta ) * fY ) * sqrt( radius );
        dx += fPosition;

	G4PrimaryVertex* vertex = new G4PrimaryVertex(dx,dt);
        G4ThreeVector momentum = fNormal*fEnergy;
	G4PrimaryParticle* photon = new G4PrimaryParticle(fOpticalPhoton, momentum.x(), momentum.y(), momentum.z());

	//rand polarization
	double phi = (G4UniformRand()*2.0-1.0)*pi;
	G4ThreeVector e1 = fNormal.orthogonal().unit();
	G4ThreeVector e2 = fNormal.unit().cross(e1);
	G4ThreeVector pol = e1*cos(phi) + e2*sin(phi);
	photon->SetPolarization(pol.x(), pol.y(), pol.z());
	photon->SetMass(0.0);

	vertex->SetPrimary(photon);
	event->AddPrimaryVertex(vertex);

	} //loop over photons 

}



void MyMiniSim::BeginOfEventAction(const G4Event* /*event*/)
{

  std::cout << "in BeginOfEventAction" << std::endl;
  fTrackingAction->SetTrackingLevel(TrackingAction::eCondensed);

}

void MyMiniSim::EndOfEventAction(const G4Event* g4ev)
{
  //this is where we see if the photon was successful
  // (i.e generated a photoelectron in the PMT bucket it entered)
  fEventCtr++;
  std::cout << "in EndOfEventAction" << std::endl;
  TVector3 pmtDirection;
     pmtDirection.SetX(0.);
     pmtDirection.SetY(0.);
     pmtDirection.SetZ(-1.); //default is negative direction, correct for this in a sec

  // Trajectory Info
  G4TrajectoryContainer* trajectories = g4ev->GetTrajectoryContainer();
  if(trajectories != NULL) 
  {
  for( size_t iTrajectory = 0; iTrajectory < trajectories->size(); iTrajectory++ ) 
    {
      PhotonTrajectory* photonTrajectory = dynamic_cast< PhotonTrajectory* >( (*trajectories)[iTrajectory] );
      if( photonTrajectory != NULL ) // Add special MCPhoton information
        {
          std::vector<DS::MCPhoton> photons = photonTrajectory->GetMCPhotons();
          for( size_t iPhoton = 0; iPhoton < photons.size(); iPhoton++ )
            {

              TVector3 inPos = photons[iPhoton].GetInPosition();
              TVector3 inDir = photons[iPhoton].GetInDirection();
              const double angle = TMath::ACos( inDir.Dot( pmtDirection ) ) * TMath::RadToDeg();
              const double xyRadius = sqrt( inPos.X() * inPos.X() + inPos.Y() * inPos.Y() );


              if( inPos.Z() < 132.0 ) continue; //ignore photons that entered PMT at 'invalid' z pos
              if(xyRadius > 137.0) continue; //ignore photons that enter bucket outisde valid bucket radius

              fIncidentPhotons->Fill(angle);
              
              if( photons[iPhoton].GetFate() == RAT::DS::MCPhoton::eReflected) fReflectedPhotons->Fill(angle); 
              if( photons[iPhoton].GetFate() == RAT::DS::MCPhoton::ePhotoelectron) fSuccessfulPhotons->Fill(angle);
            }
        }
    }
   }//if not null trajectory

 if(true){
 if(fEventCtr%100 == 0)
  {
  TH1D* SimAngResp = GetSimAngResp();
  char buffer[70];
  sprintf(buffer, "minisim_output/event_%i_%.5f.root", fEventCtr, fCurrentParam);
  std::string str (buffer);
  TFile savef(str.c_str(), "recreate");
  MyMiniSim::ScaleProperly(SimAngResp);
  SimAngResp->Write();
  fIncidentPhotons->Write();
  fSuccessfulPhotons->Write();
  fReflectedPhotons->Write();
  savef.Close();
  }
  }

}



} // namespace RAT
