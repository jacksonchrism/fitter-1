#include <RAT/DB.hh>
#include <RAT/G4Stream.hh>
#include <RAT/RunManager.hh>
#include <RAT/Log.hh>
#include <RAT/DS/Run.hh>
#include <RAT/Producer.hh>
#include <G4RunManager.hh>
#include <Randomize.hh>
#include <TFitter.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include "MyMiniSim.hh"
#include <RAT/ProcBlock.hh>
#include <RAT/PhotonThinning.hh>
#include <RAT/PMTCharge.hh>
#include <RAT/ChannelEfficiency.hh>
#include <RAT/PMTTransitTime.hh>
#include <RAT/HitPMTCollection.hh>
#include <RAT/PMTVariation.hh>
#include <RAT/DU/Utility.hh>

using namespace RAT;

//static int events_per_iteration = 100;
//static int events_per_iteration = 100000;
static int events_per_iteration = 1000000;
//static int events_per_iteration = 10000000;
std::string fname;
std::string fname2;
std::string dataset_info;
TFile* SNOdata;
TH1D* SNOpmtr;
TH1D* bestfit;

void ScaleProperly(TH1D* result)
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

void writeParams(double p0, double p1)
{
    std::ofstream dbAgingFile;
    dbAgingFile.open("/data/snoplus/home/kate/newfitter/AGED_CONCENTRATOR_PARAMS.ratdb");
    //dbAgingFile.open("./AGED_CONCENTRATOR_PARAMS.ratdb");
    if(dbAgingFile.is_open()){
    dbAgingFile << "{\n";
    dbAgingFile << "name: \"AGED_CONCENTRATOR_PARAMS\",\n";
    dbAgingFile << "valid_begin: [0,0],\n";
    dbAgingFile << "valid_end: [0,0],\n";
    char buffer[35];
    std::sprintf(buffer,"p0: %fd,\n", p0);
    dbAgingFile << buffer;
    std::sprintf(buffer,"p1: %fd,\n", p1);
    dbAgingFile << buffer;
    dbAgingFile << "}";
    dbAgingFile.close();
    }


}

void makeFile()
{ 
    //create file for saving parameter values
    std::ofstream outfile;
    outfile.open(fname.c_str(), std::ios::app);
    if(outfile.is_open()){
    std::cout << "Created file for saving stuff: " << fname.c_str()  << std::endl;
    outfile.close(); 
    }

    //create file for saving parameter values w/ bin values
    std::ofstream outfile2;
    outfile2.open(fname2.c_str(), std::ios::app);
    if(outfile2.is_open()){
    std::cout << "Created file for saving more stuff: " << fname2.c_str() << std::endl;
    outfile2 << 0.0 << '\t' << 1 << '\t' << 10 << '\t' << 20 << '\t' << 30 << '\t' << 40 << '\t' << 45 << '\t' << 50 << '\t' << 60 << std::endl;
    outfile2.close(); 
    }
}


void trackbinsvsParams(double p0, int nbins, double* bins)
{
    //this writes out a text file with parameter values
    std::cout << "trackbinsvsParams " << std::endl;
    std::ofstream thingy;
    thingy.open(fname2.c_str(), std::ios::app);
    if(thingy.is_open()){
    thingy << p0 << '\t';
        for(int i=0; i < nbins; i++){
        thingy << *(bins + i) << '\t';
	}
    thingy << std::endl;
    thingy.close();
    }
}

void trackChi2vsParams(double Chi2, double p0, double p1)
{
    //this writes out a text file with parameter values
    std::ofstream thingy;
    thingy.open(fname.c_str(), std::ios::app);
    if(thingy.is_open()){
    thingy << Chi2 << '\t' << p0 << '\t' << p1 << std::endl;
    thingy.close();
    }
}

void minuit_function(int& npar, double* gin, double& result, double* par, int flag) {

    //get minisim to run with the new aging parameters before the fitter loops again
    //double fit_border = min.GetParameter(0); //my aged concentrators used to have a border,
    //but it confuses the fitter a lot
    double fit_p0 = par[0];
    double fit_p1 = par[1];
    //std::cout << "WRITE THESE VLAUES TO DB" << std::endl;
    std::cout << "current fit p0: " << fit_p0 << std::endl;
    std::cout << "current fit p1: " << fit_p1 << std::endl;

    //write previoius iteration params to DB file for this iteration
    //writeParams(fit_p0, fit_p1);
    RAT::DB::Get()->SetD("AGED_CONCENTRATOR_PARAMS", "", "p0", fit_p0);
    RAT::DB::Get()->SetD("AGED_CONCENTRATOR_PARAMS", "", "p1", fit_p1);

    //std::cout << "CHECK THAT IT WROTE CORRECTLY" << std::endl;
    
//    RAT::DB::Get()->Load("AGED_CONCENTRATOR_PARAMS.ratdb");
    double p0;
    double p1;
    //get the degredation parameters
    DBLinkPtr agedConcentratorDB = DB::Get()->GetLink( "AGED_CONCENTRATOR_PARAMS"); 
    try
      {
        p0 = agedConcentratorDB->GetD( "p0" );
        p1 = agedConcentratorDB->GetD( "p1" );
      }
    catch( RAT::DBNotFoundError &e ) 
      { 
        RAT::Log::Die( "ConcentratorAgedOpticalModel::ConcentratorAgedOpticalModel: Cannot Find AGED_CONCENTRATOR Parameters" );
      }

   std::cout << "FROM FILE: p0 = " << p0 << std::endl;
   std::cout << "FROM FILE: p1 = " << p1 << std::endl;

    std::cout << "in fitter about to make a minisim" << std::endl;
    // run the MiniSim and extract angular response 
    RAT::MyMiniSim sim(386, fit_p0, fit_p1);
    std::cout << "a minisim is made, about to beamon" << std::endl;
    sim.BeamOn(events_per_iteration);
    TH1D* SIMpmtr = sim.GetSimAngResp();
    ScaleProperly(SIMpmtr);
    bestfit = (TH1D*)SIMpmtr->Clone("best_fit");

   //save bestfit
    char buffy[80];
    char dataset_chars[5] = {dataset_info[0], dataset_info[1], dataset_info[2], dataset_info[3], dataset_info[4]};    
    sprintf(buffy, "minisim_output/test_bestfit_%s.root",dataset_chars);
    std::string str(buffy);
    TFile bestf(str.c_str(), "recreate");
    char buffy2[80];
    sprintf(buffy2, "test_bestfit_%s_%.5f_%.5f.root",dataset_chars, p0, p1);
    std::string str2(buffy2);
    bestfit->SetName(str2.c_str());
    bestfit->SetTitle(str2.c_str());
    bestfit->Write();
    SNOpmtr->Write();
    bestf.Close();


    //calculate chi2
    double chi2 = 0;
    for(int i = 0; i < 60; i++){
	double ex = SNOpmtr->GetBinContent(i+1);
	double ob = SIMpmtr->GetBinContent(i+1);
        if(ex == 0.0) break;
	double err_ob = SIMpmtr->GetBinError(i+1);
	double err_ex = SNOpmtr->GetBinError(i+1);
	double err2 = err_ob*err_ob + err_ex*err_ex;
//	double err2 = err_ob*err_ob;

	double temp = (ob - ex) * (ob - ex) / err2;
	if(!std::isnan(temp)){
	chi2 += temp; 
	}
    }

    trackChi2vsParams(chi2, p0, p1);
    std::cout << "ch2 = " << chi2 << std::endl;

    int nbinswrite = 8;
    double bins[nbinswrite];
    bins[0] = SIMpmtr->GetBinContent(1);
    bins[1] = SIMpmtr->GetBinContent(10);
    bins[2] = SIMpmtr->GetBinContent(20);
    bins[3] = SIMpmtr->GetBinContent(30);
    bins[4] = SIMpmtr->GetBinContent(40);
    bins[5] = SIMpmtr->GetBinContent(45);
    bins[6] = SIMpmtr->GetBinContent(50);
    bins[7] = SIMpmtr->GetBinContent(60);
    trackbinsvsParams(p0, nbinswrite , bins);

    delete SIMpmtr;
    
    result = chi2;
}



int main(int argc, char* argv[])
{
   if ( argc != 2 ) // argc should be 2 for correct execution
   {   // We print argv[0]: it is the program name
      std::cout<<"usage: "<< argv[0]  << " may02, apr03, oct03, sep01 \n";
      return -1; // exit
   }

   std::string dataset = argv[1];
   dataset_info = dataset;


    // start RAT logging otherwise horrible error
    RAT::Log::Init("wwfitter.log");

    //set the seed (again). Setting it here only did not set the seed in
    // minisim, so that's why it's set in both places.
    CLHEP::HepRandom::setTheSeed(1234567);
    std::cout << "seed set" << std::endl;

    // quiet down G4 output
    SetG4coutStream(G4Stream::DETAIL);
    SetG4cerrStream(G4Stream::WARN);

    // load ratdb tables
    std::cout << "Load DB..." << std::endl;
    RAT::DB::Get()->LoadDefaults();
    RAT::DB::Get()->Load("AGED_CONCENTRATOR_PARAMS.ratdb");
    RAT::DB::Get()->SetS("CONCENTRATOR", "cRAT", "model_type", "ConcOpticsAged");
    RAT::DB::Get()->SetS("DETECTOR", "","geo_file", "empty.geo");
    RAT::DB::Get()->SetS("DETECTOR","", "pmt_info_file", "singlepmt.ratdb");

    RAT::DB::Get()->SetI("PMTCALIB", "", "use_qhs_hhp", 0);
    //set starting values(?)
    RAT::DB::Get()->SetD("AGED_CONCENTRATOR_PARAMS", "", "p0", 0.3);
    RAT::DB::Get()->SetD("AGED_CONCENTRATOR_PARAMS", "", "p1", 70.0);

    //make a file name to keep params in
    //okay so the file name comes out kind of crazy/ugly.
    //I tried somethings to make it niceer but they didn't work so whatever
    time_t rawtime;
    struct tm* timeinfo;
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    //std::string temp = asctime(timeinfo);
    //char* temp2;
    //temp.copy(temp2, 19, 0);
    char buffer[80];
    strftime(buffer,80, "/data/snoplus/home/kate/newfitter/fitOutput/params_%F-%H-%M.txt",timeinfo);
    fname = buffer;

    char buffer2[80];
    strftime(buffer2,80, "/data/snoplus/home/kate/newfitter/poppick/values_%F-%H-%M.txt", timeinfo);
    fname2 = buffer2;
    //create the files
    makeFile();
    std::cout << asctime(timeinfo) << std::endl;


    // set up the geant4 run environment
    // RAT's RunManager does most of the work for us
    RAT::RunManager* run_manager = new RAT::RunManager;
    G4RunManager* g4_run_manager = G4RunManager::GetRunManager();
    g4_run_manager->Initialize();
    DU::Utility::Get()->BeginOfRun();
    PhotonThinning::BeginOfRun();
    PMTCharge::Get()->BeginOfRun(); // Must preceed ChannelEfficiency
    ChannelEfficiency::Get()->BeginOfRun();
    PMTTransitTime::Get()->BeginOfRun();
    HitPMTCollection::Get()->BeginOfRun();
    PMTVariation::Get()->BeginOfRun();

    //SNOdata = new TFile("/data/snoplus/home/kate/newfitter/high_stats_aged_results_radiuscut.root");
    //SNOpmtr = (TH1D*) SNOdata->Get("p5");
    //SNOdata = new TFile("/data/snoplus/home/kate/newfitter/event_100_0.50000.root");
    //SNOpmtr = (TH1D*) SNOdata->Get("SimAngResponse");
    //SNOdata = new TFile("/data/snoplus/home/kate/newfitter/generation_plane_size_variation.root");
    //SNOpmtr = (TH1D*) SNOdata->Get("p5_300mm");

     if(dataset.compare("may02") == 0)
     {
     SNOdata = new TFile("/data/snoplus/home/kate/angularResponse/snodata/plots/pmtAngResp_may02new_fruns_386.root");
     }
     if(dataset.compare("oct03") == 0)
     {
     SNOdata = new TFile("/data/snoplus/home/kate/angularResponse/snodata/plots/pmtAngResp_oct03_fruns_386.root");
     }
     if(dataset.compare("sep01") == 0)
     {
     SNOdata = new TFile("/data/snoplus/home/kate/angularResponse/snodata/plots/pmtAngResp_sep01old_fruns_386.root");
     }
     if(dataset.compare("apr03") == 0)
     {
     SNOdata = new TFile("/data/snoplus/home/kate/angularResponse/snodata/plots/pmtAngResp_apr03_fruns_386.root");
     }

    SNOpmtr = (TH1D*) SNOdata->Get("AngularResponse");

    //test const aging
    //SNOdata = new TFile("../checkAR/planegen_const.root");
    //SNOpmtr = (TH1D*) SNOdata->Get("angres_386nm_agedconst_0p5_300mmgen_100000000ev");

//    std::cout << "test const aging" << std::endl;
//    SNOdata = new TFile("../checkAR/planegen_const_3.root");
//    SNOpmtr = (TH1D*) SNOdata->Get("hResponse");
//    SNOpmtr->Sumw2();
//    ScaleProperly(SNOpmtr);

    std::cout << "dataset " << dataset << " with " << events_per_iteration << " events" << std::endl;

    double fp0;
    double fp1;

    //get the degredation parameters
    DBLinkPtr agedConcentratorDB = DB::Get()->GetLink( "AGED_CONCENTRATOR_PARAMS"); 
    try
      {
        fp0 = agedConcentratorDB->GetD( "p0" );
        fp1 = agedConcentratorDB->GetD( "p1" );
      }
    catch( RAT::DBNotFoundError &e ) 
      { 
        RAT::Log::Die( "ConcentratorAgedOpticalModel::ConcentratorAgedOpticalModel: Cannot Find AGED_CONCENTRATOR Parameters" );
      }
    std::cout << "Starting minimzation..." << std::endl;
    //std::cout << "NOTE: Chi2 only up to 45deg!" << std::endl;
    // minimize
    TFitter min(2);
    min.SetFCN(minuit_function);

    min.SetParameter(0, "p0", fp0, fp0/5, 0., 1.);
    min.SetParameter(1, "p1", fp1, fp1/5, 0., 130.);
 
    //min.FixParameter(0);
    //std::cout<<"START VALUE FIXED" << std::endl;
    min.FixParameter(1);
    std::cout<<"p1 FIXED" << std::endl;

    //run the thing
    double arglist[100];
    arglist[0] = 5000;
    arglist[1] = 5.0;
    //min.ExecuteCommand("MIGRAD", arglist, 2); //MIGRAD fails miserably
    min.ExecuteCommand("SIMPLEX", arglist, 2);

    //double fit_border = min.GetParameter(0);
    double fit_p0 = min.GetParameter(0);
    double fit_p1 = min.GetParameter(1);
 //   std::cout << "Best fit degredation border: " << fit_border << std::endl;
    std::cout << "Best fit p0: " << fit_p0 << std::endl;
    std::cout << "Best fit p1: " << fit_p1 << std::endl;

    SNOdata->Close();

   //save bestfit
//    char buffy[80];
//    sprintf(buffy, "minisim_output/test_bestfit_%.5f_%.5f.root", fit_p0, fit_p1);
//    std::string str(buffy);
//    TFile bestf(str.c_str(), "recreate");
//    bestfit = (TH1D*) bestfit;
//    bestfit->Write();
//    bestf.Close();



}
