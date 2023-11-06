#include <fstream>
#include <iostream>

#include "Pythia.h"
#include "PowhegHooks.h"
#include "LHAFortran.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "Pythia8Plugins/HepMC2.h"
#include <cassert>
using namespace Pythia8;

extern "C" {
  extern struct {
    int hdecaymode;
    int py8tune;
    int nohad;
    int hepmc;
    int noMPI;
    int bstable;
    char rnd_cseed[4];
  } cpy8pars_;
}


class MyLHAupFortran : public LHAupFortran
{
  public:
  // User-written routine that does the event generation and fills hepeup.
  // We do nothing here; we fill the common block on the fortran side
  // before we start pythia and before we call each event.
  bool fillHepRup() {return true; }
  bool fillHepEup() {return true; }
};


MyLHAupFortran LHAinstance;

Pythia pythia;

// Interface for conversion from Pythia8::Event to HepMC one:
HepMC::IO_GenEvent* ptrHepMC;
HepMC::Pythia8ToHepMC toHepMC;
std::ofstream ascii_io;

bool loadHooks;
PowhegHooks *powhegHooks = NULL;

extern "C" {


  // F77 interface to pythia8 (pythia_init)
  void pythia_init_() {

    // disable the listing of the changed setting (too long)
    pythia.readString("Init:showChangedSettings = off");

    // switch on/off various levels of simulation
    pythia.readString("SpaceShower:QEDshowerByQ = off"); // From quarks.
    pythia.readString("SpaceShower:QEDshowerByL = off"); // From Leptons.
    pythia.readString("TimeShower:QEDshowerByQ = off"); // From quarks.
    pythia.readString("TimeShower:QEDshowerByL = off"); // From Leptons.

    // tune
    if(cpy8pars_.py8tune > 0) {
      std::stringstream ss;
      ss << cpy8pars_.py8tune;
      pythia.readString("Tune:pp = "+ss.str());
      cout << "pythia8F77: setting pythia tune "<<cpy8pars_.py8tune;
    }

    if(cpy8pars_.nohad == 1) {
      pythia.readString("HadronLevel:All = off");
    }

    if(cpy8pars_.noMPI == 1) {
      pythia.readString("PartonLevel:MPI = off");
    }

    // if(true) {
    //   std::stringstream ss;
    //   ss << "Pythia8StreamHepMC-" << cpy8pars_.rnd_cseed << ".dat";
    //   outfile = ss.str();
    //   cout << "Name of FIFO file for Rivet: " <<outfile<<"\n";
    //   cout << "pythia8F77: powheg initial seed: "<<cpy8pars_.rnd_cseed<<"\n";
    // }
    
    // uncomment to set pi0 stable (this avoid a plethora
    // of photons in the final state...)
    pythia.readString("111:mayDecay = off");

    // VETOED SHOWERS AND THE LIKE
    // veto (standard settings)
    //pythia.readString("SpaceShower:pTmaxMatch = 1");
    //pythia.readString("TimeShower:pTmaxMatch = 1");

    // powheg-tailored settings
    // -----> from main31.cc
    // Add further settings that can be set in the configuration file
    pythia.settings.addMode("POWHEG:nFinal",    2, true, false, 1, 10);
    pythia.settings.addMode("POWHEG:veto",      0, true, true,  0, 2);
    pythia.settings.addMode("POWHEG:vetoCount", 3, true, false, 0, 0);
    pythia.settings.addMode("POWHEG:pThard",    0, true, true,  0, 2);
    pythia.settings.addMode("POWHEG:pTemt",     0, true, true,  0, 2);
    pythia.settings.addMode("POWHEG:emitted",   0, true, true,  0, 3);
    pythia.settings.addMode("POWHEG:pTdef",     0, true, true,  0, 2);
    pythia.settings.addMode("POWHEG:MPIveto",   0, true, true,  0, 1);
    // Load configuration file
    pythia.readFile("main31.cmnd");

    // Higgs boson decay modes:
    switch(cpy8pars_.hdecaymode){
    case 0:
      pythia.readString("25:onMode = on");
      break;
    case 1:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 1 -1");
      break;
    case 2:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 2 -2");
      break;
    case 3:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 3 -3");
      break;
    case 4:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 4 -4");
      break;
    case 5:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 5 -5");
      break;
    case 6:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 6 -6");
      break;
    case 7:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 11 -11");
      break;
    case 8:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 13 -13");
      break;
    case 9:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 15 -15");
      break;
    case 10:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 24 -24");
      break;
    case 11:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 23 23");
      break;
    case 12:
      pythia.readString("25:onMode = off");
      pythia.readString("25:onIfMatch = 22 22");
      break;
    default:
      pythia.readString("25:isResonance = false");
      // -> switch off all decay channel
      pythia.readString("25:onMode = off");
    }
      
    if(cpy8pars_.bstable == 1){
      // Set B-Hadron to be stable
      pythia.particleData.mayDecay(511,false);
      pythia.particleData.mayDecay(521,false);
      pythia.particleData.mayDecay(531,false);
      pythia.particleData.mayDecay(541,false);
      pythia.particleData.mayDecay(5122,false);
      pythia.particleData.mayDecay(5132,false);
      pythia.particleData.mayDecay(5142,false);
      pythia.particleData.mayDecay(5232,false);
      pythia.particleData.mayDecay(5242,false);
      pythia.particleData.mayDecay(5332,false);
      pythia.particleData.mayDecay(5342,false);
      pythia.particleData.mayDecay(5412,false);
      pythia.particleData.mayDecay(5414,false);
      pythia.particleData.mayDecay(5422,false);
      pythia.particleData.mayDecay(5424,false);
      pythia.particleData.mayDecay(5432,false);
      pythia.particleData.mayDecay(5434,false);
      pythia.particleData.mayDecay(5442,false);
      pythia.particleData.mayDecay(5444,false);
      pythia.particleData.mayDecay(5512,false);
      pythia.particleData.mayDecay(5514,false);
      pythia.particleData.mayDecay(5522,false);
      pythia.particleData.mayDecay(5524,false);
      pythia.particleData.mayDecay(5532,false);
      pythia.particleData.mayDecay(5534,false);
      pythia.particleData.mayDecay(5542,false);
      pythia.particleData.mayDecay(5544,false);
      pythia.particleData.mayDecay(5544,false);
    }
    // Read in main settings
    int nEvent      = pythia.settings.mode("Main:numberOfEvents");
    int nError      = pythia.settings.mode("Main:timesAllowErrors");
    // Read in POWHEG settings
    int vetoMode    = pythia.settings.mode("POWHEG:veto");
    int MPIvetoMode = pythia.settings.mode("POWHEG:MPIveto");
    loadHooks  = (vetoMode > 0 || MPIvetoMode > 0);

    if (loadHooks) {
      // Set ISR and FSR to start at the kinematical limit
      if (vetoMode > 0) {
	pythia.readString("SpaceShower:pTmaxMatch = 2");
	pythia.readString("TimeShower:pTmaxMatch = 2");
      }
      // Set MPI to start at the kinematical limit
      if (MPIvetoMode > 0) {
	pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
      }
      powhegHooks = new PowhegHooks();
      pythia.setUserHooksPtr(powhegHooks);
    } else {
      // do this (wimpy shower) if hooks are not loaded.
      // (which can be obtained setting POWHEG:veto=0)
      cout<<"powheg hooks not loaded. Will use wimpy shower"<<endl;
      pythia.readString("SpaceShower:pTmaxMatch = 1");
      pythia.readString("TimeShower:pTmaxMatch = 1");
    }

    if(cpy8pars_.hepmc > 0) {

      std::stringstream ss;
      ss << "Pythia8StreamHepMC-" << cpy8pars_.rnd_cseed << ".dat";
      std::string tmp = ss.str();
      const char * outfile = tmp.c_str();
      cout << "Name of FIFO file for Rivet: " <<outfile<<"\n";
      cout << "pythia8F77: powheg initial seed: "<<cpy8pars_.rnd_cseed<<"\n";

      //ptrHepMC = new HepMC::IO_GenEvent(std::cout);
      ascii_io.open(outfile);
      ptrHepMC = new HepMC::IO_GenEvent(ascii_io);
      //toHepMC.set_crash_on_problem();
    }

    pythia.setLHAupPtr(&LHAinstance);
    pythia.init();

    //pythia.LHAeventList();
    //pythia.info.list();
    //pythia.process.list();
    //pythia.event.list();

    cout<< " END OF pythia_init "<<endl;

  }

  void pythia_next_(int & iret){
    //    cout<<" in pythia_next"<<endl;
    iret = pythia.next();
    //cout<<" scalup is: "<<infoPtr->scalup()<<endl;
    if(iret == 0) {
      // bad event
      if (pythia.info.atEndOfFile()) {
	cout << "Info: end of Les Houches file reached" << endl;
	iret=-1;
	return ;
      }
    } else {
      // good event
      //cout<<" EVENT SIZE IS "<<pythia.event.size()<<" "<<iret<<endl;
      //cout<<" weight pythia.info.weight() "<<pythia.info.weight()<<endl;
      //cout<<"iret = "<<iret<<endl;

      if(ptrHepMC) {
	// Construct new empty HepMC event and fill it.
	// Units will be as chosen for HepMC build, but can be changed
	// by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
	HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
	toHepMC.fill_next_event( pythia , hepmcevt );
	// FIX BELOW to allow multiple weights
	//HepMC::WeightContainer& weights = hepmcevt->weights();
	//weights.clear();
	//weights["000"] = _lastWeight.central();
	// Write out the HepMC event.
	*ptrHepMC << hepmcevt;
	delete hepmcevt;
      }
    }
  }

  void pythia_to_hepevt_(const int &nmxhep, int & nhep, int * isthep,
			 int * idhep,
			 int  (*jmohep)[2], int (*jdahep)[2],
			 double (*phep)[5], double (*vhep)[4]) {
    nhep = pythia.event.size();
    if(nhep>nmxhep) {cout << "too many particles!" ; exit(-1); }
    for (int i = 0; i < pythia.event.size(); ++i) {
      *(isthep+i) = pythia.event[i].statusHepMC();
      *(idhep+i) = pythia.event[i].id();
      // All pointers should be increased by 1, since here we follow
      // the c/c++ convention of indeces from 0 to n-1, but in fortran
      // they are from 1 to n.
      (*(jmohep+i))[0] = pythia.event[i].mother1() + 1;
      (*(jmohep+i))[1] = pythia.event[i].mother2() + 1;
      (*(jdahep+i))[0] = pythia.event[i].daughter1() + 1;
      (*(jdahep+i))[1] = pythia.event[i].daughter2() + 1;
      (*(phep+i))[0] = pythia.event[i].px();
      (*(phep+i))[1] = pythia.event[i].py();
      (*(phep+i))[2] = pythia.event[i].pz();
      (*(phep+i))[3] = pythia.event[i].e();
      (*(phep+i))[4] = pythia.event[i].m();
    }
    // override mother of very first event, set to 0
    *(jmohep)[0] = 0 ;
    *(jmohep)[1] = 0 ;

  }

  void pythia_stat_() {
    pythia.stat();
    if (powhegHooks)
      delete powhegHooks;
    if (ptrHepMC)
      delete ptrHepMC;
    if (ascii_io.is_open())
      ascii_io.close();
  }

}


