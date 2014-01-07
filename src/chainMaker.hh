//#if !defined(__CINT__) || defined(__MAKECINT__)
#ifndef CHAINMAKER_HH
#define CHAINMAKER_HH

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <TString.h>
#include <assert.h>

#include "/afs/cern.ch/user/h/habdalla/Delphes-3.0.10/classes/DelphesClasses.h"
#include "/afs/cern.ch/user/h/habdalla/Delphes-3.0.10/modules/Delphes.h"
#include "/afs/cern.ch/user/h/habdalla/Delphes-3.0.10/external/ExRootAnalysis/ExRootTreeReader.h"

using namespace std;

// this function creates adds root files to a TChain passed to it by reference. 
// To make this function handle more samples, need to add a conditional for the catalog_code.
// This function is currently optimized to handle the workflow in which samples are read from /eos/store
void chainMaker(TChain *chain, TString catalog_code, Bool_t regenerate){
  
  // create TString for where in /eos/cms/store the samples are located
  TString storeDir;
  if (catalog_code == "HHToGGBB"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4/140PileUp/HHToGGBB_14TeV";}
  else if (catalog_code == "BB-4p-0-300-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/BB-4p-0-300-v1510_14TEV";}
  else if (catalog_code == "BB-4p-1300-2100-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/BB-4p-1300-2100-v1510_14TEV";}
  else if (catalog_code == "BB-4p-2100-100000"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/BB-4p-2100-100000_14TEV";}
  else if (catalog_code == "BB-4p-300-700-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/BB-4p-300-700-v1510_14TEV";}
  else if (catalog_code == "BB-4p-700-1300-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/BB-4p-700-1300-v1510_14TEV";}
  else if (catalog_code == "BBB-4p-0-600-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/BBB-4p-0-600-v1510_14TEV";}
  else if (catalog_code == "BBB-4p-1300-100000-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/BBB-4p-1300-100000-v1510_14TEV";}
  else if (catalog_code == "BBB-4p-600-1300-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/BBB-4p-600-1300-v1510_14TEV";}
  else if (catalog_code == "ttB-4p-0-900-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/ttB-4p-0-900-v1510_14TEV";}
  else if (catalog_code == "ttB-4p-1600-2500-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4v2/140PileUp/ttB-4p-1600-2500-v1510_14TEV";}
  else if (catalog_code == "ttB-4p-2500-100000-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4/140PileUp/ttB-4p-2500-100000-v1510_14TEV";}
  else if (catalog_code == "ttB-4p-900-1600-v1510"){storeDir = "/store/group/phys_higgs/upgrade/PhaseII/Configuration4/140PileUp/ttB-4p-900-1600-v1510_14TEV";}
  
  // prepare string to pass to the bash script that will generate the text file containing the root files in the sample specified
  char buffer[200];
  string arg1, arg2;
  arg1 = storeDir;
  arg2 = catalog_code;

  // print information to console
  cout << endl << "-------------------------------------------------" << endl;
  cout << "CHAINMAKER INFORMATION BELOW: " << endl;
  cout << endl <<  "storeDir: " << storeDir << endl;
  string mystring = "/afs/cern.ch/user/h/habdalla/chainfiles/local/listRootFiles.sh " + arg1 + " " + arg2;
  // cout << endl << "mystring: "  << mystring << endl << endl;
  
  // NOTE THIS STEP REQUIRES THAT YOU HAVE "LOGGED INTO" AFS USING kinit AND aklog
  // unless the user explicitly asks the text file containing root files of sample not being regenerated,
  // call the bash script in afs that will create this file
  if (regenerate) system(mystring.c_str()); 

  // create input stream
  ifstream ifs;

  // open and assert the text file created above using listRootFiles.sh
  ifs.open("/afs/cern.ch/user/h/habdalla/chainfiles/local"+TString("/")+catalog_code+TString(".txt"));

  // print to console
  cout << endl << "trying to open: /afs/cern.ch/user/h/habdalla/chainfiles/local" << TString("/") << catalog_code << TString(".txt") << endl << endl;
  assert(ifs.is_open());
  
  string line;

  // loop through the lines of the text file created above, which file containing the file names of the sample
  while (getline(ifs, line)) {
    if (line[0] == ' ') continue;
    else {
      stringstream ss(line);
      TString filename;
      ss >> filename;

      //cout << "Adding file to chain: " << storeDir+TString("/")+filename << endl;
      
      // can add the file using the xrootd server since these samples are at CERN
      // NOTE: xrootd.unl.edu automagically redirects to CERN
      // if runnign at CERN. use eoscms.cern.ch (faster than xrootd if you're already at CERN!)
      chain->Add(TString("root://eoscms.cern.ch//eos/cms")+storeDir+TString("/")+filename);
   
    }
  }
  
  // close the input stream
  ifs.close();

  // print information to console
  cout << endl << endl << "Chain successfully created.";
  if (regenerate) cout << "from newly generated list of root files";
    else cout << " from already existing list of root files. ";
  cout << "proceeding to run" << endl << endl;
  cout << endl << "-------------------------------------------------" << endl;
}

#endif