
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <array>
#include <utility>
#include <sstream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <valarray>


#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TString.h"


#include <TMath.h>
#include <TNtuple.h>
#include <TFile.h>
#include <Rtypes.h>

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"

int PDSrestructure(){
    
    int NB=3;
    int NCPMT=7;
    int NC=NCPMT+1;
    int MAXSAMPLES=2100;
    
    TString file_name;
    
    unsigned int          event_number_g;
    int           computer_secIntoEpoch_g;
    long long         computer_nsIntoSec_g;
    unsigned int          gps_nsIntoSec_g;
    unsigned int          gps_secIntoDay_g;
    unsigned short        gps_daysIntoYear_g;
    unsigned short        gps_Year_g;
    unsigned short        gps_ctrlFlag_g;
    unsigned int          digitizer_size_one_g;
    unsigned int          digitizer_size_two_g;
    unsigned int          digitizer_size_three_g;
    unsigned int         digitizer_chMask_one_g[NC];
    unsigned int          digitizer_chMask_two_g[NC];
    unsigned int         digitizer_chMask_three_g[NC];
    unsigned int          digitizer_evNum_one_g;
    unsigned int          digitizer_evNum_two_g;
    unsigned int          digitizer_evNum_three_g;
    unsigned int          digitizer_time_one_g;
    unsigned int         digitizer_time_two_g;
    unsigned int         digitizer_time_three_g;
    unsigned short       digitizer_waveforms_one_g[NC][MAXSAMPLES];
    unsigned short        digitizer_waveforms_two_g[NC][MAXSAMPLES];
    unsigned short         digitizer_waveforms_three_g[NC][MAXSAMPLES];
    unsigned int          nDigitizers_g;
    unsigned int          nChannels_g;
    unsigned int          nSamples_g;
    unsigned int          nData_g;
    
    TFile *f = new TFile("PDS_all_LowInten.root","RECREATE");
    
    TTree *digit_one = new TTree("digitizer1","first digitizer data");
    digit_one->Branch("event_number_global",&event_number_g);
    digit_one->Branch("file_name",&file_name);
    digit_one->Branch("computer_secIntoEpoch",&computer_secIntoEpoch_g);
    digit_one->Branch("computer_nsIntoSec",&computer_nsIntoSec_g);
    digit_one->Branch("gps_nsIntoSec",&gps_nsIntoSec_g);
    digit_one->Branch("gps_secIntoDay",&gps_secIntoDay_g);
    digit_one->Branch("gps_daysIntoYear",&gps_daysIntoYear_g);
    digit_one->Branch("gps_Year",&gps_Year_g);
    digit_one->Branch("gps_ctrlFlag",&gps_ctrlFlag_g);
    digit_one->Branch("digitizer_size",&digitizer_size_one_g);
    digit_one->Branch("digitizer_chMask",&digitizer_chMask_one_g,"digitizer_chMask_one_g[8]/i");
    digit_one->Branch("digitizer_evNum",&digitizer_evNum_one_g);
    digit_one->Branch("digitizer_time",&digitizer_time_one_g);
    
    digit_one->Branch("digitizer_waveforms",&digitizer_waveforms_one_g,"digitizer_waveforms_one_g[8][2100]/s");
    
    digit_one->Branch("nDigitizers",&nDigitizers_g);
    digit_one->Branch("nChannels",&nChannels_g);
    digit_one->Branch("nSamples",&nSamples_g);
    digit_one->Branch("nData",&nData_g);
    
    TTree *digit_two = new TTree("digitizer2","second digitizer data");
    digit_two->Branch("event_number_global",&event_number_g);
    digit_two->Branch("file_name",&file_name);
    digit_two->Branch("computer_secIntoEpoch",&computer_secIntoEpoch_g);
    digit_two->Branch("computer_nsIntoSec",&computer_nsIntoSec_g);
    digit_two->Branch("digitizer_size",&digitizer_size_two_g);
    digit_two->Branch("digitizer_chMask",&digitizer_chMask_two_g,"digitizer_chMask_two_g[8]/i");
    digit_two->Branch("digitizer_evNum",&digitizer_evNum_two_g);
    digit_two->Branch("digitizer_time",&digitizer_time_two_g);
   digit_two->Branch("digitizer_waveforms",&digitizer_waveforms_two_g,"digitizer_waveforms_two_g[8][2100]/s");
    digit_two->Branch("nDigitizers",&nDigitizers_g);
    digit_two->Branch("nChannels",&nChannels_g);
    digit_two->Branch("nSamples",&nSamples_g);
    digit_two->Branch("nData",&nData_g);
    TTree *digit_three = new TTree("digitizer3","third digitizer data");
    digit_three->Branch("event_number_global",&event_number_g);
    digit_three->Branch("file_name",&file_name);
    digit_three->Branch("computer_secIntoEpoch",&computer_secIntoEpoch_g);
    digit_three->Branch("computer_nsIntoSec",&computer_nsIntoSec_g);
    digit_three->Branch("digitizer_size",&digitizer_size_three_g);
    digit_three->Branch("digitizer_chMask",&digitizer_chMask_three_g,"digitizer_chMask_three_g[8]/i");
    digit_three->Branch("digitizer_evNum",&digitizer_evNum_three_g);
    digit_three->Branch("digitizer_time",&digitizer_time_three_g);
   digit_three->Branch("digitizer_waveforms",&digitizer_waveforms_three_g,"digitizer_waveforms_three_g[8][2100]/s");
    digit_three->Branch("nDigitizers",&nDigitizers_g);
    digit_three->Branch("nChannels",&nChannels_g);
    digit_three->Branch("nSamples",&nSamples_g);
    digit_three->Branch("nData",&nData_g);
    
    const char *dirname = "/Users/sergey/Desktop/PDS/PDS_beamtime_lowintensity_runs/";
    const char *ext = ".root";
    int totalEnrty=0;
    
    TString pwd(gSystem->pwd());
    TSystemDirectory dir(dirname, dirname);
    TList *files = dir.GetListOfFiles();
    files->Sort();
    if (files) {
        TSystemFile* file;
        TString fname;
        TIter next(files);
        gSystem->cd(pwd.Data());
        while ( (file = (TSystemFile*)next()) ) {
            fname = file->GetName();
            
            if (!file->IsDirectory() && fname.EndsWith(ext)) {
                std::cout<<fname<<std::endl;
                TFile* Tfile = new TFile(dirname+fname,"READ");
                TTree* pmtTree = (TTree*)Tfile->Get("pmt_tree");
                unsigned int          event_number;
                 int           computer_secIntoEpoch;
                long long        computer_nsIntoSec;
               unsigned int          gps_nsIntoSec;
               unsigned int          gps_secIntoDay;
               unsigned short        gps_daysIntoYear;
               unsigned short        gps_Year;
               unsigned short        gps_ctrlFlag;
               unsigned int          digitizer_size[NB];
               unsigned int          digitizer_chMask[NB][NC];
               unsigned int          digitizer_evNum[NB];
               unsigned int          digitizer_time[NB];
               unsigned short        digitizer_waveforms[NB][NC][MAXSAMPLES];

               unsigned int          nDigitizers;
               unsigned int          nChannels;
               unsigned int          nSamples;
               unsigned int          nData;

                TBranch        *b_digitizer_waveforms;   //!

                
                pmtTree->SetBranchAddress("event_number", &event_number, &b_event_number);
                pmtTree->SetBranchAddress("computer_secIntoEpoch", &computer_secIntoEpoch, &b_computer_secIntoEpoch);
                pmtTree->SetBranchAddress("computer_nsIntoSec", &computer_nsIntoSec, &b_computer_nsIntoSec);
                pmtTree->SetBranchAddress("gps_nsIntoSec", &gps_nsIntoSec, &b_gps_nsIntoSec);
                pmtTree->SetBranchAddress("gps_secIntoDay", &gps_secIntoDay, &b_gps_secIntoDay);
                pmtTree->SetBranchAddress("gps_daysIntoYear", &gps_daysIntoYear, &b_gps_daysIntoYear);
                pmtTree->SetBranchAddress("gps_Year", &gps_Year, &b_gps_Year);
                pmtTree->SetBranchAddress("gps_ctrlFlag", &gps_ctrlFlag, &b_gps_ctrlFlag);
                pmtTree->SetBranchAddress("digitizer_size", digitizer_size, &b_digitizer_size);
                pmtTree->SetBranchAddress("digitizer_chMask", digitizer_chMask, &b_digitizer_chMask);
                pmtTree->SetBranchAddress("digitizer_evNum", digitizer_evNum, &b_digitizer_evNum);
                pmtTree->SetBranchAddress("digitizer_time", digitizer_time, &b_digitizer_time);
                pmtTree->SetBranchAddress("digitizer_waveforms", digitizer_waveforms, &b_digitizer_waveforms);
                pmtTree->SetBranchAddress("nDigitizers", &nDigitizers, &b_nDigitizers);
                pmtTree->SetBranchAddress("nChannels", &nChannels, &b_nChannels);
                pmtTree->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
                pmtTree->SetBranchAddress("nData", &nData, &b_nData);
                int nentries = pmtTree->GetEntries();
                totalEnrty+=nentries;
                std::cout<<nentries<<std::endl;
                for(int i=0;i<nentries;++i){
                    pmtTree->GetEntry(i);
                    file_name = fname;
                    event_number_g=event_number;
                
                    computer_secIntoEpoch_g=computer_secIntoEpoch;
                    computer_nsIntoSec_g=computer_nsIntoSec        ;
                    
                    gps_nsIntoSec_g= gps_nsIntoSec;
                    gps_secIntoDay_g =         gps_secIntoDay;
                    gps_daysIntoYear_g =        gps_daysIntoYear;
                    gps_Year_g=        gps_Year;
                    gps_ctrlFlag_g =        gps_ctrlFlag;
                    
                    digitizer_size_one_g   =   digitizer_size[0];
                    digitizer_size_two_g   =   digitizer_size[1];
                    digitizer_size_three_g   =   digitizer_size[2];
                    
                    for(int j=0;j<NC;++j){
                        digitizer_chMask_one_g[j]=digitizer_chMask[0][j];
                        digitizer_chMask_two_g[j]=digitizer_chMask[1][j];
                        digitizer_chMask_three_g[j]=digitizer_chMask[2][j];
                        for(int k=0;k<MAXSAMPLES;++k){
                            digitizer_waveforms_one_g[j][k]=digitizer_waveforms[0][j][k];
                           digitizer_waveforms_two_g[j][k]=digitizer_waveforms[1][j][k];
                           digitizer_waveforms_three_g[j][k]=digitizer_waveforms[2][j][k];
                        }
                    }
                    
                    digitizer_evNum_one_g     =    digitizer_evNum[0];
                    digitizer_evNum_two_g     =    digitizer_evNum[1];
                    digitizer_evNum_three_g     =    digitizer_evNum[2];
                    digitizer_time_one_g    =     digitizer_time[0];
                    digitizer_time_two_g    =    digitizer_time[1];
                    digitizer_time_three_g    =     digitizer_time[2];
                    
                    nDigitizers_g =          nDigitizers;
                    nDigitizers_g =          nChannels;
                    nDigitizers_g =          nSamples;
                    nDigitizers_g =          nData;
                    digit_one->Fill();
                    digit_two->Fill();
                    digit_three->Fill();
                   /*     event_number_g=0;
                              computer_secIntoEpoch_g=0;
                    computer_nsIntoSec_g=0;
                  gps_nsIntoSec_g=0;
                    gps_secIntoDay_g=0;
                    gps_daysIntoYear_g=0;
                    gps_Year_g=0;
                    gps_ctrlFlag_g=0;
                    digitizer_size_one_g=0;
                    digitizer_size_two_g=0;
                    digitizer_size_three_g=0;
                    std::fill_n(digitizer_chMask_one_g,NC,-999);
                  std::fill_n(digitizer_chMask_two_g,NC,-999);
                  std::fill_n(digitizer_chMask_three_g,NC,-999);
                    digitizer_evNum_one_g=0;
                              digitizer_evNum_two_g=0;
                              digitizer_evNum_three_g=0;
                              digitizer_time_one_g=0;
                             digitizer_time_two_g=0;
                             digitizer_time_three_g=0;
                    
                    for(int j=0;j<NC;++j){
                        for(int k=0;k<MAXSAMPLES;++k){
                            digitizer_waveforms_one_g[j][k]=-999;
                            digitizer_waveforms_two_g[j][k]=-999;
                            digitizer_waveforms_three_g[j][k]=-999;
                        }
                    }
                    
                    
                              nDigitizers_g=0;
                              nChannels_g=0;
                              nSamples_g=0;
                              nData_g=0;*/
                }
                delete pmtTree;
                delete Tfile;
            }
        }
    }
    delete files;
    f->Write();
    f->Close();
    std::cout<<"allentries="<<totalEnrty<<std::endl;
    return 0;
}

