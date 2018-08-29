
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

int FindStop(double mean, double stdDev,std::vector<double>* wf,int peakSample){
    if(peakSample>2075)return peakSample+3;
    int stop=-999;
    for(int l=1;l<25;++l){
        if(fabs(((*wf)[peakSample+l]-mean))<stdDev || (*wf)[peakSample+l]-mean>0){
            stop=peakSample+l;
            break;
        }
    }
    if(stop==-999){stop=peakSample+3;}
    return stop;
}

int FindStart(double mean, double stdDev,std::vector<double>* wf,int peakSample){
    if(peakSample<25 )return peakSample-1;
    int start=-999;
    for(int l=1;l<25;++l){
        if(fabs(((*wf)[peakSample-l]-mean))<stdDev || (*wf)[peakSample-l]-mean>0 ){
            start=peakSample-l;
            break;
        }
    }
    if(start==-999){start=peakSample-2;}
    return start;
}

void FindHits(TTree* input,TTree* output,int digitN){
    int NB=3;
    int NCPMT=7;
    int NC=NCPMT+1;
    int MAXSAMPLES=2100;
    

    TString* file_name(new TString);
    unsigned int          event_number_g;
    int           computer_secIntoEpoch_g;
    long long         computer_nsIntoSec_g;
    unsigned int          digitizer_time_one_g;
     int                   digitizer_dup_one_g;
    unsigned short       digitizer_waveforms_one_g[NC][MAXSAMPLES];
    
    TBranch        *b_digitizer_waveforms;   //!
    
    input->SetBranchAddress("event_number_global",&event_number_g);
    input->SetBranchAddress("file_name",&file_name);
    input->SetBranchAddress("computer_secIntoEpoch",&computer_secIntoEpoch_g);
    input->SetBranchAddress("computer_nsIntoSec",&computer_nsIntoSec_g);
    input->SetBranchAddress("digitizer_time",&digitizer_time_one_g);
    input->SetBranchAddress("digitizer_waveforms",digitizer_waveforms_one_g);
    input->SetBranchAddress("digitizer_duplicate",&digitizer_dup_one_g);
    
    std::vector<int> pmtNum;
    std::vector<double> pmtCharge;
    std::vector<int> pmtTimeStart;
    std::vector<int> pmtTimeStop;
    std::vector<int> pmtTimePeak;
    

    unsigned short RF_ADC[2100];
    std::vector<int> RF_timeStart;
    std::vector<int> RF_timeStop;

    
    TString file_name_c;
    unsigned int          event_number_g_c;
    int           computer_secIntoEpoch_g_c;
    long long         computer_nsIntoSec_g_c;
     int                   digitizer_dup_one_g_c;
    unsigned int          digitizer_time_one_g_c;
    
    output->Branch("event_number_global",&event_number_g_c);
    output->Branch("file_name",&file_name_c);
    output->Branch("computer_secIntoEpoch",&computer_secIntoEpoch_g_c);
    output->Branch("computer_nsIntoSec",&computer_nsIntoSec_g_c);
    output->Branch("digitizer_duplicate",&digitizer_dup_one_g_c);
    
    output->Branch("digitizer_time",&digitizer_time_one_g_c);
    
    output->Branch("hit_pmt",&pmtNum);
    output->Branch("hit_charge",&pmtCharge);
    output->Branch("hit_timeStart",&pmtTimeStart);
    output->Branch("hit_timeStop",&pmtTimeStop);
    output->Branch("hit_timePeak",&pmtTimePeak);
    
    
    output->Branch("RF_waveform",&RF_ADC,"RF_ADC[2100]/i");
    output->Branch("RF_timeStart",&RF_timeStart);
    output->Branch("RF_timeStop",&RF_timeStop);

    //calibration constants
    
    double calibConst[3][7];

    calibConst[0][0] = 16.067827;
    calibConst[0][1] = 16.525127;
    calibConst[0][2] = 14.188653;
    calibConst[0][3] = 14.454693;
    calibConst[0][4] = 18.382082;
    calibConst[0][5] = 16.412492;
    calibConst[0][6] = 13.077328;
    
    calibConst[1][0] = 15.618638;
    calibConst[1][1] = 13.337854;
    calibConst[1][2] = 14.853365;
    calibConst[1][3] = 13.151699;
    calibConst[1][4] = 13.937958;
    calibConst[1][5] = 9.0;
    calibConst[1][6] = 9.3;
    
    calibConst[2][0] = 11.376641;
    calibConst[2][1] = 12.508981;
    calibConst[2][2] = 12.759682;
    calibConst[2][3] = 10.936271;
    calibConst[2][4] = 17.772873;
    calibConst[2][5] = 17.732332;
    calibConst[2][6] = 12.473879;
    
    double peakCut[3][7];
    
    peakCut[0][0]=8;
    peakCut[0][1]=9;
    peakCut[0][2]=9;
    peakCut[0][3]=7;
    peakCut[0][4]=6;
    peakCut[0][5]=6;
    peakCut[0][6]=8;
    
    peakCut[1][0]=5;
    peakCut[1][1]=5;
    peakCut[1][2]=6;
    peakCut[1][3]=5;
    peakCut[1][4]=6;
    peakCut[1][5]=5;
    peakCut[1][6]=5;
    
    peakCut[2][0]=4;
    peakCut[2][1]=5;
    peakCut[2][2]=5;
    peakCut[2][3]=5;
    peakCut[2][4]=8;
    peakCut[2][5]=8;
    peakCut[2][6]=5;
    
    int entry1=input->GetEntries();
    std::cout<<entry1<<std::endl;
    int n=1000;
    for(int i=0;i<383454;++i){
       // 383454 is end of low int
        //total number of entries is 446141
        input->GetEntry(i);
        event_number_g_c=event_number_g;
        file_name_c=(*file_name);
       computer_secIntoEpoch_g_c= computer_secIntoEpoch_g;
        computer_nsIntoSec_g_c=computer_nsIntoSec_g;
        digitizer_time_one_g_c=digitizer_time_one_g;
        
        for(int pmt=0;pmt<7;++pmt){
        double mean=0;
        for(int k=0;k<MAXSAMPLES;++k){
            mean+=(double)digitizer_waveforms_one_g[pmt][k];
        }
        mean/=MAXSAMPLES;
        
        std::vector<std::pair<int,double>> peaks;
        for(int k=1;k<MAXSAMPLES-1;++k){
            if((int)digitizer_waveforms_one_g[pmt][k]<(int)digitizer_waveforms_one_g[pmt][k+1] && (int)digitizer_waveforms_one_g[pmt][k]<(int)digitizer_waveforms_one_g[pmt][k-1] && fabs((double)digitizer_waveforms_one_g[pmt][k]-mean)>peakCut[digitN][pmt] && ((double)digitizer_waveforms_one_g[pmt][k]-mean)<0){
                if(fabs((double)(digitizer_waveforms_one_g[pmt][k]-digitizer_waveforms_one_g[pmt][k-2])/calibConst[digitN][pmt])>0.4 && digitizer_waveforms_one_g[pmt][k]<digitizer_waveforms_one_g[pmt][k-2]){
                    peaks.push_back(std::make_pair(k,digitizer_waveforms_one_g[pmt][k]));
                }
                
            }
        }
        
        std::vector<double> wf ;
        for(int k=0;k<MAXSAMPLES;++k){
            wf.push_back((double)digitizer_waveforms_one_g[pmt][k]);
        }
            
        for(std::size_t j=0;j<peaks.size();){
            
            std::vector<std::pair<int,double>> closepeaks;
            closepeaks.push_back(peaks[j]);
            int shift=0;
            for(std::size_t k=j+1;k<peaks.size();++k){
                int dt = peaks[k].first-peaks[k-1].first;
                if(dt<6){
                    closepeaks.push_back(peaks[k]);
                    shift++;
                }else{break;}
            }
            j=j+1+shift;
            
            if((int)closepeaks.size()>1){
                int dist = closepeaks[1].first-closepeaks[0].first;
                if(dist>2)dist/=2;
                int start = FindStart(mean,peakCut[digitN][pmt]/3,&wf,closepeaks[0].first);
                double local = (double)digitizer_waveforms_one_g[pmt][start];
                int stop = closepeaks[0].first+dist;
                double charge = (-1)*((double)closepeaks[0].second-local)/calibConst[digitN][pmt];
                pmtNum.push_back(pmt);
                pmtCharge.push_back(charge);
                pmtTimeStart.push_back(start);
                pmtTimeStop.push_back(stop);
                pmtTimePeak.push_back(closepeaks[0].first);

                for(std::size_t k=1;k<closepeaks.size();++k){
                    
                    int gap=0;
                    int dt = closepeaks[k].first-closepeaks[k-1].first;
                    if(dt>2){gap=dt/2;}else{gap=dt;}
                    if(k==(int)closepeaks.size()-1){
                        if(dt==2){
                            int start = closepeaks[k].first;
                            int stop = FindStop(mean,peakCut[digitN][pmt]/3,&wf,closepeaks[0].first);
                            double charge = (-1)*((double)closepeaks[k].second-digitizer_waveforms_one_g[pmt][closepeaks[k].first-1])/calibConst[digitN][pmt];
                            pmtNum.push_back(pmt);
                            pmtCharge.push_back(charge);
                            pmtTimeStart.push_back(start);
                            pmtTimeStop.push_back(stop);
                            pmtTimePeak.push_back(closepeaks[k].first);
   
                            }
                        if(dt>2){
                            double local=-9999;
                            for(int s=1;s<dt;++s){
                                local=std::max(local,(double)digitizer_waveforms_one_g[pmt][closepeaks[k].first-s]);
                            }
                            int start = closepeaks[k].first-1;
                            int stop = FindStop(mean,peakCut[digitN][pmt]/3,&wf,closepeaks[0].first);
                            double charge = (-1)*((double)closepeaks[k].second-local)/calibConst[digitN][pmt];
                            pmtNum.push_back(pmt);
                            pmtCharge.push_back(charge);
                            pmtTimeStart.push_back(start);
                            pmtTimeStop.push_back(stop);
                            pmtTimePeak.push_back(closepeaks[k].first);

                            
                        }
                    }else{
                        int gap_up=0;
                        int dt_up = closepeaks[k+1].first-closepeaks[k].first;
                        if(dt_up>2){gap_up=dt_up/2;}else{gap_up=dt_up;}
                        if(dt==2){
                            int start = closepeaks[k].first;
                            int stop = closepeaks[k].first+gap_up;
                            double charge = (-1)*((double)closepeaks[k].second-(double)digitizer_waveforms_one_g[pmt][closepeaks[k].first-1])/calibConst[digitN][pmt];
                            pmtNum.push_back(pmt);
                            pmtCharge.push_back(charge);
                            pmtTimeStart.push_back(start);
                            pmtTimeStop.push_back(stop);
                            pmtTimePeak.push_back(closepeaks[k].first);

                            }
                        if(dt>2){
                            int start = closepeaks[k].first-1;
                            double local=-9999;
                            for(int s=1;s<dt;++s){
                                local=std::max(local,(double)digitizer_waveforms_one_g[pmt][closepeaks[k].first-s]);
                            }
                            int stop = closepeaks[k].first+gap_up;
                            double charge = (-1)*((double)closepeaks[k].second-local)/calibConst[digitN][pmt];
                            pmtNum.push_back(pmt);
                            pmtCharge.push_back(charge);
                            pmtTimeStart.push_back(start);
                            pmtTimeStop.push_back(stop);
                            pmtTimePeak.push_back(closepeaks[k].first);

                            
                        }
                    }
                }
            }else{
                int start = FindStart(mean,peakCut[digitN][pmt]/3,&wf,closepeaks[0].first);
                double local = (double)digitizer_waveforms_one_g[pmt][start];
                int stop = FindStop(mean,peakCut[digitN][pmt]/3,&wf,closepeaks[0].first);
                double charge = (-1)*(closepeaks[0].second-local)/calibConst[digitN][pmt];
                pmtNum.push_back(pmt);
                pmtCharge.push_back(charge);
                pmtTimeStart.push_back(start);
                pmtTimeStop.push_back(stop);
                pmtTimePeak.push_back(closepeaks[0].first);

            }
            closepeaks.clear();
            
        }
    }

        //GET RF
       int pmt=7;
        for(int k=0;k<MAXSAMPLES;++k){
            RF_ADC[k]=digitizer_waveforms_one_g[pmt][k];
        }
        double mean =0;
        for(int k=0;k<MAXSAMPLES;++k){
            mean+=(double)digitizer_waveforms_one_g[pmt][k];
        }
        mean/=MAXSAMPLES;
        double rms=0;
        for(int k=0;k<MAXSAMPLES;++k){
            rms+=pow((double)digitizer_waveforms_one_g[pmt][k]-mean,2);
        }
        rms=sqrt(rms/MAXSAMPLES);
        int rfStart=-1;
        int rfStop=-1;
        int shift=0;
        for(int k=0;k<MAXSAMPLES;++k){
            double wf = (double)digitizer_waveforms_one_g[pmt][k]-mean;
            if(wf<0 && fabs(wf)>3*rms){
                if(shift==0){rfStart=k;rfStop=-1;}
                shift++;
                
            }
            if(shift>0 && fabs(wf)<3*rms){
                rfStop=k-1;
                int rfL = rfStop-rfStart;
                if(rfL>1){
                RF_timeStart.push_back(rfStart);
                RF_timeStop.push_back(rfStop);
                }
                rfStart=-1;
                shift=0;
            }
        }
        if((int)RF_timeStart.size()==0 || (int)RF_timeStop.size()==0){
            RF_timeStart.push_back(-1);
            RF_timeStop.push_back(-1);
        }

        
        output->Fill();
        RF_timeStart.clear();
        RF_timeStop.clear();
        pmtNum.clear();
        pmtCharge.clear();
        pmtTimeStart.clear();
        pmtTimeStop.clear();
        pmtTimePeak.clear();
        
        
    }

    delete file_name;
}

int PDSfindhits(){
   
    TFile *f = new TFile("PDS_all_LowInten.root","READ");
    
    TTree *digit_one = (TTree*)f->Get("digitizer1");
    
    TTree *digit_two = (TTree*)f->Get("digitizer2");
    
    TTree *digit_three = (TTree*)f->Get("digitizer3");
    
    TFile *f_new = new TFile("PDS_all_LowInten_withHits.root","RECREATE");
    
    TTree *digit_one_new = new TTree("digitizer1","first digitizer data");
   
    TTree *digit_two_new = new TTree("digitizer2","second digitizer data");
    
    
    TTree *digit_three_new = new TTree("digitizer3","third digitizer data");
    
    FindHits(digit_one,digit_one_new,1);
    FindHits(digit_two,digit_two_new,2);
    FindHits(digit_three,digit_three_new,3);
    
    
    
    f_new->Write();


    return 0;
}

