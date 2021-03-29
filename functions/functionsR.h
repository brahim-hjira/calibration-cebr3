/**
 * @author      : o.allegrini@ipnl.in2p3.fr
 * @file        : main
 * @created     : jeudi oct. 15, 2020 16:46:04 CEST
 * @comment     : This program analyses the data recorded in the CeBr3 detectors of the classic temporal camera and the data from the detector trigger (diamond hodoscope)
 */
// TEMPORAL CAMERA LIBRARY
#include "quadratorevt.h"

//C++ LIBRARIES
#include <iostream>
#include <filereader.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <string>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <cstring>
#include <cstdlib>
#include <cstddef>
#include <cmath>
#include <dirent.h> //class used to manipulate and traverse directories
#include <stdio.h>
#include <mutex>
#include <vector>
#include <array>
#include <algorithm> //class used for the max_element funtion

//ROOT LIBRARIES
#include "TLatex.h"
#include "TROOT.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TSpectrum.h"
#include "TF1.h"
#include "TFile.h"
#include "TApplication.h"
#include "TRint.h"
#include "TAxis.h"
#include "TAttLine.h"
#include "TTimer.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TGWindow.h"
#include "TGClient.h"
#include "TPaveText.h"
#include "TPaveLabel.h"
#include "TMath.h"
#include "TH3.h"


/** SplitFilename function
 * @brief: diviser le chemin complet du nom de fichier et de l'extension pour ajouter le dépôt de sortie dans le chemin
 * @return: le chemin complet du fichier d'analyse des résultats (sans extension)
 * @note:
 * */
std::string SplitFilename (const std::string& str, std::string repository)
{
  std::string file;
  std::string interfile;
  std::size_t rm_path = str.find_last_of("/\\");
  file = str.substr(0,rm_path);
  file = (const std::string&)file;
  rm_path = file.find_last_not_of("/\\");
  file = file.substr(0,rm_path);
  rm_path = file.find_last_of("/\\");
  file = file.substr(0,rm_path);
  file+=repository;

  rm_path = str.find_last_of("/\\");
  interfile = str.substr(rm_path+1);
  interfile = (const std::string&)interfile;
  interfile = interfile.substr(0,rm_path);
  std::size_t rm_ext = interfile.find_last_of(".");
  interfile = interfile.substr(0,rm_ext);

  file+=interfile;
  return file;
}


/** sort_in_time function
 * @brief: trier le TuileEvt d'un TuileDam à partir d'un des arguments de la structure TuileEvt
 * @return: vecteur trié de TuileEvt basé sur le numéro de frame et le timestamp pour les événements ayant le même numéro de frame
 * */
struct sort_function
{
    inline bool operator() (const TuileEvt& struct1, const TuileEvt& struct2)
    {
        return ((struct1.htimestamp) < (struct2.htimestamp)) ||
               (((struct1.htimestamp) == (struct2.htimestamp)) && ((struct1.ltimestamp&0x00FFFFFF) > (struct2.ltimestamp&0x00FFFFFF)));
    }
};

/** Global analysis function
 * @brief: Cette fonction récupère les événements qui ont un numéro de trame correspondant (coïncidence temporelle)
 * @return:réseau de vecteurs contenant : [0] Valeurs CRT [1] Données énergétiques de la tuile A [2] Données énergétiques de la tuile B [3] Température de la tuile A [4] Température de la tuile B [5] Position X de la tuile A [6] Position Y de la tuile A [7] Position Z de la tuile A [8] [9] [10] Positions de la tuile B
 * @note: Améliorations futures : Enregistrement supplémentaire des positions x, y, z après les différentes options de filtrage
 * */

std::array<std::vector<double>,11> Global_analysis (std::vector<TuileEvt> TileA, std::vector<TuileEvt> TileB){
    int it_TileA = 0, it_TileB = 0, CRT_true = 0, CRT = 0, energy_TileA =0, energy_TileB = 0, frameB = 0, temper_TileA = 0, temper_TileB = 0, posX_TileA = 0, posX_TileB = 0, posY_TileA = 0, posY_TileB = 0, posZ_TileA = 0, posZ_TileB = 0;
    bool real_value = false;
    std::array<std::vector<double>,11> record_data;
    while(((TileA.at(it_TileA+1).htimestamp) == (TileA.at(it_TileA).htimestamp)) || ((TileA.at(it_TileA).htimestamp) < (TileB.at(it_TileB).htimestamp))){
        it_TileA++;
    }
    while(((TileB.at(it_TileB).htimestamp) <= (TileA.at(it_TileA).htimestamp)) && ((it_TileA < (TileA.size()-1)) && (it_TileB < (TileB.size()-1)))){
        if ((TileB.at(it_TileB).htimestamp) == (TileA.at(it_TileA).htimestamp)){
            if(CRT==0 && real_value == false){
                frameB = it_TileB;
                CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF)));
                CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF);
                energy_TileA = TileA.at(it_TileA).photons;
                energy_TileB = TileB.at(frameB).photons;
                temper_TileA = TileA.at(it_TileA).temper;
                temper_TileB = TileB.at(frameB).temper;
                posX_TileA = TileA.at(it_TileA).mX;
                posX_TileB = TileB.at(frameB).mX;
                posY_TileA = TileA.at(it_TileA).mY;
                posY_TileB = TileB.at(frameB).mY;
                posZ_TileA = TileA.at(it_TileA).mZ;
                posZ_TileB = TileB.at(frameB).mZ;
                real_value = true;
            }
            else if(real_value == true){
                if(std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(it_TileB).ltimestamp&0x00FFFFFF))) <= CRT){
                    frameB = it_TileB;
                    CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF)));
                    CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF);
                    energy_TileA = TileA.at(it_TileA).photons;
                    energy_TileB = TileB.at(frameB).photons;
                    temper_TileA = TileA.at(it_TileA).temper;
                    temper_TileB = TileB.at(frameB).temper;
                    posX_TileA = TileA.at(it_TileA).mX;
                    posX_TileB = TileB.at(frameB).mX;
                    posY_TileA = TileA.at(it_TileA).mY;
                    posY_TileB = TileB.at(frameB).mY;
                    posZ_TileA = TileA.at(it_TileA).mZ;
                    posZ_TileB = TileB.at(frameB).mZ;
                }
            }
        }
        it_TileB++;
        while ((TileB.at(it_TileB).htimestamp) > (TileA.at(it_TileA).htimestamp) && ((it_TileA < (TileA.size()-1)) && (it_TileB < (TileB.size()-1)))){
            if (TileA.at(it_TileA).htimestamp == TileB.at(frameB).htimestamp){
                if(std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF))) < CRT){
                    CRT = std::abs((int)((TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF)));
                    CRT_true = (TileA.at(it_TileA).ltimestamp&0x00FFFFFF) - (TileB.at(frameB).ltimestamp&0x00FFFFFF);
                    energy_TileA = TileA.at(it_TileA).photons;
                    energy_TileB = TileB.at(frameB).photons;
                    temper_TileA = TileA.at(it_TileA).temper;
                    temper_TileB = TileB.at(frameB).temper;
                    posX_TileA = TileA.at(it_TileA).mX;
                    posX_TileB = TileB.at(frameB).mX;
                    posY_TileA = TileA.at(it_TileA).mY;
                    posY_TileB = TileB.at(frameB).mY;
                    posZ_TileA = TileA.at(it_TileA).mZ;
                    posZ_TileB = TileB.at(frameB).mZ;
                }
                if (((TileA.at(it_TileA+1).htimestamp) > (TileA.at(it_TileA).htimestamp)) && (it_TileA+1 <= TileA.size())){
                    record_data[0].push_back((double)CRT_true*5/256);
                    record_data[1].push_back(energy_TileA);
                    record_data[2].push_back(energy_TileB);
                    record_data[3].push_back(temper_TileA);
                    record_data[4].push_back(temper_TileB);
                    record_data[5].push_back(posX_TileA);
                    record_data[6].push_back(posX_TileB);
                    record_data[7].push_back(posY_TileA);
                    record_data[8].push_back(posY_TileB);
                    record_data[9].push_back(posZ_TileA);
                    record_data[10].push_back(posZ_TileB);
                    CRT = 0;
                    CRT_true = 0;
                    energy_TileA = 0;
                    energy_TileB = 0;
                    temper_TileA = 0;
                    temper_TileB = 0;
                    posX_TileA = 0; posY_TileA = 0; posZ_TileA =0;
                    posX_TileB = 0; posY_TileB = 0; posZ_TileB =0;
                    real_value = false;
                }
            }
            it_TileA++;
        }
    }
    return record_data;
}

/** CRT filter function
 * @brief: Cette fonction calcule le pic CRT et la ligne de base dus au bruit électrique afin de filtrer les données qui correspondent aux événements physiques..
 * @return: tableau de vecteurs contenant ([0] CRT [1] Energy TileA [2] Energy TileB) d'événements physiques sous le pic.
 * @note: La fonction Global_analysis est appelée à fournir les données ayant un frame de coïncidence dans les deux tuiles qui sont comparées
 **/
    std::array<std::vector<double>,11> CRT_filter (std::vector<TuileEvt> TileA, std::vector<TuileEvt> TileB){
    std::array<std::vector<double>,11> Raw_data;
    std::array<std::vector<double>,11> Clear_data;
    double CRT_peak;
    TH1F* Hist_Raw_data = new TH1F("","",83886080,-41943040,41943040);
    Raw_data = Global_analysis(TileA, TileB);

    //Seek for the CRT peak==Rechercher le pic de la CRT
    for(int i=0; i<Raw_data[0].size(); i++){
        Hist_Raw_data->Fill(Raw_data[0].at(i));
    }
    CRT_peak = Hist_Raw_data->GetBinCenter(Hist_Raw_data->GetMaximumBin());
    //std::cout<<CRT_peak<<std::endl;  //Use it to find the CRT peak and redefine the histogram range==Utilisez-le pour trouver le pic du TRC et redéfinir la plage de l'histogramme

    //Data filtering (window +/- 5 ns which also corresponds to +/- 256 timestamp)
    for(int i=0; i<Raw_data[0].size(); i++){
        if((Raw_data[0].at(i) >= CRT_peak-5) && (Raw_data[0].at(i) <= CRT_peak+5)){
            Clear_data[0].push_back(Raw_data[0].at(i));
            Clear_data[1].push_back(Raw_data[1].at(i));
            Clear_data[2].push_back(Raw_data[2].at(i));
            Clear_data[3].push_back(Raw_data[3].at(i));
            Clear_data[4].push_back(Raw_data[4].at(i));
            Clear_data[5].push_back(Raw_data[5].at(i));
            Clear_data[6].push_back(Raw_data[6].at(i));
            Clear_data[7].push_back(Raw_data[7].at(i));
            Clear_data[8].push_back(Raw_data[8].at(i));
            Clear_data[9].push_back(Raw_data[9].at(i));
            Clear_data[10].push_back(Raw_data[10].at(i));
        }
    }
    return Clear_data;
}

/** Nrj Temper filter function
 * @brief: Cette fonction applique un filtre énergétique et/ou des filtres de température sur les événements physiques, centrés sur le pic maximum de +/- 10% pour le filtre énergétique et de +/- 0,25 degré Celsius pour la température
 * @return: réseau de vecteurs contenant ([0] CRT [1] Tuile d'énergieA [2] Tuile d'énergieB) des événements physiques sous le pic.
 * @note: La fonction CRT_filter est appelée pour fournir les données d'entrée
 **/
std::array<std::vector<double>,11> Nrj_Temper_filter(bool nrg_filtering, bool tmp_filtering, int energy_peak_TileA, int energy_peak_TileB, int temper_peak_TileA, int temper_peak_TileB, std::vector<TuileEvt> TileA, std::vector<TuileEvt> TileB){
    std::array<std::vector<double>,11> input_data;
    std::array<std::vector<double>,11> record_data;
    input_data = CRT_filter(TileA, TileB);
    for(int i=0; i<input_data[0].size(); i++){
        ///////////////// energy filtering only ///////////////////
        if(nrg_filtering == true && tmp_filtering == false){
            if((energy_peak_TileA != -1) && (energy_peak_TileB != -1)){
                //if((input_data[1].at(i) >= 1.1*energy_peak_TileA) && (input_data[2].at(i) >= 1.1*energy_peak_TileB)){ //==> Filter to study the Compton events of 1275 keV gammas from 22Na source
                if((input_data[1].at(i) >= energy_peak_TileA*0.9) && (input_data[1].at(i) <= energy_peak_TileA*1.1) && (input_data[2].at(i) >= energy_peak_TileB*0.9) && (input_data[2].at(i) <= energy_peak_TileB*1.1)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if((energy_peak_TileA == -1) && (energy_peak_TileB != -1)){
                if((input_data[2].at(i) >= energy_peak_TileB*0.9) && (input_data[2].at(i) <= energy_peak_TileB*1.1)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if((energy_peak_TileB == -1) && (energy_peak_TileA != -1)){
                if((input_data[1].at(i) >= energy_peak_TileA*0.9) && (input_data[1].at(i) <= energy_peak_TileA*1.1)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
        }
        ////////////////// temperature filtering only ///////////////////
        else if(nrg_filtering == false && tmp_filtering == true){
            if((temper_peak_TileA != -1) && (temper_peak_TileB != -1)){
                if(((input_data[3].at(i) >= temper_peak_TileA-25) && (input_data[3].at(i) <= temper_peak_TileA+25)) && ((input_data[4].at(i) >= temper_peak_TileB-25) && (input_data[4].at(i) <= temper_peak_TileB*+25))){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if(temper_peak_TileA == -1){
                if((input_data[4].at(i) >= temper_peak_TileB-25) && (input_data[4].at(i) <= temper_peak_TileB+25)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if(temper_peak_TileB == -1){
                if((input_data[3].at(i) >= temper_peak_TileA-25) && (input_data[3].at(i) >= temper_peak_TileA+25)){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
        }
        ////////////////// energy and temperature filtering ///////////////////
        else if(nrg_filtering == true && tmp_filtering == true){
            if((energy_peak_TileA != -1) && (energy_peak_TileB != -1)){
                if(((input_data[1].at(i) >= energy_peak_TileA*0.9) && (input_data[1].at(i) <= energy_peak_TileA*1.1)) && ((input_data[2].at(i) >= energy_peak_TileB*0.9) && (input_data[2].at(i) <= energy_peak_TileB*1.1)) && ((input_data[3].at(i) >= temper_peak_TileA-25) && (input_data[3].at(i) <= temper_peak_TileA+25)) && ((input_data[4].at(i) >= temper_peak_TileB-25) && (input_data[4].at(i) <= temper_peak_TileB+25))){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if(energy_peak_TileA == -1){
                if(((input_data[2].at(i) >= energy_peak_TileB*0.9) && (input_data[2].at(i) <= energy_peak_TileB*1.1)) && ((input_data[4].at(i) >= temper_peak_TileB-25) && (input_data[4].at(i) <= temper_peak_TileB+25))){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
            else if(energy_peak_TileB == -1){
                if(((input_data[1].at(i) >= energy_peak_TileA*0.9) && (input_data[1].at(i) <= energy_peak_TileA*1.1)) && ((input_data[3].at(i) >= temper_peak_TileA-25) && (input_data[3].at(i) <= temper_peak_TileA+25))){
                    record_data[0].push_back(input_data[0].at(i));
                    record_data[1].push_back(input_data[1].at(i));
                    record_data[2].push_back(input_data[2].at(i));
                    record_data[3].push_back(input_data[3].at(i));
                    record_data[4].push_back(input_data[4].at(i));
                    record_data[5].push_back(input_data[5].at(i));
                    record_data[6].push_back(input_data[6].at(i));
                    record_data[7].push_back(input_data[7].at(i));
                    record_data[8].push_back(input_data[8].at(i));
                    record_data[9].push_back(input_data[9].at(i));
                    record_data[10].push_back(input_data[10].at(i));
                }
            }
        }

    }
    return record_data;
}

/** 1D_DepthProfile function
 * @brief: Reconstruction de la 2D avec tous les événements dans une tuile sans filtre et/ou sélection du CRT
 * @return: TH1F de la profondeur d'interaction du gamma dans le détecteur.
 * @note: Principalement utilisé pour le développement et le débogage. La reconstruction avec sélection et filtres CRT est plus appropriée.
 * @note: Cette fonction n'utilise pas les coordonnées X Y Z obtenues par vectorisation. La reconstruction 3D potentielle est l'impossible
 * */

 TH1F* DepthProfile (std::vector<double> PosX, std::vector<double> PosY, std::vector<double> PosZ, char dimension){
    float dimX;
    int nbins;
    if((dimension == ('X' | 'x')) || (dimension == ('Y' | 'y'))){
        dimX = 31;
        nbins = 8;
    }
    else if(dimension == ('Z' | 'z')){
        dimX = 50;
        nbins = 8;
    }

    TH1F* histo_1d = new TH1F("profile_1D","profile_1D", nbins, 0, dimX);
    histo_1d->SetTitle(";Position (mm);Counts");
    histo_1d->GetXaxis()->SetLabelFont(22);
    histo_1d->GetXaxis()->SetTitleFont(22);
    histo_1d->GetYaxis()->SetLabelFont(22);
    histo_1d->GetYaxis()->SetTitleFont(22);
    histo_1d->SetLineWidth(2);
    if(dimension == ('X' | 'x')){
        for(int i=0; i<PosX.size(); i++){
            histo_1d->Fill((PosX.at(i)-500)/(8000/31));
        }
    }
    else if(dimension == ('Y' | 'y')){
        for(int i=0; i<PosY.size(); i++){
            histo_1d->Fill((PosY.at(i)-500)/(8000/31));
        }
    }
    else if(dimension == ('Z' | 'z')){
        for(int i=0; i<PosZ.size(); i++){
            histo_1d->Fill((PosZ.at(i)/100));
        }
    }
    return histo_1d;
}

/** 2D map function
 * @brief:  Reconstruction de la 2D avec tous les événements dans une tuile sans filtre et/ou sélection du CRT
 * @return: TH2D[8][8] TH2D [8] [8] contenant la somme des photons par pixel pour l'ensemble de l'acquisition.
 * @note: Principalement utilisé pour le développement et le débogage. La reconstruction avec sélection et filtres CRT est plus appropriée..
 * @note: Cette fonction n'utilise pas les coordonnées X Y Z obtenues par vectorisation. La reconstruction 3D potentielle est l'impossible
 * */

 TH2D* map2D (std::vector<double> PosX, std::vector<double> PosY, int dim){
    TH2D* map2D = new TH2D("2D_Map", "2D_Map", dim, 0, 3.1, dim, 0, 3.1);
    map2D->SetTitle(";X (cm);Y (cm)");
    map2D->GetXaxis()->SetLabelFont(22);
    map2D->GetXaxis()->SetTitleFont(22);
    map2D->GetYaxis()->SetLabelFont(22);
    map2D->GetYaxis()->SetTitleFont(22);
    //uint32_t map2D [dim][dim];

    //loop on data
    for(int i =0; i< PosX.size(); i++){
        map2D->Fill((PosX.at(i)-500)/(8000/3.1), (PosY.at(i)-500)/(8000/3.1));
    }
    return map2D;
}

 /** 3D map function
  * @brief: Reconstruction de la 2D avec tous les événements dans une tuile sans filtre et/ou sélection du CRT
  * @return: TH3F [8][8][3] containing the sum of the photons per pixel for the whole acquisition.
  * @note: contenant la somme des photons par pixel pour l'ensemble de l'acquisition..
  * @note: Cette fonction n'utilise pas les coordonnées X Y Z obtenues par vectorisation. La reconstruction 3D potentielle est l'impossible
  * */

  TH3F* map3D (std::vector<double> PosX, std::vector<double> PosY, std::vector<double> PosZ, int dim){
     TH3F* map3D = new TH3F("3D_Map", "3D_Map", dim, 0, 3.1, dim, 0, 3.1, 18, 0, 5);
     map3D->SetTitle(";X (cm);Y (cm);Z (cm)");
     map3D->GetXaxis()->SetLabelFont(22);
     map3D->GetXaxis()->SetTitleFont(22);
     map3D->GetYaxis()->SetLabelFont(22);
     map3D->GetYaxis()->SetTitleFont(22);
     map3D->GetZaxis()->SetLabelFont(22);
     map3D->GetZaxis()->SetTitleFont(22);

     //loop on data
     for(int i =0; i< PosX.size(); i++){
         map3D->Fill((PosX.at(i)-500)/(8000/3.1), (PosY.at(i)-500)/(8000/3.1), PosZ.at(i)/(2000/2));
     }
     return map3D;
 }

  /** E_TileA_vs_E_TileB
   * @brief: Reconstruction de la 2D avec tous les événements dans une tuile sans filtre et/ou sélection du CRT
   * @return: TH2D[8][8] contenant la somme des photons par pixel pour l'ensemble de l'acquisition
   * @note: Principalement utilisé pour le développement et le débogage. La reconstruction avec sélection et filtres CRT est plus appropriée.
   * @note: Cette fonction n'utilise pas les coordonnées X Y Z obtenues par vectorisation. La reconstruction 3D potentielle est l'impossible
   * */

   TH2D* E_TileA_vs_E_TileB (std::vector<double> E_TileA, std::vector<double> E_TileB){
      TH2D* ETAvsETB = new TH2D("", "", 1200, 0, 12000, 1200, 0, 12000);
      ETAvsETB->SetTitle(";E_TileA (a.u);E_TileB (a.u)");
      ETAvsETB->GetXaxis()->SetLabelFont(22);
      ETAvsETB->GetXaxis()->SetTitleFont(22);
      ETAvsETB->GetYaxis()->SetLabelFont(22);
      ETAvsETB->GetYaxis()->SetTitleFont(22);
      //uint32_t map2D [dim][dim];

      //loop on data
      for(int i =0; i< E_TileA.size(); i++){
          ETAvsETB->Fill(E_TileA.at(i), E_TileB.at(i));
      }
      return ETAvsETB;
  }

 /** PrintEnergyPeaks function
  * @brief: Afficher les positions des pics d'énergie dans un tableau et l'imprimer dans la coquille
  * @return: string for pretty printing
  * */
 inline std::string PrintEnergyPeaks(double* position, int nb_peak)
 {
   std::string result;
   int width   = 11;
   char corner = '+';
   char hline  = '-';
   char vline  = '|';
   char sep    = ' ';

   auto hLineSep = [&]() {
     for(int i=0; i<nb_peak+1;i++){
       result.push_back(corner);
       result.append(width,hline);
     }
     result.push_back(corner);
     result.push_back('\n');
   };
   auto hLineNumber = [&](){
       result.push_back(vline);
       std::string zero = std::to_string(-1);
       int zeroSize = zero.size();
       result.append((width-zeroSize)/2 + (zeroSize+1)%2,sep);
       result.append(zero);
       result.append((width-zeroSize)/2,sep);

       for(int i=0;i<nb_peak;i++){
           std::string num = std::to_string(i+1);
           int numSize = num.size();
           result.push_back(vline);
           result.append((width-numSize)/2 + (numSize+1)%2,sep);
           result.append(num);
           result.append((width-numSize)/2,sep);
       }
       result.push_back(vline);
       result.push_back('\n');
   };

   auto hLinePos = [&](){
       result.push_back(vline);
       std::string none = "none";
       int noneSize = none.size();
       result.append((width-noneSize)/2 + (noneSize+1)%2,sep);
       result.append(none);
       result.append((width-noneSize)/2,sep);

       for(int i=0;i<nb_peak;i++){
           std::string pos = std::to_string((int)position[i]);
           int posSize = pos.size();
           result.push_back(vline);
           result.append((width-posSize)/2 + (posSize+1)%2,sep);
           result.append(pos);
           result.append((width-posSize)/2,sep);
       }
       result.push_back(vline);
       result.push_back('\n');
   };

   //printing
     hLineSep();
     hLineNumber();
     hLineSep();
     hLinePos();
     hLineSep();
   return result;
 }
