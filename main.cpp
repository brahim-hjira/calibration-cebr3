/**
 * @author      : b.hjira@ip2i.in2p3.fr
 * @file        : main
 * @created     : Monday, February 15th, 2021 16:46:04 CEST
 * @comment     : Ce programme analyse les données enregistrées dans les détecteurs CeBr3 de la caméra temporelle classique et sortie une carte de correction
 */
// TEMPORAL CAMERA LIBRARY
#include "quadratorevt.h"
#include <filereader.h>

//FUNCTIONS
#include <functions.h>

//C++ LIBRARIES
#include <iostream>
#include <unistd.h>
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
#include <dirent.h> //classe utilisée pour manipuler et parcourir les répertoires
#include <stdio.h>
#include <mutex>
#include <vector>
#include <array>
#include <algorithm> //classe utilisée pour la fonction max_element
#include <bits/stdc++.h>
#include <chrono>

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
#include "TSystem.h"
#include "TStyle.h"
#include "TGWindow.h"
#include "TGClient.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TH3.h"
#include <TGraph2D.h>
#include <TRandom.h>
using namespace std;

int main(int argc, char **argv)
{
    /////////////////// Starting ROOT in the c++ program /////////////////
    TApplication *theApp = new TApplication("App", &argc, argv);
    if (gROOT->IsBatch())
    {
        fprintf(stderr, "%s: cannot run in batch mode\n", argv[0]);
        return 1;
    }
    ///////////////////////////save the outputs////////////////////////

    string nomFichier("/home/hjira/calibration/output/output_root/calibration.root");
    TFile monFlux(nomFichier.c_str(), "RECREATE");
    ///////////////////////////variable declaration////////////////////
    int n ;// use in for loop==> to choose all acquisitions or one
    int choice_tile = 2;
    int m;
    int number_pixels_axis;
    int peak_found = 0;
    double *peaks_position;
    int p; //position of canal

    

    cout << "Enter the number of pixels per axis" << endl;
    cin >> number_pixels_axis;
    int number_pixels = pow(number_pixels_axis, 2);
    cout << "*** Number_pixels=***" << number_pixels << endl;

    TH1F *Photons_Spectrum_Tile2 = new TH1F("Photons_Spectrum_Tile2", "Photons_Spectrum_Tile2", 1000, 0, 20000);
    TSpectrum *spectrum = new TSpectrum();

    FileReader reader;
    FileReader reader2;
    TuileDam readTile;
    std::vector<TuileEvt> Tile;
    std::vector<TuileImg> Tile_Img;
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++print_number_events_pixel_acquisition++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    string filename;
    std::cout<<"Enter1 <path/filename.type>"<<std::endl;
    cin>>filename;
    QString Qfilename = QString::fromStdString(filename); //Filereader variables uses QString input ==> This line converts string to QString

    if (reader.openReadFile(Qfilename) == FileReader::OPENNING_ERROR)
    {
    cout << "Recherche du fichier: " << filename << " pas de fichier" << std::endl;
    return -1;
    }

    while (reader.getNextEvent(&readTile) == FileReader::OK)
    {

    if(readTile.tuile.tuile == choice_tile)
    {
    Tile.push_back(readTile.tuile);
    Tile_Img.push_back(readTile.image);
    }
    else
    {
    cout << "ERREUR D'IDENTIFICATION DE LA TUILE: " << readTile.tuile.tuile << std::endl;
    return 0;
    }
    }
    cout << "Number of events in tile " << choice_tile << " : " << Tile.size() << std::endl;

    for (int i = 0; i < Tile.size(); i++)
    {
    Photons_Spectrum_Tile2->Fill(Tile.at(i).photons);
    }

    peak_found = spectrum->Search(Photons_Spectrum_Tile2,10,"",0.05); //Search the peaks in the spectrum and displays them in the histogram
    peaks_position = spectrum->GetPositionX();
    cout<<PrintEnergyPeaks(peaks_position,peak_found)<<endl;
    cout<<"Choose the position of the energy filter window . Press 0 if you don't want apply filter"<<endl;
    cin>>p;
    print_photons_Spectrum_filter(Tile);
    print_number_htimestamp(Tile);
    print_mean_temperature(Tile);
    //--------------------------------------------------------------------------------------
    gStyle-> SetOptStat(000);
    TCanvas* Map_tuile_event= new TCanvas("Map_tuile_event","Map_tuile_event", 1400, 1000);
    TCanvas* Global=new TCanvas("","", 1600, 1000);
    TH2D* Map2D_14_die = new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));
    TH2D* Map2D_16_die = new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));
    //
    TH1D* Map1D_14_die_x = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));
    TH1D* Map1D_14_die_y = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));
    //
    TH1D* Map1D_16_die_x = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));
    TH1D* Map1D_16_die_y = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));
    float x_min=500;
    float x_max=8500;
    float y_min=0.;
    float y_max=31.;
    float a=(y_max-y_min)/(8000.);
    float b=y_min-a*x_min;
    float x=0;
    float y=0;
    float h=(x_max-x_min)/(number_pixels_axis);

for(int k=0;k<Tile.size();k++)//TileA.size()//13627070, 2e6
{				     //100000//463934

   if(p==0)
    {
	x=Tile.at(k).mX;
	y=Tile.at(k).mY;
	x=a*x+b;
	y=a*y+b;
	Map2D_14_die->Fill(x,y) ;
	Map1D_14_die_x->Fill(x) ;
	Map1D_14_die_y->Fill(y) ;

    }else{

       if ((Tile.at(k).photons >= p*(0.95))&&(Tile.at(k).photons <= p*1.05))// if ((k>=3850*0.9)&&(k<=3850*1.1))
       {                                       //p*(0.95))&&(TileA.at(k).photons <= p*1.05))
						      //5250)&&(TileA.at(k).photons <= 6775))

		x=Tile.at(k).mX;
		y=Tile.at(k).mY;
		x=a*x+b;
		y=a*y+b;
		Map2D_14_die->Fill(x,y) ;
		//Map1D_14_die_x->Fill(x) ;
		Map1D_14_die_y->Fill(y) ;
		Map1D_14_die_x=Map2D_14_die->ProjectionX();
       }
    }
    //--------------------
   if (Tile.at(k).nbDie>=16)
   {
       if(p==0)
	{
	   x=Tile.at(k).mX;
	   y=Tile.at(k).mY;
	   x=a*x+b;
	   y=a*y+b;
	   Map2D_16_die->Fill(x,y) ;
	   Map1D_16_die_x->Fill(x) ;
	   Map1D_16_die_y->Fill(y) ;

	}else{

	   if ((Tile.at(k).photons >= p*(0.95))&&(Tile.at(k).photons <= p*1.05))// if ((k>=3850*0.9)&&(k<=3850*1.1))
	   {                                       //p*(0.95))&&(TileA.at(k).photons <= p*1.05))
							  //5250)&&(TileA.at(k).photons <= 6775))

		    x=Tile.at(k).mX;
		    y=Tile.at(k).mY;
		    x=a*x+b;
		    y=a*y+b;
		    Map2D_16_die->Fill(x,y) ;
		    Map1D_16_die_x->Fill(x) ;
		    Map1D_16_die_y->Fill(y) ;
	   }
   }
   }
}
//------------------
Map2D_14_die->Draw("COLZ");
Map2D_14_die->SetTitle(";X [mm];Y [mm]");
Map2D_14_die->GetXaxis()->SetLabelFont(22);
Map2D_14_die->GetXaxis()->SetTitleFont(22);
Map2D_14_die->GetYaxis()->SetLabelFont(22);
Map2D_14_die->GetYaxis()->SetTitleFont(22);
Map2D_14_die->SetLineWidth(3);

Map1D_14_die_x->Draw("hist");
Map1D_14_die_x->SetTitle(";X [mm]; Events");
Map1D_14_die_x->GetXaxis()->SetLabelFont(22);
Map1D_14_die_x->GetXaxis()->SetTitleFont(22);
Map1D_14_die_x->GetYaxis()->SetLabelFont(22);
Map1D_14_die_x->GetYaxis()->SetTitleFont(22);
Map1D_14_die_x->SetLineWidth(4);
Map1D_14_die_x->SetLineColor(kBlue);

Map1D_14_die_y->Draw("hist");
Map1D_14_die_y->SetTitle(";Y [mm]; Events");
Map1D_14_die_y->GetXaxis()->SetLabelFont(22);
Map1D_14_die_y->GetXaxis()->SetTitleFont(22);
Map1D_14_die_y->GetYaxis()->SetLabelFont(22);
Map1D_14_die_y->GetYaxis()->SetTitleFont(22);
Map1D_14_die_y->SetLineWidth(4);
Map1D_14_die_y->SetLineColor(kBlue);


Map2D_14_die->Write("Map2D_14_die");
Map1D_14_die_x->Write("Map1D_14_die_x");
Map1D_14_die_y->Write("Map1D_14_die_y");
//------------------
Map2D_16_die->Draw("COLZ");
Map2D_16_die->SetTitle(";X [mm];Y [mm]");
Map2D_16_die->GetXaxis()->SetLabelFont(22);
Map2D_16_die->GetXaxis()->SetTitleFont(22);
Map2D_16_die->GetYaxis()->SetLabelFont(22);
Map2D_16_die->GetYaxis()->SetTitleFont(22);
Map2D_16_die->SetLineWidth(3);

Map1D_16_die_x->Draw("hist");
Map1D_16_die_x->SetTitle(";X [mm]; Events");
Map1D_16_die_x->GetXaxis()->SetLabelFont(22);
Map1D_16_die_x->GetXaxis()->SetTitleFont(22);
Map1D_16_die_x->GetYaxis()->SetLabelFont(22);
Map1D_16_die_x->GetYaxis()->SetTitleFont(22);
Map1D_16_die_x->SetLineWidth(4);
Map1D_16_die_x->SetLineColor(kBlue);

Map1D_16_die_y->Draw("hist");
Map1D_16_die_y->SetTitle(";Y [mm]; Events");
Map1D_16_die_y->GetXaxis()->SetLabelFont(22);
Map1D_16_die_y->GetXaxis()->SetTitleFont(22);
Map1D_16_die_y->GetYaxis()->SetLabelFont(22);
Map1D_16_die_y->GetYaxis()->SetTitleFont(22);
Map1D_16_die_y->SetLineWidth(4);
Map1D_16_die_y->SetLineColor(kBlue);


Map2D_16_die->Write("Map2D_16_die");
Map1D_16_die_x->Write("Map1D_16_die_x");
Map1D_16_die_y->Write("Map1D_16_die_y");

//Map11D16->GetYaxis()->SetTitle("N#circ events");
//Map11D16->GetXaxis()->SetTitle("N#circ pixels");



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++print_number_events_pixel_CDC++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
string filename2;
cout<<"Enter a number > 0 to open a correction map or 0 to quit program "<<endl;
cin>>m;
if (m==0)
{ cout<<"end program"<<endl;
}else{
    //-----------------------------
    std::cout<<"Enter2 <path/filename2.type>"<<std::endl;
    cin>>filename2;
    QString Qfilename2 = QString::fromStdString(filename2); //Filereader variables uses QString input ==> This line converts string to QString

    if (reader2.openReadFile(Qfilename2) == FileReader::OPENNING_ERROR)
    {
    cout << "Recherche du fichier: " << filename2 << " pas de fichier" << std::endl;
    return -1;
    }

    while (reader2.getNextEvent(&readTile) == FileReader::OK)
    {

    if(readTile.tuile.tuile == choice_tile)
    {
    Tile.push_back(readTile.tuile);
    Tile_Img.push_back(readTile.image);
    }
    else
    {
    cout << "ERREUR D'IDENTIFICATION DE LA TUILE: " << readTile.tuile.tuile << std::endl;
    return 0;
    }
    }
    cout << "Number of events in tile " << choice_tile << " : " << Tile.size() << std::endl;

    for (int i = 0; i < Tile.size(); i++)
    {
    Photons_Spectrum_Tile2->Fill(Tile.at(i).photons);
    }

    peak_found = spectrum->Search(Photons_Spectrum_Tile2,10,"",0.05); //Search the peaks in the spectrum and displays them in the histogram
    peaks_position = spectrum->GetPositionX();
    cout<<PrintEnergyPeaks(peaks_position,peak_found)<<endl;
    cout<<"Choose the position of the energy filter window . Press 0 if you don't want apply filter"<<endl;
    cin>>p;
    //---------------------------------
    TH2D* correction_map2D_14_die = new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));
    TH2D* correction_map2D_16_die = new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));

    float V_min_14;
    float V_min_16;
    float V_max_14;
    float V_max_16;

for(int k=0;k<Tile.size();k++)//TileA.size()//13627070, 2e6
{				     //100000//463934

   if(p==0)
    {
	x=Tile.at(k).mX;
	y=Tile.at(k).mY;
	x=a*x+b;
	y=a*y+b;
	correction_map2D_14_die->Fill(x,y) ;
    }else{

       if ((Tile.at(k).photons >= p*(0.95))&&(Tile.at(k).photons <= p*1.05))// if ((k>=3850*0.9)&&(k<=3850*1.1))
       {                                       //p*(0.95))&&(TileA.at(k).photons <= p*1.05))
						      //5250)&&(TileA.at(k).photons <= 6775))
	   x=Tile.at(k).mX;
	   y=Tile.at(k).mY;
	   x=a*x+b;
	   y=a*y+b;
	   correction_map2D_14_die->Fill(x,y) ;
       }
    }
//----------------------
   if (Tile.at(k).nbDie>=16)
   {
       if(p==0)
	{
	   x=Tile.at(k).mX;
	   y=Tile.at(k).mY;
	   x=a*x+b;
	   y=a*y+b;
	   correction_map2D_16_die->Fill(x,y) ;
	}else{

	   if ((Tile.at(k).photons >= p*(0.95))&&(Tile.at(k).photons <= p*1.05))// if ((k>=3850*0.9)&&(k<=3850*1.1))
	   {

	       x=Tile.at(k).mX;
	       y=Tile.at(k).mY;
	       x=a*x+b;
	       y=a*y+b;
	       correction_map2D_16_die->Fill(x,y) ;

	   }
   }
   }
}
//--------------------
V_max_14 = correction_map2D_14_die->GetMaximum();
correction_map2D_14_die->Scale(1/V_max_14);
V_max_16 = correction_map2D_16_die->GetMaximum();
correction_map2D_16_die->Scale(1/V_max_16);

V_min_14 = correction_map2D_14_die->GetMinimum(0.0);
V_min_16 = correction_map2D_16_die->GetMinimum(0.0);

for(int i=0;i<correction_map2D_14_die->GetNbinsX(); i++)
  {
	   for(int j=0;j<correction_map2D_14_die->GetNbinsY(); j++)
	   {
	       if(correction_map2D_14_die->GetBinContent(i,j)==0)
	       {
		   correction_map2D_14_die->SetBinContent(i,j,V_min_14);
	       }
	       if(correction_map2D_16_die->GetBinContent(i,j)==0)
	       {
		   correction_map2D_16_die->SetBinContent(i,j,V_min_16);
	       }
  }
}

correction_map2D_14_die->Draw("COLZ");
correction_map2D_14_die->SetTitle(";X [mm];Y [mm]");
correction_map2D_14_die->GetXaxis()->SetLabelFont(22);
correction_map2D_14_die->GetXaxis()->SetTitleFont(22);
correction_map2D_14_die->GetYaxis()->SetLabelFont(22);
correction_map2D_14_die->GetYaxis()->SetTitleFont(22);
correction_map2D_14_die->SetLineWidth(3);

correction_map2D_16_die->Draw("COLZ");
correction_map2D_16_die->SetTitle(";X [mm];Y [mm]");
correction_map2D_16_die->GetXaxis()->SetLabelFont(22);
correction_map2D_16_die->GetXaxis()->SetTitleFont(22);
correction_map2D_16_die->GetYaxis()->SetLabelFont(22);
correction_map2D_16_die->GetYaxis()->SetTitleFont(22);
correction_map2D_16_die->SetLineWidth(3);

correction_map2D_14_die->Write("corr_map2D_14_die");
correction_map2D_16_die->Write("corr_map2D_16_die");
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++corrected_map++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TH2D* corrected_map_14= new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));
TH2D* corrected_map_16= new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));
//
TH1D* corrected_map_14_x = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));
TH1D* corrected_map_14_y = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));
//
TH1D* corrected_map_16_x = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));
TH1D* corrected_map_16_y = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));

corrected_map_14->Divide(Map2D_14_die,correction_map2D_14_die);
corrected_map_16->Divide(Map2D_14_die,correction_map2D_16_die);

corrected_map_14->Draw("COLZ");
corrected_map_16->Draw("COLZ");

corrected_map_14->Write("corrected_map_14");
corrected_map_16->Write("corrected_map_16");

corrected_map_14_x=corrected_map_14->ProjectionX();
corrected_map_14_y=corrected_map_14->ProjectionY();
//
corrected_map_16_x=corrected_map_16->ProjectionX();
corrected_map_16_y=corrected_map_16->ProjectionY();
//
corrected_map_14_x->Draw("hist");
corrected_map_14_x->SetTitle(";X [mm]; Events");
corrected_map_14_x->GetXaxis()->SetLabelFont(22);
corrected_map_14_x->GetXaxis()->SetTitleFont(22);
corrected_map_14_x->GetYaxis()->SetLabelFont(22);
corrected_map_14_x->GetYaxis()->SetTitleFont(22);
corrected_map_14_x->SetLineWidth(4);
corrected_map_14_x->SetLineColor(kBlue);

corrected_map_14_y->Draw("hist");
corrected_map_14_y->SetTitle(";Y [mm]; Events");
corrected_map_14_y->GetXaxis()->SetLabelFont(22);
corrected_map_14_y->GetXaxis()->SetTitleFont(22);
corrected_map_14_y->GetYaxis()->SetLabelFont(22);
corrected_map_14_y->GetYaxis()->SetTitleFont(22);
corrected_map_14_y->SetLineWidth(4);
corrected_map_14_y->SetLineColor(kBlue);
//
corrected_map_16_x->Draw("hist");
corrected_map_16_x->SetTitle(";X [mm]; Events");
corrected_map_16_x->GetXaxis()->SetLabelFont(22);
corrected_map_16_x->GetXaxis()->SetTitleFont(22);
corrected_map_16_x->GetYaxis()->SetLabelFont(22);
corrected_map_16_x->GetYaxis()->SetTitleFont(22);
corrected_map_16_x->SetLineWidth(4);
corrected_map_16_x->SetLineColor(kBlue);

corrected_map_16_y->Draw("hist");
corrected_map_16_y->SetTitle(";Y [mm]; Events");
corrected_map_16_y->GetXaxis()->SetLabelFont(22);
corrected_map_16_y->GetXaxis()->SetTitleFont(22);
corrected_map_16_y->GetYaxis()->SetLabelFont(22);
corrected_map_16_y->GetYaxis()->SetTitleFont(22);
corrected_map_16_y->SetLineWidth(4);
corrected_map_16_y->SetLineColor(kBlue);

corrected_map_14_x->Write("corrected_map_14_x");
corrected_map_14_y->Write("corrected_map_14_y");
corrected_map_16_x->Write("corrected_map_16_x");
corrected_map_16_y->Write("corrected_map_16_y");

}

    return 0;
}
