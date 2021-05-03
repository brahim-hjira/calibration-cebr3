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
#include <TGraph2D.h>
#include <TRandom.h>
using namespace std;

void print_number_events_pixel_acquisition(std::vector<TuileEvt> TileA, int number_pixels_axis,int p)
{

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

     //TH1D* Map1D = new TH1D(" ", " ",number_pixels_axis,0. ,(31. ));
     float x_min=500;
     float x_max=8500;
     float y_min=0.;
     float y_max=31.;
     float a=(y_max-y_min)/(8000.);
     float b=y_min-a*x_min;
     float x=0;
     float y=0;
     float h=(x_max-x_min)/(number_pixels_axis);


//------------------------------------------------------------------
for(int k=0;k<TileA.size();k++)//TileA.size()//13627070, 2e6
{				     //100000//463934

    if(p==0)
     {
	 x=TileA.at(k).mX;
	 y=TileA.at(k).mY;
	 x=a*x+b;
	 y=a*y+b;
	 Map2D_14_die->Fill(x,y) ;
	 Map1D_14_die_x->Fill(x) ;
	 Map1D_14_die_y->Fill(y) ;

     }else{

	if ((TileA.at(k).photons >= p*(0.95))&&(TileA.at(k).photons <= p*1.05))// if ((k>=3850*0.9)&&(k<=3850*1.1))
	{                                       //p*(0.95))&&(TileA.at(k).photons <= p*1.05))
						       //5250)&&(TileA.at(k).photons <= 6775))

		 x=TileA.at(k).mX;
		 y=TileA.at(k).mY;
		 x=a*x+b;
		 y=a*y+b;
		 Map2D_14_die->Fill(x,y) ;
		 Map1D_14_die_x->Fill(x) ;
		 Map1D_14_die_y->Fill(y) ;
	}

     }
//--------------------------------------------------------
    if (TileA.at(k).nbDie>=16)
    {
	if(p==0)
	 {
	    x=TileA.at(k).mX;
	    y=TileA.at(k).mY;
	    x=a*x+b;
	    y=a*y+b;
	    Map2D_16_die->Fill(x,y) ;
	    Map1D_16_die_x->Fill(x) ;
	    Map1D_16_die_y->Fill(y) ;

	 }else{

	    if ((TileA.at(k).photons >= p*(0.95))&&(TileA.at(k).photons <= p*1.05))// if ((k>=3850*0.9)&&(k<=3850*1.1))
	    {                                       //p*(0.95))&&(TileA.at(k).photons <= p*1.05))
							   //5250)&&(TileA.at(k).photons <= 6775))

		     x=TileA.at(k).mX;
		     y=TileA.at(k).mY;
		     x=a*x+b;
		     y=a*y+b;
		     Map2D_16_die->Fill(x,y) ;
		     Map1D_16_die_x->Fill(x) ;
		     Map1D_16_die_y->Fill(y) ;
	    }
    }

    }
}
//--------------------------------------------------------------------------
Map2D_14_die->Draw("COLZ");
Map2D_14_die->SetTitle(";X [mm];Y [mm]");
Map2D_14_die->GetXaxis()->SetLabelFont(22);
Map2D_14_die->GetXaxis()->SetTitleFont(22);
Map2D_14_die->GetYaxis()->SetLabelFont(22);
Map2D_14_die->GetYaxis()->SetTitleFont(22);
Map2D_14_die->SetLineWidth(3);

Map1D_14_die_x->Draw("histo ");
Map1D_14_die_x->SetTitle(";X [mm]; Events");
Map1D_14_die_x->GetXaxis()->SetLabelFont(22);
Map1D_14_die_x->GetXaxis()->SetTitleFont(22);
Map1D_14_die_x->GetYaxis()->SetLabelFont(22);
Map1D_14_die_x->GetYaxis()->SetTitleFont(22);
Map1D_14_die_x->SetLineWidth(3);
Map1D_14_die_x->SetLineColor(kBlue);

Map1D_14_die_y->Draw("histo ");
Map1D_14_die_y->SetTitle(";Y [mm]; Events");
Map1D_14_die_y->GetXaxis()->SetLabelFont(22);
Map1D_14_die_y->GetXaxis()->SetTitleFont(22);
Map1D_14_die_y->GetYaxis()->SetLabelFont(22);
Map1D_14_die_y->GetYaxis()->SetTitleFont(22);
Map1D_14_die_y->SetLineWidth(3);
Map1D_14_die_y->SetLineColor(kBlue);


Map2D_14_die->Write("Map2D_14_die");
Map1D_14_die_x->Write("Map1D_14_die_x");
Map1D_14_die_y->Write("Map1D_14_die_y");
//--------------------------------------------------------------------------
Map2D_16_die->Draw("COLZ");
Map2D_16_die->SetTitle(";X [mm];Y [mm]");
Map2D_16_die->GetXaxis()->SetLabelFont(22);
Map2D_16_die->GetXaxis()->SetTitleFont(22);
Map2D_16_die->GetYaxis()->SetLabelFont(22);
Map2D_16_die->GetYaxis()->SetTitleFont(22);
Map2D_16_die->SetLineWidth(3);

Map1D_16_die_x->Draw("histo ");
Map1D_16_die_x->SetTitle(";X [mm]; Events");
Map1D_16_die_x->GetXaxis()->SetLabelFont(22);
Map1D_16_die_x->GetXaxis()->SetTitleFont(22);
Map1D_16_die_x->GetYaxis()->SetLabelFont(22);
Map1D_16_die_x->GetYaxis()->SetTitleFont(22);
Map1D_16_die_x->SetLineWidth(3);
Map1D_16_die_x->SetLineColor(kBlue);

Map1D_16_die_y->Draw("histo ");
Map1D_16_die_y->SetTitle(";Y [mm]; Events");
Map1D_16_die_y->GetXaxis()->SetLabelFont(22);
Map1D_16_die_y->GetXaxis()->SetTitleFont(22);
Map1D_16_die_y->GetYaxis()->SetLabelFont(22);
Map1D_16_die_y->GetYaxis()->SetTitleFont(22);
Map1D_16_die_y->SetLineWidth(3);
Map1D_16_die_y->SetLineColor(kBlue);


Map2D_16_die->Write("Map2D_16_die");
Map1D_16_die_x->Write("Map1D_16_die_x");
Map1D_16_die_y->Write("Map1D_16_die_y");

 //Map11D16->GetYaxis()->SetTitle("N#circ events");
 //Map11D16->GetXaxis()->SetTitle("N#circ pixels");

}
 //==========================================
  void print_number_events_pixel_map_of_correction(std::vector<TuileEvt> TileA, int number_pixels_axis,int p)
 {

      gStyle-> SetOptStat(000);
      TCanvas* Map_tuile_event= new TCanvas("Map_tuile_event","Map_tuile_event", 1400, 1000);
      TCanvas* Global=new TCanvas("","", 1600, 1000);
      TH2D* correction_map2D_14_die = new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));
      TH2D* correction_map2D_16_die = new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));
      float x_min=500;
      float x_max=8500;
      float y_min=0.;
      float y_max=31.;
      float a=(y_max-y_min)/(8000.);
      float b=y_min-a*x_min;
      float x=0;
      float y=0;
      float h=(x_max-x_min)/(number_pixels_axis);
      float V_min_14;
      float V_min_16;
      float V_max_14;
      float V_max_16;

 //------------------------------------------------------------------
 for(int k=0;k<TileA.size();k++)//TileA.size()//13627070, 2e6
 {				     //100000//463934

     if(p==0)
      {
	  x=TileA.at(k).mX;
	  y=TileA.at(k).mY;
	  x=a*x+b;
	  y=a*y+b;
	  correction_map2D_14_die->Fill(x,y) ;
      }else{

	 if ((TileA.at(k).photons >= p*(0.95))&&(TileA.at(k).photons <= p*1.05))// if ((k>=3850*0.9)&&(k<=3850*1.1))
	 {                                       //p*(0.95))&&(TileA.at(k).photons <= p*1.05))
							//5250)&&(TileA.at(k).photons <= 6775))

	     x=TileA.at(k).mX;
	     y=TileA.at(k).mY;
	     x=a*x+b;
	     y=a*y+b;
	     correction_map2D_14_die->Fill(x,y) ;
	 }

      }
 //--------------------------------------------------------
     if (TileA.at(k).nbDie>=16)
     {
	 if(p==0)
	  {
	     x=TileA.at(k).mX;
	     y=TileA.at(k).mY;
	     x=a*x+b;
	     y=a*y+b;
	     correction_map2D_16_die->Fill(x,y) ;
	  }else{

	     if ((TileA.at(k).photons >= p*(0.95))&&(TileA.at(k).photons <= p*1.05))// if ((k>=3850*0.9)&&(k<=3850*1.1))
	     {

		 x=TileA.at(k).mX;
		 y=TileA.at(k).mY;
		 x=a*x+b;
		 y=a*y+b;
		 correction_map2D_16_die->Fill(x,y) ;

	     }
     }

     }
 }
 //--------------------------------------------------------------------------
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


 //--------------------------------------------------------------------------
 correction_map2D_16_die->Draw("COLZ");
 correction_map2D_16_die->SetTitle(";X [mm];Y [mm]");
 correction_map2D_16_die->GetXaxis()->SetLabelFont(22);
 correction_map2D_16_die->GetXaxis()->SetTitleFont(22);
 correction_map2D_16_die->GetYaxis()->SetLabelFont(22);
 correction_map2D_16_die->GetYaxis()->SetTitleFont(22);
 correction_map2D_16_die->SetLineWidth(3);

 correction_map2D_14_die->Write("corr_map2D_14_die");
 correction_map2D_16_die->Write("corr_map2D_16_die");

 }
//-------------------------------------------------
  void print_number_events_pixel_corrected_map(int number_pixels_axis)
 {
      TH2D* Map2D_14;
      TH2D* Map2D_16;

      TH2D* CDC_14;
      TH2D* CDC_16;

      TH2D* corrected_map2D_14_die;
      TH2D* corrected_map2D_16_die;


      //------------------------------------------------
	      string rootfilename1 = "/home/hjira/calibration/output/output_root/calibration.root";
	      TFile* f1 = new TFile(rootfilename1.c_str());
	      f1->ls();

	      if (f1->IsZombie()) {
		     std::cout << "Error opening file" << std::endl;
	      }
	      else{
		  std::cout<<"File is opened"<< std::endl;
	      }
      //-----------------------------------------------------

	      Map2D_14 = (TH2D*)f1->Get("Map2D_14_die");
	      Map2D_16 = (TH2D*)f1->Get("Map2D_16_die");
	      CDC_14 = (TH2D*)f1->Get("corr_map2D_14_die");
	      CDC_16 = (TH2D*)f1->Get("corr_map2D_16_die");
//--------------------------------------
TH2D* ratio1= new TH2D(" ", " ", number_pixels_axis, 0. ,(31. ), number_pixels_axis,0. ,(31. ));

	      ratio1->Divide(Map2D_14,CDC_14);
	     // corrected_map2D_14_die->Draw("colz");
//ratio1->Draw("colz");
//ratio1->Write("CC");

	      //TH2D h3 = (*Map2D_14_die)/(*correction_map2D_14_die);
	      //corrected_map2D_14_die->Write("corrected_map2D_14_die");
	       //TH2D *h3 = (TH2D*)Map2D_14_die->Clone("h3");
	       //h3->Divide(Map2D_14_die,correction_map2D_14_die);
	       //h3->Draw("colz");
}

 /** PrintEnergyPeaks function
   * @brief: Display the energy peak positions in a table and print it in the shell
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
//--------------------------------------------------
  void print_photons_Spectrum_filter(std::vector<TuileEvt> TileA)// p=position
   {
    TCanvas *ESpectrum_Canvas = new TCanvas("ESpectrum", "ESpectrum", 1600, 1000);
    //TH1F *Photons_Spectrum_pic = new TH1F(" Photons_Spectrum_pic", " ", 1000,  p*0.95, p*1.05);// p*0.95, p*1.05
    TH1F *Photons_Spectrum= new TH1F("Photons_Spectrum ", " ",1000 , 0, 12000);

			      for(int k=0;k<TileA.size();++k)
			       {
				  Photons_Spectrum->Fill(TileA.at(k).photons);

			       }

			     /// if(p==0){
				  //ESpectrum_Canvas->cd();
				  Photons_Spectrum->Draw("histo ");
				  Photons_Spectrum->SetTitle(";Channel Numbers; Events");
				  Photons_Spectrum->GetXaxis()->SetLabelFont(22);
				  Photons_Spectrum->GetXaxis()->SetTitleFont(22);
				  Photons_Spectrum->GetYaxis()->SetLabelFont(22);
				  Photons_Spectrum->GetYaxis()->SetTitleFont(22);
				  Photons_Spectrum->SetLineWidth(2);
				  Photons_Spectrum->SetLineColor(kBlue);
				  Photons_Spectrum->Write("Photons_Spectrum");
				  //ESpectrum_Canvas->Write("Photons_canvas");

   }
































  void print_number_htimestamp(std::vector<TuileEvt> TileA)
{

    int c=0;
    for(int k=0;k<TileA.size()-1;k++)
      {

			   if (TileA.at(k+1).htimestamp>TileA.at(k).htimestamp)
			   {
			   c++;

			   //

			   //cout<<c<<endl;
			   }

     }
     //cout<<c<<endl;
    cout<<"number_htimestamp = "<<c<<endl;

}
//---------------------------------
  void print_mean_temperature(std::vector<TuileEvt> TileA)
  {
      double Tsum=0;
      vector<double> Tarr;
      int size=TileA.size();

      for(int k=0;k<size;k++)
      {
	  Tsum+=(double)(TileA.at(k).temper/100);
      }

  cout<<"Tmoy= "<<Tsum/TileA.size()<<endl;

  }


