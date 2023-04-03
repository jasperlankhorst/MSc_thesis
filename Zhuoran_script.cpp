#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TColor.h"
#include <vector>
#include "TPad.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TMath.h"

#include "TLatex.h"
#include "TText.h"
#include "TSystem.h"
#include "TStyle.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>
#include <ios>
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TGraphQQ.h"
#include "TransformTool/TransformTool/HistoTransform.h"

using namespace std;
bool g_doSystematics = false;

//////////////////////////////////////////
//scale factor provide (if additoinal scale factors should be applied to the input histograms)
class scaleFactorProvider {  
private:
  std::map<std::string, double> m_scaleFactors;
public:
  scaleFactorProvider(std::string fileName);
  double getScaleFactor(std::string sampleName);
};

//////////////////////////////////////////
//process structure with all the needed informations for plotting
struct process{
  vector<string> v_dsids;
  string name;
  string tex_name;
  Int_t color;
  TH1D *hist;
  TH1D *hist_stack;
  TH1D *h_rat_sys_up;
  TH1D *h_rat_sys_down;
  TH1D *h_rat_sys_modelling_up;
  TH1D *h_rat_sys_modelling_down;
};

//////////////////////////////////////////
//declare the functions needed
//////////////////////////////////////////
void binsHist(TH1* H_1);
TH1* getDataMCComparison(TH1 *h_data, TH1 *h_MC);
TH1* getRatioStatUncertaintyBand(TH1 *h_MC, int mode);
void createStackPlots(string INPUT_FILE_NAME, string year, string out, string s_sample1, vector<process> v_processes, string s_region, string s_pTV, string s_lepChannel, string s_ctag, string s_varName, string s_legX, string s_binWidth, scaleFactorProvider* sfp = NULL, TH1D* h_chi2P=NULL, TH1D* h_ksP=NULL);
void fill_processHistos(vector<process> &v_processes, TFile *f, string s_region, string s_pTV, string s_lepChannel, string s_ctag, string s_varName, scaleFactorProvider* sfp = NULL);
void fill_stackedHistos(vector<process> &v_processes);
void selectionPlots_stacked_syst( string INPUT_FILE_NAME, string year, string out, bool applySFs);
void define_processes(vector<process> &v_processes);
void save_yields(TH1* hdat, vector<process> &v_processes, string name);
void addSystematics(TH1* hnom, TH1* halt, TH1* hup, TH1* hdown);
void fillSystematicHistograms(TFile* f, TString name, TH1* hnom, TH1* hup, TH1* hdown, TH1* hMup, TH1* hMdown);
TH1D* getTotalSystematics(vector<process> &v_processes, int ud, TH1* h_add);
TGraphAsymmErrors* getErrorBand(TH1*hup, TH1* hdown);
TGraphAsymmErrors* getTotalUncertaintyBand(TH1D* htot, TH1D* hup, TH1D* hdown);
TH1* HistSignificance(vector<process> &v_process, int sigIdx);
TH1* HistSoverB(vector<process> &v_processes, int sigIdx);
double ComputeSignificance(vector<process> &v_process, int sigIdx);
double ComputeSigOverBkg(vector<process> &v_processes, int sigIdx);
double ComputeSigOverBkg_blind(vector<process> &v_processes, int sigIdx);

//gROOT->ProcessLine(".L TransformTool/Root/HistoTransform.cxx+");

//////////////////////////////////////////
void selectionPlots_stacked_syst( string INPUT_FILE_NAME, string year, string out, bool applySFs){

  //gROOT->ProcessLine(".include TransformTool/");
  //gROOT->ProcessLine(".L TransformTool/Root/HistoTransform.cxx+");
  gSystem->Load("libTransformTool.so");

  if( year!="all" &&  year!="2015_2016" &&  year!="2017" &&  year!="2018" ){
    std::cerr << "<!> The year argument is not one of the allowed values. Please check." << std::endl;
    return;
  }

  //BM: Histograms to save the p-values of the chi2 and the KS test
  TH1D* h_chi2P = new TH1D("h_chi2P", "", 10, 0, 1);
  TH1D* h_ksP = new TH1D("h_ksP", "", 10, 0, 1);

  //vector to strore process information
  vector<process> v_processes;
  define_processes(v_processes); 
  
  vector<string> v_chan;
  v_chan.push_back("elmu");

  vector<string> v_reg;
  v_reg.push_back("R2"); // SR 2jets 
  v_reg.push_back("R3"); // SR 3jets
  v_reg.push_back("R4");

  vector<string> v_pTV;
  v_pTV.push_back("Low");
  v_pTV.push_back("Mid");
  //v_pTV.push_back("High");
  v_pTV.push_back("Mid+High");
  v_pTV.push_back("Inclusive");
  

  //v_Ctag_CR.push_back("NoRequire");

  vector<string> v_Ctag;
  //v_Ctag_SR.push_back("default"); // no ctagging requirements
  string ctag_temp = "";
  v_Ctag.push_back("TT");
  v_Ctag.push_back("TL+LT");
  //v_Ctag.push_back("TN+NT");
  

  //additional scale factors for the histograms per sample per region? 
  std::map<std::string, scaleFactorProvider*> v_sfprov;
  for(auto reg: v_reg){
    //if(applySFs){
      for(auto ctag: v_Ctag)v_sfprov.insert(std::pair<std::string, scaleFactorProvider*>(reg, new scaleFactorProvider(std::string("sfs_S" + reg + "_" + ctag + ".txt"))));
    //}else{
    //  for(auto ctag: v_Ctag)v_sfprov.insert(std::pair<std::string, scaleFactorProvider*>(reg, NULL));
    //}

  }

  for(auto lep: v_chan){
    for(auto reg: v_reg){
      for(auto pTV: v_pTV){
        std::cout << "check region: " << reg << std::endl;
        for(auto Ctag: v_Ctag){
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "minDRBjets","min #Delta R(c, j)", "", v_sfprov[reg], h_chi2P, h_ksP);
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "pTBB","p_{T}^{cc}", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "mBB","m_{cc} [GeV]", "", v_sfprov[reg], h_chi2P, h_ksP);
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "dRBB","#Delta R(cc)", "", v_sfprov[reg], h_chi2P, h_ksP);
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "MET","MET [GeV]", "", v_sfprov[reg], h_chi2P, h_ksP);
          /*createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "mB1","m_{c1} [GeV]", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "mB2","m_{c2} [GeV]", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "softMET","p_{T}^{miss,st} [GeV]", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "MEff","m_{eff} [GeV]", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "dEtaBB","#Delta #eta(cc)", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "dPhiBB","#Delta #phi(cc)", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "dPhiVBB","#Delta #phi(V,H)", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "dEtaVBB","#Delta #eta(V,H)", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "pTB1","p_{T}^{c1}", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "pTB2","p_{T}^{c2}", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "bin_bTagB1", "bin_bTagB1", "", v_sfprov[reg], h_chi2P, h_ksP);
          createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "bin_bTagB2", "bin_bTagB2", "", v_sfprov[reg], h_chi2P, h_ksP);*/
          // if(pTV!="High")createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "VHccBDT", "VHccBDT", "", v_sfprov[reg], h_chi2P, h_ksP);
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "nJetPt2030", "nJetPt2030", "", v_sfprov[reg], h_chi2P, h_ksP);
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "sumPtAddJets", "sumPtAddJets", "", v_sfprov[reg], h_chi2P, h_ksP);
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "mBBJ","m_{ccj}", "", v_sfprov[reg], h_chi2P, h_ksP);
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "mBBJFSR01","m_{ccj, 01}", "", v_sfprov[reg], h_chi2P, h_ksP);
          //createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag, "mBBJFSR02","m_{ccj, 02}", "", v_sfprov[reg], h_chi2P, h_ksP);
        }
        /*if(reg=="R3" || reg=="R4"){
          for(auto Ctag1: v_Ctag){
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "pTJ3","p_{T}^{j3}", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "etaJ3","#eta(j3)", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "mBBJ","m_{ccj}", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "mB1J3","m_{c1,j3}", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "mB2J3","m_{c2,j3}", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "dRBBJ3","#Delta R(H,j3)", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "dRB1J3","#Delta R(c1,j3)", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "dRB2J3","#Delta R(c2,j3)", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag1, "bin_bTagJ3", "bin_bTagJ3", "", v_sfprov[reg], h_chi2P, h_ksP);
          }
        }
        if(reg=="R4"){
          for(auto Ctag2: v_Ctag){
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag2, "mBBJ4","m_{ccj4}", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag2, "dRB1J4","#Delta R(c1,j4)", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag2, "dRB2J4","#Delta R(c2,j4)", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag2, "mB1J4","m_{c1,j4}", "", v_sfprov[reg], h_chi2P, h_ksP);
            createStackPlots(INPUT_FILE_NAME, year, out, "data", v_processes, reg, pTV, lep, Ctag2, "mB2J4","m_{c2,j4}", "", v_sfprov[reg], h_chi2P, h_ksP);
          }
        }*/
      }
    }
  }  

  //BM: Finally, plot the chi2 and ks values
  TCanvas *cpval = new TCanvas("cpval","", 900, 600);
  h_chi2P->GetXaxis()->SetTitle("p-value");
  h_chi2P->GetYaxis()->SetTitle("plots/bin");
  h_chi2P->SetLineWidth(2);
  h_chi2P->SetLineColor(TColor::GetColor("#E03A3E"));
  h_ksP->SetLineWidth(2);
  h_ksP->SetLineColor(TColor::GetColor("#3598DB"));
  double maximum = h_chi2P->GetMaximum();
  if (h_ksP->GetMaximum() > maximum) maximum = h_ksP->GetMaximum();
  h_chi2P->GetYaxis()->SetRangeUser(0, maximum*1.45);
  h_chi2P->Draw("hist");
  h_ksP->Draw("hist same");
  TLatex l_LS;
  l_LS.SetNDC();
  l_LS.SetTextColor(1);
  l_LS.SetTextSize(0.045);
  l_LS.DrawLatex(0.16, 0.87, "#font[72]{ATLAS} #font[42]{Internal}");
  l_LS.SetTextSize(0.035);
  string s_label = "p-value summary";
  //  s_label += std::to_string(i)+"-lepton channel";
  l_LS.DrawLatex(0.16, 0.81, s_label.c_str());
  TLegend *LR = new TLegend(0.16, 0.76, 0.65, 0.70);
  LR->SetNColumns(2);
  LR->SetFillColor(0);
  LR->SetLineWidth(0);
  LR->SetBorderSize(1);
  LR->SetLineColor(kWhite);
  LR->AddEntry(h_chi2P, "#chi^{2} test", "l");
  LR->AddEntry(h_ksP, "KS test", "l");
  LR->Draw();
  string s_saveName = out + "/C_PValueSummary.pdf";
  cpval->SaveAs(s_saveName.c_str());

}

/////////////////////////////////////////////
//   Helper functions
/////////////////////////////////////////////

//////////////////////////////////////////
//add all the relevant processes and information about them
void define_processes(vector<process> &v_processes){
  
  v_processes.push_back({{"Zcc", "Zbb"},"Zhf", "Z+hf", TColor::GetColor("#3F65C5"), NULL, NULL});
  v_processes.push_back({{"Zl"},"Zlf", "Z+lf", TColor::GetColor("#A4CBFA"), NULL, NULL});
  v_processes.push_back({{"Zcl", "Zbl", "Zbc"},"Zmf", "Z+mf", TColor::GetColor("#7197F8"), NULL, NULL});
  
  v_processes.push_back({{"Wl"},"Wlf", "W+lf", TColor::GetColor("#B1FCA3"), NULL, NULL});
  v_processes.push_back({{"Wcl", "Wbl", "Wbc"},"Wmf", "W+mf", TColor::GetColor("#80C972"), NULL, NULL});
  v_processes.push_back({{"Wcc", "Wbb"},"Whf", "W+hf", TColor::GetColor("#529742"), NULL, NULL});
  v_processes.push_back({{"stopsbc", "stopsbb", "stopsbl", "stoptbb", "stoptbc", "stoptbl", "stopscc", "stopscl", "stopsl", "stoptcc", "stoptcl", "stoptl"}, "stopst", "stop(s/t)", TColor::GetColor("#F2E0C2"), NULL, NULL});
  v_processes.push_back({{"ttbarbl", "ttbarbc", "ttbarbb", "stopWtDSbb", "stopWtDSbc", "stopWtDSbl"},"top(b)","top(b)", TColor::GetColor("#F7CE46"), NULL, NULL}); // need to check
  v_processes.push_back({{"ttbarcc", "ttbarcl", "ttbarl", "stopWtDScc", "stopWtDScl", "stopWtDSl"},"topother","top(other)", TColor::GetColor("#FEFE61"), NULL, NULL});//need to check
  v_processes.push_back({{"ZZbkg", "WZlephadbkg", "ggZZbkg", "WW", "ggWW", "WZhadlep"},"VVbkg","VZ(Bkg)+VW", TColor::GetColor("#CCCCCC"), NULL, NULL});  
  v_processes.push_back({{"ZZbb", "WZlephadbb", "ggZZbb"},"VZbb","VZ(#rightarrow b#bar{b})", TColor::GetColor("#999999"), NULL, NULL});
  v_processes.push_back({{"WZlephadcc", "ZZcc", "ggZZcc"},"VZcc","VZ(#rightarrow c#bar{c})", TColor::GetColor("#666666"), NULL, NULL});
  v_processes.push_back({{"ggZllH125", "qqZllH125", "ggZvvH125", "qqZvvH125", "qqWlvH125"},"VHbb","VHb#bar{b}", TColor::GetColor("#5C0E93"), NULL, NULL});

  // check for signal : cc version
  v_processes.push_back({{"ggZllH125cc", "qqZllH125cc", "ggZvvH125cc", "qqZvvH125cc", "qqWlvH125cc"},"signal","Signal", TColor::GetColor("#E03A3E"), NULL, NULL});
}


//////////////////////////////////////////
void binsHist(TH1* H_1){
  H_1->SetStats(0);
  /*
  //UF treatement
  H_1->SetBinContent(1, H_1->GetBinContent(0)+H_1->GetBinContent(1));
  H_1->SetBinError(1, sqrt(pow(H_1->GetBinError(0),2)+pow(H_1->GetBinError(1),2)));
  H_1->SetBinContent(0, 0 );
  H_1->SetBinError(0, 0);
  //OF treatement
  H_1->SetBinContent(H_1->GetNbinsX(), H_1->GetBinContent(H_1->GetNbinsX())+H_1->GetBinContent(H_1->GetNbinsX()+1));
  H_1->SetBinError(H_1->GetNbinsX(), sqrt(pow(H_1->GetBinError(H_1->GetNbinsX()),2)+pow(H_1->GetBinError(H_1->GetNbinsX()+1),2)));
  H_1->SetBinContent(H_1->GetNbinsX()+1, 0);
  H_1->SetBinError(H_1->GetNbinsX()+1, 0);*/
  H_1->ClearUnderflowAndOverflow();
  //Rebinning
  TString hname = H_1->GetName();

  if(hname.Contains("mBB")) H_1->Rebin(5);
  if(hname.Contains("pTV")) H_1->Rebin(4);
  if(hname.Contains("dRB")) H_1->Rebin(10); // Previously 4
  if(hname.Contains("minDR")) H_1->Rebin(10);
  if(hname.Contains("MET")) H_1->Rebin(4);
  if(hname.Contains("dEtaBB")) H_1->Rebin(10); // Previously 4
  if(hname.Contains("dPhiBB")) H_1->Rebin(10); // Previously 4
  if(hname.Contains("dEtaVBB")) H_1->Rebin(8); // new
  if(hname.Contains("dPhiVBB")) H_1->Rebin(2); // new
  if(hname.Contains("pTB1")) H_1->Rebin(4);
  if(hname.Contains("pTBB")) H_1->Rebin(4);
  if(hname.Contains("pTB2")) H_1->Rebin(10);
  if(hname.Contains("pTJ3")) H_1->Rebin(5); //previously 5
  //if(hname.Contains("mBBJ")) H_1->Rebin(2); // new
  if(hname.Contains("MEff")) H_1->Rebin(5);
  if(hname.Contains("etaJ3")) H_1->Rebin(8);
  if(hname.Contains("sumPtAddJets")) H_1->Rebin(10);

  if(hname.Contains("mB1")) H_1->Rebin(5);
  if(hname.Contains("mB2")) H_1->Rebin(5);

  if(hname.Contains("MEff")) H_1->GetXaxis()->SetRangeUser(200, 1000);
  if(hname.Contains("softMET")) H_1->GetXaxis()->SetRangeUser(0,40);
  if(hname.Contains("BDT"))H_1->GetXaxis()->SetRangeUser(-1,1);
  if(hname.Contains("dRB"))H_1->GetXaxis()->SetRangeUser(0.4, 3.2);
  if(hname.Contains("minDR"))H_1->GetXaxis()->SetRangeUser(0.4, 3.2);
  if(hname.Contains("dEtaBB"))H_1->GetXaxis()->SetRangeUser(0, 3.2);
  if(hname.Contains("pTBB"))H_1->GetXaxis()->SetRangeUser(50,550);
  if(hname.Contains("MET"))H_1->GetXaxis()->SetRangeUser(0,400);
  //if(hname.Contains("BDT"))H_1->Rebin(10);



}//binsHist


//////////////////////////////////////////
TH1* getDataMCComparison(TH1 *h_data, TH1 *h_MC) {

  TH1 *h_comp = (TH1*)h_data->Clone();
  h_comp->SetDirectory(0);
  h_comp->Reset();

  for(int bin=1; bin<=h_comp->GetNbinsX(); bin++) {

    double nData = h_data->GetBinContent(bin);
    double eData = h_data->GetBinError(bin);
    double nMC = h_MC->GetBinContent(bin);

    if(nData > 1e-6 && eData > 1e-6 && nMC > 1e-6) {
      double nComp = nData / nMC;
      double eComp = eData / nMC;
      h_comp->SetBinContent(bin, nComp);
      h_comp->SetBinError(bin, eComp);
    }
  }
  return h_comp;
}

TH1* HistSignificance(vector<process> &v_processes, int sigIdx){
  TH1 *h_comp = (TH1*)v_processes.at(sigIdx).hist->Clone();
  h_comp->SetDirectory(0);
  h_comp->Reset();

  double Significance = 0.0;
  double sig_temp = 0.0;
  double nSig = 0.0;
  double nBkg = 0.0;

  for(int bin = 1; bin<=h_comp->GetNbinsX(); bin++){
    nSig = v_processes.at(sigIdx).hist->GetBinContent(bin);
    nBkg = v_processes.at(v_processes.size()-1).hist_stack->GetBinContent(bin) - nSig;
    
    sig_temp = 0.0;
    if(nBkg>1e-6){
      Significance += 2 * ((nSig+nBkg)*log(1+nSig/nBkg)-nSig);
      sig_temp = sqrt(2 * ((nSig+nBkg)*log(1+nSig/nBkg)-nSig));
    }
    h_comp->SetBinContent(bin, sig_temp*100.0);
    sig_temp = h_comp->GetBinContent(bin);
    std::cout<< "Check significance bin " << bin << ", value: " << sig_temp <<endl;
  }
  double final_Sig = sqrt(Significance);
  return h_comp;
}

TH1* HistSoverB(vector<process> &v_processes, int sigIdx){
  TH1 *h_comp = (TH1*)v_processes.at(sigIdx).hist->Clone();
  h_comp->SetDirectory(0);
  h_comp->Reset();

  double SoverB = 0.0;
  double nSig = 0.0;
  double nBkg = 0.0;

  for(int bin = 1; bin<=h_comp->GetNbinsX(); bin++){
      nSig = v_processes.at(sigIdx).hist->GetBinContent(bin);
      nBkg = v_processes.at(v_processes.size()-1).hist_stack->GetBinContent(bin) - nSig;
      
      SoverB = 0.0;
      if(nBkg>1e-6)SoverB = nSig / nBkg;
      h_comp->SetBinContent(bin, SoverB*300.0);
      SoverB = h_comp->GetBinContent(bin);
      std::cout<< "Check SoverB bin " << bin << ", value: " << SoverB <<endl;
  }
  return h_comp;
}

double ComputeSignificance(vector<process> &v_processes, int sigIdx){
  TH1 *h_comp = (TH1*)v_processes.at(sigIdx).hist->Clone();
  h_comp->SetDirectory(0);
  h_comp->Reset();

  double Significance = 0.0;
  double sig_temp = 0.0;
  double nSig = 0.0;
  double nBkg = 0.0;

  for(int bin = 1; bin<=h_comp->GetNbinsX(); bin++){
    nSig = v_processes.at(sigIdx).hist->GetBinContent(bin);
    nBkg = v_processes.at(v_processes.size()-1).hist_stack->GetBinContent(bin) - nSig;
    sig_temp = 0.0;
    if(nBkg>1e-6){
      Significance += 2 * ((nSig+nBkg)*log(1+nSig/nBkg)-nSig);
      sig_temp = sqrt(2 * ((nSig+nBkg)*log(1+nSig/nBkg)-nSig));
    }
    //h_comp->Fill(bin, sig_temp);
  }
  double final_Sig = sqrt(Significance);
  return final_Sig;
}

double ComputeSigOverBkg(vector<process> &v_processes, int sigIdx){
  TH1 *h_comp = (TH1*)v_processes.at(sigIdx).hist->Clone();
  h_comp->SetDirectory(0);
  h_comp->Reset();

  double SoverB = 0.0;
  double sig_temp = 0.0;
  double nSig = 0.0;
  double nBkg = 0.0;
  nSig = v_processes.at(sigIdx).hist->Integral();
  nBkg = v_processes.at(v_processes.size()-1).hist_stack->Integral() - nSig;
  if(nBkg>0)SoverB = nSig/nBkg;
  return SoverB;
}

double ComputeSigOverBkg_blind(vector<process> &v_processes, int sigIdx){
  TH1 *h_comp = (TH1*)v_processes.at(sigIdx).hist->Clone();
  h_comp->SetDirectory(0);
  //h_comp->Reset();

  double SoverB = 0.0;
  double sig_temp = 0.0;
  double nSig = 0.0;
  double nBkg = 0.0;

  nSig = v_processes.at(sigIdx).hist->Integral(3,9);
  nBkg = v_processes.at(v_processes.size()-1).hist_stack->Integral(3, 9) - nSig;

  if(nBkg>0)SoverB = nSig/nBkg;
  return SoverB;
}


//////////////////////////////////////////
TH1* getRatioStatUncertaintyBand(TH1 *h_MC, int mode){

  TString name = h_MC->GetName();
  if(mode == 0){
    name.Append("MCStat_up");
  }else{
    name.Append("MCStat_down");
  }

  TH1 *h_comp = (TH1*)h_MC->Clone();
  h_comp->SetDirectory(0);
  h_comp->Reset();

  for(int bin=1; bin<=h_comp->GetNbinsX(); bin++) {

    double nMC = h_MC->GetBinContent(bin);
    double eMC = h_MC->GetBinError(bin);

    if(nMC > 1e-6) {
      double diff = fabs( eMC/nMC );
      if( mode == 1 ) h_comp->SetBinContent(bin, 1 + diff);
      else h_comp->SetBinContent(bin, 1 - diff);
    }else{
      h_comp->SetBinContent(bin, 1.0);
    }
  }
  return h_comp;

}

//////////////////////////////////////////
void createStackPlots( string INPUT_FILE_NAME, string year, string out, string s_sample1, vector<process> v_processes, string s_region, string s_pTV, string s_lepChannel, string s_ctag, string s_varName, string s_legX, string s_binWidth, scaleFactorProvider* sfp, TH1D* h_chi2P, TH1D* h_ksP){

  TGaxis::SetMaxDigits(4);
  TH1::SetDefaultSumw2();
  
  //TO EDIT: Input file name
  string s_inputFileName = INPUT_FILE_NAME;
  //open the right input file
  TFile* f_input = new TFile(s_inputFileName.c_str(),"READ"); 
  TH1D* h_data;
  TH1D* h_SigBkg1;
  TH1D* h_SigBkg2;

  //this is for getting the first histogram, usually data
  string s_tmp_lepName = "";
  if(s_lepChannel == "emu"){
    s_tmp_lepName = "elmu";
  }else{
    s_tmp_lepName = s_lepChannel;
  }
  
  string s_h1_name = "hist_";
  s_h1_name.append(s_sample1);
  s_h1_name.append("_");
  s_h1_name.append(s_tmp_lepName);
  s_h1_name.append("_");
  s_h1_name.append(s_region);
  s_h1_name.append("_");
  string s_h1_name2 = s_h1_name;
  string s_h1_name_temp = s_h1_name;
  s_h1_name.append(s_pTV);
  s_h1_name.append("_");
  s_h1_name.append(s_varName);
  string s_h1_name_save = s_h1_name;
  s_h1_name.append("_");
  s_h1_name.append(s_ctag);

  std::map<int, string> SingleJet{{0, "N"}, {1, "L"}, {2, "T"},{3, "B"}, };

  vector<string> v_ctag_temp;
  int i_temp = 0;
  string ctag_temp1 = "";
  int count_T, count_T_leading, count_L, count_L_leading;
  int tag_j1,tag_j2,tag_j3;
  vector<string> v_ptv_temp;
  
  for(int Label_i=0; Label_i<64; Label_i++){
    if(Label_i%4 != 0 && s_region == "R2")continue;
    count_T = 0;
    count_L = 0;
    tag_j1 = Label_i>>4;
    if(tag_j1>2.5)continue;
    else if(tag_j1>1.5)count_T++;
    else if(tag_j1>0.5)count_L++;

    tag_j2 = (Label_i - tag_j1*16)>>2;
    if(tag_j2>2.5)continue;
    else if(tag_j2>1.5)count_T++;
    else if(tag_j2>0.5)count_L++;

    tag_j3 = Label_i - tag_j1*16 - tag_j2*4;
    if(tag_j3>2.5)continue;
    else if(tag_j3>1.5)count_T++;
    else if(tag_j3>0.5)count_L++;

    ctag_temp1 = std::to_string(Label_i);
    if(count_T>1.5){
      if(s_ctag=="TT")v_ctag_temp.push_back(ctag_temp1);
    }else if(count_T>0.5 && count_L>0.5){
      if(s_ctag=="TL+LT")v_ctag_temp.push_back(ctag_temp1);
    }else if(count_T>0.5){
      if(s_ctag=="TN+NT")v_ctag_temp.push_back(ctag_temp1);
    }
  }
  if(s_pTV=="Inclusive"){
    v_ptv_temp.push_back("Low");
    v_ptv_temp.push_back("Mid");
    v_ptv_temp.push_back("High");
    //v_ptv_temp.push_back(s_pTV);
  }else if(s_pTV=="Mid+High"){
    v_ptv_temp.push_back("Mid");
    v_ptv_temp.push_back("High");
  }else{
    v_ptv_temp.push_back(s_pTV);
  }

  int flag_hist=0;
  string yLabel = "events";
  for(auto ctag_temp: v_ctag_temp){
    for(auto ptv_temp: v_ptv_temp){
      TH1D* h1;
      //just for check EventWeight
      //h1 = new TH1D("h1", "h1", 300, -200.0, 100.0);
      s_h1_name_temp = s_h1_name2;
      s_h1_name_temp.append(ptv_temp);
      s_h1_name_temp.append("_");
      s_h1_name_temp.append(s_varName);
      s_h1_name_temp.append("_");
      s_h1_name_temp.append(ctag_temp);
      //which histograms to get
      //do everything with the data histogram
      std::cout << "Clone: " << s_h1_name_temp << std::endl;
      if(f_input->Get(s_h1_name_temp.c_str())==NULL)continue;
      h1 = (TH1D *)f_input->Get(s_h1_name_temp.c_str())->Clone();
      std::cout << "Finished cloneing " << std::endl;
      h1->SetDirectory(0);
      binsHist(h1);
      h1->SetMarkerStyle(8);
      h1->SetMarkerColor(kBlack);
      h1->SetMarkerSize(1.0);
      //h1->GetXaxis()->SetRangeUser(0,9);

    
      if( !(s_binWidth == "") ){
        yLabel.append("/");
        yLabel.append(s_binWidth);
      }

      h1->GetXaxis()->SetTitle(s_legX.c_str());
      h1->GetYaxis()->SetTitle(yLabel.c_str());
      h1->GetXaxis()->SetLabelOffset(999);
      h1->GetXaxis()->SetLabelSize(0);
      std::cout << "Finished h1 " << std::endl;
      if(flag_hist==0){
        h_data = (TH1D *)h1->Clone(s_h1_name.c_str());
        h_SigBkg1 = (TH1D *)h1->Clone(s_h1_name.c_str());
        flag_hist = 1;
      }else{
        h_SigBkg1->Add(h1);
        h_data->Add(h1);

      }
      std::cout << "Finished data" << std::endl;
      // MC
      fill_processHistos(v_processes, f_input, s_region, ptv_temp, s_lepChannel, ctag_temp, s_varName, sfp);
      if(v_processes.at(0).hist == NULL){
        cout << "weird things happen here" << endl;
      }
    }
  }
  if(s_pTV=="Inclusicve")s_pTV="150GeV+";
  else if(s_pTV=="High")s_pTV="400GeV+";
  else if(s_pTV=="Mid")s_pTV="250-400GeV";
  else if(s_pTV=="Low")s_pTV="150-250GeV";
  else if(s_pTV=="Mid+High")s_pTV="250GeV+";
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerColor(kBlack);
  h_data->SetMarkerSize(1.0);
  //h_data->SetLineWidth(0);
  h_data->GetXaxis()->SetTitle(s_legX.c_str());
  h_data->GetYaxis()->SetTitle(yLabel.c_str());
  h_data->GetXaxis()->SetLabelOffset(999);
  h_data->GetXaxis()->SetLabelSize(0);

  //now do the drawing in a stacked way - we need to plot reversed and also add always all the other histograms to let it look like it is stacked
  fill_stackedHistos(v_processes);
  //where is the signal?
  int sigIdx = -1;
  int VHbbIdx = -1;
  for(int i = v_processes.size()-1; i >= 0; i--){
    if( v_processes.at(i).name == "signal" ) sigIdx = i;
    if(v_processes.at(i).name == "VHbb") VHbbIdx = i;
  }
  if(sigIdx < 0) std::cout << "Couldn't find a signal to do the blinding. Please check." << std::endl;
  if(VHbbIdx < 0) std::cout << "Couldn't find a VHbb. Please check." << std::endl;
  
  // TODO: just a test for VHccBDT trafoD
  HistoTransform histoTrafo;
  histoTrafo.transformBkgBDTs = true;
  histoTrafo.trafoSixY = 5.;
  histoTrafo.trafoSixZ = 5.;
  int trafo_method = 6;
  float trafo_maxUnc = 1.;
  vector<int> trafo_bins = histoTrafo.getRebinBins(v_processes.at(v_processes.size()-1).hist_stack, v_processes.at(sigIdx).hist, trafo_method, trafo_maxUnc);

  TCanvas *c1 = new TCanvas(s_varName.c_str(),"",750,850);
  c1->cd();

  TPad *P_1 = new TPad((s_varName+"_1").c_str(), (s_varName+"_1").c_str(), 0, 0, 1, 0.3);
  TPad *P_2 = new TPad((s_varName+"_2").c_str(), (s_varName+"_2").c_str(), 0, 0.3, 1, 1);
  TPad *P_3 = new TPad((s_varName+"_3").c_str(), (s_varName+"_3").c_str(), 0, 0, 1, 0.3);
  P_1->Draw();
  P_2->Draw();
  P_1->SetBottomMargin(0.3);
  P_1->SetTopMargin(0.02);
  P_1->SetRightMargin(0.05);
  P_1->SetGridy();
  P_2->SetTopMargin(0.05);
  P_2->SetBottomMargin(0.02);
  P_2->SetRightMargin(0.05);

  P_3->SetFillStyle(4000);
  P_3->SetFrameFillStyle(0);
  P_3->SetTopMargin(0.02);
  P_3->SetBottomMargin(0.3);
  P_3->SetRightMargin(0.05);

  P_2->cd();

  //set the max value to either data or the stacked backgrounds
  if(s_varName == "VHccBDT"){
    histoTrafo.rebinHisto(h_data,&trafo_bins,false,false);
  }
  double d_maximum = h_data->GetMaximum();
  if(d_maximum<v_processes.at(v_processes.size()-1).hist_stack->GetMaximum()){
    d_maximum = v_processes.at(v_processes.size()-1).hist_stack->GetMaximum();
  }

  //Blinding time \o/
  //unblinded signal
  TH1D* h_unblind = (TH1D*) h_data->Clone( (std::string(h_data->GetName())+"_unblind").c_str() );
  //fill the blinding map
  for(int i = 1; i < h_data->GetNbinsX()+1; i++){
    bool blind = false;
    if(s_varName == "mBB" || s_varName.find("mBBJ") != std::string::npos ){ //mBB between 70 and 140
      if( (h_data->GetXaxis()->GetBinCenter(i) + h_data->GetXaxis()->GetBinWidth(i)*0.5 > 70) && (h_data->GetXaxis()->GetBinCenter(i) - h_data->GetXaxis()->GetBinWidth(i)*0.5 < 140) ) blind = true;
    }
    if(blind){
      h_data->SetBinContent(i, 0.0);
      h_data->SetBinError(i,0.0);
    }
  }

  h_data->GetYaxis()->SetRangeUser(0,d_maximum*1.70);
  
  //h_data->Draw("e");
  
  
  
  //now to the reverse drawing to make the plot look stacked
  for(int i = v_processes.size()-1; i >= 0; i--){
    if(s_varName == "VHccBDT"){
        histoTrafo.rebinHisto(v_processes.at(i).hist_stack,&trafo_bins, false, false);
        histoTrafo.rebinHisto(v_processes.at(i).hist, &trafo_bins, false, false);
        histoTrafo.rebinHisto(v_processes.at(i).h_rat_sys_up, &trafo_bins, false, false);
        histoTrafo.rebinHisto(v_processes.at(i).h_rat_sys_down, &trafo_bins, false, false);
        histoTrafo.rebinHisto(v_processes.at(i).h_rat_sys_modelling_up, &trafo_bins, false, false);
        histoTrafo.rebinHisto(v_processes.at(i).h_rat_sys_modelling_down, &trafo_bins, false, false);
    }

    if(s_varName == "ptv_all") v_processes.at(i).hist_stack->GetYaxis()->SetRangeUser(0.1,d_maximum*20);
    else v_processes.at(i).hist_stack->GetYaxis()->SetRangeUser(0,d_maximum*1.70);  
    
    std::cout << "Check Bin Error: Stack " << i << std::endl;
    for(int j = 1; j < v_processes.at(i).hist_stack->GetNbinsX()+1; j++){
        double temperr = v_processes.at(i).hist_stack->GetBinError(j);
        double tempcont = v_processes.at(i).hist_stack->GetBinContent(j);
        double tempcenter = v_processes.at(i).hist_stack->GetBinCenter(j);
        std::cout << "Bin " << j << ", Center: " << tempcenter << ", Content: "<< tempcont << ", Error: " << temperr << std::endl;
    }
    v_processes.at(i).hist_stack->Draw("hist f same");
  }
  // New blinding time for VHccBDT!
  float temp_sigyields = v_processes.at(sigIdx).hist->Integral();
  if(s_varName == "VHccBDT"){
    for(int i = 1; i < h_data->GetNbinsX()+1; i++){
        bool blind_VHccBDT = false;
        if(v_processes.at(sigIdx).hist->Integral(1, i)>=0.5*temp_sigyields)blind_VHccBDT = true;
        if(blind_VHccBDT){
            h_data->SetBinContent(i, 0.0);
            h_data->SetBinError(i, 0.0);
        }     
    }
  }

  // additionally draw signal & VHbb
  v_processes.at(sigIdx).hist->Scale(300);
  v_processes.at(sigIdx).hist->SetLineColor(v_processes.at(sigIdx).color);
  v_processes.at(sigIdx).hist->SetFillStyle(0);
  v_processes.at(sigIdx).hist->SetMarkerSize(0);
  v_processes.at(sigIdx).hist->SetLineWidth(2);
  v_processes.at(sigIdx).hist->Draw("hist same");

  v_processes.at(VHbbIdx).hist->Scale(5);
  v_processes.at(VHbbIdx).hist->SetLineColor(v_processes.at(VHbbIdx).color);
  v_processes.at(VHbbIdx).hist->SetFillStyle(0);
  v_processes.at(VHbbIdx).hist->SetMarkerSize(0);
  v_processes.at(VHbbIdx).hist->SetLineWidth(2);
  v_processes.at(VHbbIdx).hist->Draw("hist same");

  
  
  //Calculate chi2 and KS values to judge compatibility
  double res[h_unblind->GetNbinsX()];
  double chi2 = h_unblind->Chi2Test(v_processes.at(v_processes.size()-1).hist_stack, "UW");
  //double chi2 = h1->Chi2Test(v_processes.at(v_processes.size()-1).hist_stack, "WW", res);
  double KS = h_unblind->KolmogorovTest(v_processes.at(v_processes.size()-1).hist_stack);

  //BM: Fill histograms
  if(h_chi2P != NULL) h_chi2P->Fill(chi2);
  if(h_ksP != NULL) h_ksP->Fill(KS);


  P_2->RedrawAxis();

  std::string data_legend = "";
  if(year=="2015_2016") data_legend = "#sqrt{s} = 13 TeV, 36.2 fb^{-1} (2015/6)"; //2015+2016
  else if(year=="2017") data_legend = "#sqrt{s} = 13 TeV, 44.3 fb^{-1} (2017)";   //2017
  else if(year=="2018") data_legend = "#sqrt{s} = 13 TeV, 58.5 fb^{-1} (2018)";   //2018
  else if(year=="all") data_legend = "#sqrt{s} = 13 TeV, 139 fb^{-1} (Run 2)";    //all Run 2 data
  //Legend with all backgrounds
  TLatex l_LS;
  l_LS.SetNDC();
  l_LS.SetTextColor(1);
  l_LS.SetTextSize(0.045);
  //ATLASLabel(0.2,0.86,"Internal",1)
  l_LS.DrawLatex(0.16, 0.87, "#font[72]{ATLAS} #font[42]{Internal}");
  l_LS.SetTextSize(0.035);
  l_LS.DrawLatex(0.16, 0.81, data_legend.c_str());

  //std::string add_message = "";
  if(strcmp(s_region.c_str(),"R2")==0){
    l_LS.DrawLatex(0.16, 0.76, "p_{T}^{V} #geq 150 GeV, 2 jets");
  } else if(strcmp(s_region.c_str(),"R3")==0){
    l_LS.DrawLatex(0.16, 0.76, "p_{T}^{V} #geq 150 GeV, 3 jets");
  } else if(strcmp(s_region.c_str(),"R4")==0){
    l_LS.DrawLatex(0.16, 0.76, "p_{T}^{V} #geq 150 GeV, 4jets");
  }else if(strcmp(s_region.c_str(),"CR2")==0){
    l_LS.DrawLatex(0.16, 0.76, "top CR, 3 jets");
  }

  string Significance_filename = out + "/Significance_all_SR_" + s_varName + ".txt";
  string SoverB_filename = out + "/SoverB_all_SR_" + s_varName + ".txt";
  string SoverB_blind_filename = out + "/SoverB_all_blind_SR_" + s_varName + ".txt";
  string signal_yields_filename = out + "/SigYields_all_SR_" + s_varName + ".txt";
  string background_yields_filename = out + "/BkgYields_all_SR_" + s_varName + ".txt";
  std::stringstream stream_ctag;
  std::stringstream stream_SoverB;
  if(s_varName == "mBB" || s_varName == "mBB_AllSig" || s_varName == "VHccBDT"){
    v_processes.at(sigIdx).hist->Scale(1/300.);
    double Significance = ComputeSignificance(v_processes, sigIdx);
    double SoverB = ComputeSigOverBkg(v_processes, sigIdx);
    double SoverB_blind = ComputeSigOverBkg_blind(v_processes, sigIdx);
    double signal_yields = v_processes.at(sigIdx).hist->Integral();
    double bkg_yields = v_processes.at(v_processes.size()-1).hist_stack->Integral()-signal_yields;
    std::ofstream log_Significance(Significance_filename, std::ios_base::app | std::ios_base::out);
    log_Significance << s_region << "\t" << s_pTV << "\t" << s_ctag << "\t" << Significance <<"\n";
    log_Significance.close();
    std::ofstream log_SoverB(SoverB_filename, std::ios_base::app | std::ios_base::out);
    log_SoverB << s_region << "\t" << s_pTV << "\t" << s_ctag << "\t" << SoverB <<"\n";
    log_SoverB.close();
    std::ofstream log_SoverB_blind(SoverB_blind_filename, std::ios_base::app | std::ios_base::out);
    log_SoverB_blind << s_region << "\t" << s_pTV << "\t" << s_ctag << "\t" << SoverB_blind << "\n";
    log_SoverB_blind.close();
    std::ofstream log_SigYields(signal_yields_filename, std::ios_base::app | std::ios_base::out);
    log_SigYields << s_region << "\t" << s_pTV << "\t" << s_ctag << "\t" << signal_yields <<"\n";
    log_SigYields.close();
    std::ofstream log_BkgYields(background_yields_filename, std::ios_base::app | std::ios_base::out);
    log_BkgYields << s_region << "\t" << s_pTV << "\t" << s_ctag << "\t" << bkg_yields <<"\n";
    log_BkgYields.close();
    v_processes.at(sigIdx).hist->Scale(300);

    stream_ctag << s_ctag << ", " << s_pTV << ", " << std::setprecision(4) << Significance;
    if(s_varName == "VHccBDT")stream_SoverB << "(300 #times) S/B = " << std::setprecision(4) << SoverB*300;
    else stream_SoverB << "(300 #times) S/B = " << std::setprecision(4) << SoverB*300 << ", (blinded) " << std::setprecision(4) << SoverB_blind*300;
  }

  if(strcmp(s_ctag.c_str(),"default")==0){
    l_LS.DrawLatex(0.16, 0.66, "No c-tag requirement");
  }else if(strcmp(s_ctag.c_str(),"1Ctag")==0){
    l_LS.DrawLatex(0.16, 0.66, "1 c-tag");
  }else{
    if(s_varName == "mBB" || s_varName == "mBB_AllSig" || s_varName == "VHccBDT"){
      l_LS.DrawLatex(0.16, 0.66, stream_ctag.str().c_str());
      l_LS.DrawLatex(0.16, 0.61, stream_SoverB.str().c_str());
    }else{
      stream_ctag << s_ctag << ", " << s_pTV;
      l_LS.DrawLatex(0.16, 0.66, stream_ctag.str().c_str());
    }
  }
  
  
  

  std::stringstream stream; 
  stream << "#chi^{2}: " << std::setprecision(3) << chi2 << "  KS: " << std::setprecision(3) << KS;
  l_LS.DrawLatex(0.16, 0.71, stream.str().c_str());

  //Ratio Pad 
  P_1->cd();
  P_1->SetGridy();

  TH1F* h_ratio = (TH1F*) getDataMCComparison(h_data,v_processes.at(v_processes.size()-1).hist_stack);
  
  h_ratio->GetXaxis()->SetTitleSize(0.1);
  h_ratio->GetXaxis()->SetTitleOffset(1.35);
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetLabelOffset(0.03);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetTitleOffset(0.55);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetRangeUser(0.4, 1.6);
  h_ratio->GetYaxis()->SetTitle("data/pred.");
  h_ratio->GetXaxis()->SetTitle(s_legX.c_str());

  h_ratio->SetLineColor(kBlack);
  h_ratio->SetLineWidth(2);
  h_ratio->SetMarkerStyle(8);
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerSize(1.0);

  h_ratio->GetYaxis()->SetNdivisions(3,6,0,kTRUE);
  h_ratio->GetYaxis()->SetLabelOffset(999);  
  
  h_ratio->Draw("ep2");

  
  //draw the statistical uncertainty band of the MC thingy in the ratio plot
  //TODO: remove for high pTV
  TH1D* h_stat_up; 
  TH1D* h_stat_down;
  TH1D* h_rat_sys_up;
  TH1D* h_rat_sys_down;
  TH1D* h_rat_sys_model_up;
  TH1D* h_rat_sys_model_down;
  TGraphAsymmErrors* band_sys;
  TGraphAsymmErrors* band_sys_model;
  if(s_pTV!="400GeV+"){
  
    h_stat_up = (TH1D*) getRatioStatUncertaintyBand(v_processes.at(v_processes.size()-1).hist_stack, 0);
    h_stat_down = (TH1D*) getRatioStatUncertaintyBand(v_processes.at(v_processes.size()-1).hist_stack, 1);
  
    //draw systematic uncertainties
    // Remove for R2, 42, 46, 62, 63
  
    h_rat_sys_up = getTotalSystematics(v_processes, 0, h_stat_up);
    h_rat_sys_down = getTotalSystematics(v_processes, 1, h_stat_down);
    h_rat_sys_model_up = getTotalSystematics(v_processes, 2, h_rat_sys_up);
    h_rat_sys_model_down = getTotalSystematics(v_processes, 3, h_rat_sys_down);

    band_sys = getErrorBand(h_rat_sys_up, h_rat_sys_down);
    band_sys_model = getErrorBand(h_rat_sys_model_up, h_rat_sys_model_down);

    band_sys->SetLineColor(kWhite);
    band_sys_model->SetLineColor(kWhite);
    band_sys->SetFillColorAlpha(TColor::GetColor("#3598DB"), 1.0);  
    band_sys_model->SetFillColorAlpha(TColor::GetColor("#AED6F1"), 1.0);

    band_sys_model->Draw("e2 same");
    band_sys->Draw("e2 same");
 
    h_ratio->Draw("ep2 same");
  }

  //Ratio plot legend
  /*
  TLegend *LR = new TLegend(0.15, 0.33, 0.45, 0.48);
  LR->SetNColumns(2);
  LR->SetFillColor(0);
  LR->SetLineWidth(0);
  LR->SetTextSize(0.08);
  LR->SetBorderSize(1);
  LR->SetLineColor(kWhite);
  //LR->AddEntry(band_sys, "exp. sys.", "f");
  //LR->AddEntry(band_sys_model, "total sys.", "f");
  LR->Draw();*/

  P_1->RedrawAxis();

  c1->Modified();

  TLine *L_One;
  L_One = new TLine(h_data->GetXaxis()->GetXmin(), 1, h_data->GetXaxis()->GetXmax(), 1);

  //TLine *L_One = new TLine(80, 1, 280, 1);
  L_One->SetLineWidth(2);
  L_One->SetLineColor(1);
  L_One->SetLineStyle(2);
  L_One->Draw();

  TText t;
  t.SetTextAlign(32);
  t.SetTextSize(0.1);
  t.SetTextAngle(0);
  t.SetNDC();
  t.DrawText(0.114,0.93,"1.5");
  t.DrawText(0.114,0.65,"1.0");
  t.DrawText(0.114,0.35,"0.5");

  c1->cd();
  P_3->Draw();
  P_3->cd();
  TH1D *h_SigBkg;
  TH1D *h_SoverB;
  if(s_varName == "mBB" || s_varName == "mBB_AllSig" || s_varName == "VHccBDT"){
    v_processes.at(sigIdx).hist->Scale(1/300.);
    h_SigBkg = (TH1D *)HistSignificance(v_processes, sigIdx);
    h_SoverB = (TH1D *)HistSoverB(v_processes, sigIdx);
    v_processes.at(sigIdx).hist->Scale(300.);
    
    double tempSoverB = 0.0;
    string tempSoverB_filename = out + "/Significance_singlePlot_" + s_region + "_" + s_varName + "_" + s_pTV + "_" + s_ctag + ".txt";
    for(int bin = 1; bin<=h_SigBkg->GetNbinsX(); bin++){
        std::ofstream templog_SoverB(tempSoverB_filename, std::ios_base::app | std::ios_base::out);
        tempSoverB = h_SigBkg->GetBinContent(bin) /100.0;
        templog_SoverB << bin << "\t" << tempSoverB <<"\n";
        templog_SoverB.close();    
    }
    h_SigBkg->SetLineColor(kRed);
    h_SoverB->SetLineColor(kGreen+2);
    h_SigBkg->GetYaxis()->SetRangeUser(0., h_SigBkg->GetMaximum()*1.5);
    h_SigBkg->GetYaxis()->SetAxisColor(kRed);
    h_SigBkg->GetYaxis()->SetLabelColor(kRed);
    h_SigBkg->GetYaxis()->SetLabelSize(0.1);
    h_SigBkg->GetYaxis()->SetNdivisions(5,5,0, kTRUE);
    h_SigBkg->GetXaxis()->SetNdivisions(0,0,0, kTRUE);
  }

  //Ratio plot legend
  TLegend *LR = new TLegend(0.15, 0.80, 0.50, 0.95);
  LR->SetNColumns(2);
  LR->SetFillColor(0);
  LR->SetFillStyle(0);
  LR->SetLineWidth(0);
  LR->SetTextSize(0.08);
  LR->SetBorderSize(0);
  LR->SetLineColor(kWhite);
  if(s_varName != "mBB" && s_varName!="mBB_AllSig" && s_pTV!="400GeV+" && s_varName!="VHccBDT")LR->AddEntry(band_sys, "exp. sys.", "f");
  //if(s_varName != "mBB" && s_varName!="mBB_AllSig" && s_varName!="VHccBDT") LR->AddEntry(band_sys, "exp. sys.", "f");
  if (s_varName == "mBB" || s_varName=="mBB_AllSig" || s_varName == "VHccBDT"){
    LR->AddEntry(h_SigBkg, "100#times significance","l");
    LR->AddEntry(h_SoverB, "300#times s/b","l");
    h_SigBkg->Draw("Y+ same");
    h_SoverB->Draw("Y+ same");
  }
  LR->Draw("Same");
  c1->Modified();

  //Now draw the total uncertainty band on top of the top panel
  // remove for R2, 42, 46, 62, 63!!!
  P_2->cd();
  // TODO: remove for high pTV
  TGraphAsymmErrors *ge;
  if(s_pTV!="400GeV+"){
    ge = getTotalUncertaintyBand(v_processes.at(v_processes.size()-1).hist_stack, h_rat_sys_model_up, h_rat_sys_model_down);
    gStyle->SetHatchesSpacing(1.0);
    gStyle->SetHatchesLineWidth(2);
  
    ge->SetFillColor(kGray + 1);
    ge->SetFillStyle(3354);
    ge->SetLineColor(kGray + 1);
    ge->SetLineWidth(0);
    ge->GetXaxis()->SetLabelOffset(999);
    ge->GetXaxis()->SetLabelSize(0);
    ge->Draw("e2 same");
  }

  //draw the data again to be on top of everything
  h_data->GetYaxis()->SetRangeUser(0,d_maximum*1.70);
  h_data->Draw("e same");

  P_2->RedrawAxis();

  //  TLegend *L1 = new TLegend(0.510, 0.733, 0.735, 0.909);
  TLegend *L1 = new TLegend(0.510, 0.689, 0.735, 0.909);
  L1->SetFillColor(0);
  L1->SetLineWidth(0);
  L1->SetTextSize(0.035);
  L1->SetBorderSize(1);
  L1->SetLineColor(kWhite);

  TLegend *L2 = new TLegend(0.705, 0.689, 0.930, 0.909);
  L2->SetFillColor(0);
  L2->SetLineWidth(0);
  L2->SetTextSize(0.035);
  L2->SetBorderSize(1);
  L2->SetLineColor(kWhite);

  gStyle -> SetLegendFont(42);

  L1->AddEntry(h_data,"Data","lp");
  //L2->AddEntry(ge,"Uncertainty","f");

  unsigned int n_processes = v_processes.size();
  //v_processes.at(VHbbIdx).hist->SetFillStyle(1001);
  //v_processes.at(VHbbIdx).hist->SetFillColor(v_processes.at(VHbbIdx).color);
  for(unsigned int i = 0; i < n_processes; i++){
    if(i < (n_processes/2)){
      if(i==sigIdx){L1->AddEntry(v_processes.at(i).hist, "VH(cc) #times 300", "l");}
      //else{L1->AddEntry(v_processes.at(i).hist,v_processes.at(i).tex_name.c_str(),"f");}
      else if(i==VHbbIdx){
          //v_processes.at(i).hist->SetFillStyle(0);
          L2->AddEntry(v_processes.at(i).hist, "VH(bb) #times 5", "l");
      }
      else{L1->AddEntry(v_processes.at(i).hist,v_processes.at(i).tex_name.c_str(),"f");}
    }
    if(i >= (n_processes/2)){
      if(i==sigIdx){L1->AddEntry(v_processes.at(i).hist, "VH(cc) #times 300", "l");}
      //else {L2->AddEntry(v_processes.at(i).hist,v_processes.at(i).tex_name.c_str(),"f");}
      else if(i==VHbbIdx){
          //v_processes.at(i).hist->SetFillStyle(0);
          L2->AddEntry(v_processes.at(i).hist, "VH(bb) #times 5", "l");
      }
      else{L2->AddEntry(v_processes.at(i).hist,v_processes.at(i).tex_name.c_str(),"f");}
    }
  }

  L1->Draw();
  L2->Draw();


  c1->Modified();

  gSystem -> mkdir( out.c_str(), true ); 
  string s_saveName = out + "/";
  s_saveName.append(s_varName);
  s_saveName.append("/");
  gSystem -> mkdir(s_saveName.c_str(), true);
  s_saveName.append("C_");
  s_saveName.append(s_lepChannel);
  s_saveName.append("_");
  s_saveName.append(s_region);
  s_saveName.append("_");
  s_saveName.append(s_pTV);
  s_saveName.append("_");
  s_saveName.append(s_varName);
  s_saveName.append("_");
  s_saveName.append(s_ctag);
  s_saveName.append("_fullSelection.pdf");
  c1->SaveAs(s_saveName.c_str());

 /* //Save QQ plot
//  cqq->SaveAs("qqplot.pdf");
  string qq_saveName = out + "/qq_plot/QQ_";
  qq_saveName.append(s_region);
  qq_saveName.append("_");
  qq_saveName.append(s_varName);
  qq_saveName.append(".pdf");
  cqq->SaveAs(qq_saveName.c_str());

  //Save residual plot
  string res_saveName = out + "/res/res_";
  res_saveName.append(s_region);
  res_saveName.append("_");
  res_saveName.append(s_varName);
  res_saveName.append(".pdf");
  cres->SaveAs(res_saveName.c_str());

*/

  std::cout << "Saving yields as well..." << std::endl;
  //save_yields(h_data, v_processes, s_saveName);

  if(s_pTV!="")delete h_rat_sys_up, h_rat_sys_down, h_rat_sys_model_up, h_rat_sys_model_down, band_sys, band_sys_model;
  //delete h_rat_sys_up, h_rat_sys_down, h_rat_sys_model_up, h_rat_sys_model_down, band_sys, band_sys_model;
  delete h_data, h_SigBkg, h_unblind, h_ratio;
  delete f_input;
  

}

//////////////////////////////////////////
void fill_processHistos(vector<process> &v_processes, TFile *f, string s_region, string s_pTV, string s_lepChannel, string s_ctag, string s_varName, scaleFactorProvider* sfp){

  for(auto &pr: v_processes){
    for(auto &id: pr.v_dsids){
      std::vector<std::string> v_lepChannel;
      if(strcmp(s_lepChannel.c_str(), "emu") == 0){
	      v_lepChannel.push_back("el");
	      v_lepChannel.push_back("mu");
      }else{
	      v_lepChannel.push_back(s_lepChannel);
      }
      for(auto &lc: v_lepChannel){

	      string s_hName = "hist_" + id + "_" + lc + "_" + s_region + "_" + s_pTV + "_" + s_varName + "_" + s_ctag; 
	      if(id==-1){
	        s_hName = "hist_QCD_" + lc + "_" + s_region + "_" + s_pTV + "_" + s_varName + "_" + s_ctag; 
	      }
	      TH1D* htemp;
        
        htemp = (TH1D*) f->Get(s_hName.c_str());
	      if(htemp == NULL){
	        cout << "Histogram " << s_hName << " is not available!" << endl;
	        continue;
	      }
	      htemp->Clone();
	      htemp->SetDirectory(0);

  	    binsHist(htemp);
  
  	    //histogram rescaling needs to be done here
  	    if(sfp != NULL){
  	      double sf = sfp->getScaleFactor(id);
  	      cout << "Scaling histogram" << s_hName << " by " << sf << endl;
    	    htemp->Scale(sf);
    	  }

    	  if(pr.hist == NULL){
    	    pr.hist = htemp;
    	  }else{
    	    pr.hist->Add(htemp);
    	  }
      }
    }
    
    if(pr.hist == NULL){
      cout << "Process " << pr.name << " not available for " << s_region << ", " << s_lepChannel << ", " << s_varName << ", " << s_ctag << endl;
      continue;
    }

    //set the fill color and line color, only needed if plotted non-stacked
    pr.hist->SetLineColor(kBlack);
    pr.hist->SetLineWidth(2);
    pr.hist->SetFillColor(pr.color);
    pr.hist->SetLabelOffset(999);
    pr.hist->SetLabelSize(0);

    //now copy the nominal hist and make empty for the systematics
    TString tmp_hname = pr.hist->GetName();
    pr.h_rat_sys_up = (TH1D*) pr.hist->Clone((tmp_hname+"_sys_up").Data());
    pr.h_rat_sys_up->Reset();
    pr.h_rat_sys_down = (TH1D*) pr.hist->Clone((tmp_hname+"_sys_down").Data());
    pr.h_rat_sys_down->Reset();
    pr.h_rat_sys_modelling_up = (TH1D*) pr.hist->Clone((tmp_hname+"_sys_mod_up").Data());
    pr.h_rat_sys_modelling_up->Reset();
    pr.h_rat_sys_modelling_down = (TH1D*) pr.hist->Clone((tmp_hname+"_sys_mod_down").Data());
    pr.h_rat_sys_modelling_down->Reset();

    //fill the systematic histograms
    for(auto &id: pr.v_dsids){
      std::vector<std::string> v_lepChannel;
      if(strcmp(s_lepChannel.c_str(), "emu") == 0){
	      v_lepChannel.push_back("el");
	      v_lepChannel.push_back("mu");
      }else{
	      v_lepChannel.push_back(s_lepChannel);
      }
      for(auto &lc: v_lepChannel){

      	TString s_hName = "hist_" + id + "_" + lc + "_" + s_region + "_" + s_pTV + "_" + s_varName + "_" + s_ctag; 
      	if(id==-1){
      	  s_hName = "hist_QCD_" + lc + "_" + s_region + "_" + s_pTV + "_" + s_varName + "_" + s_ctag; 
      	}
      
      	fillSystematicHistograms(f, s_hName, pr.hist, pr.h_rat_sys_up, pr.h_rat_sys_down, pr.h_rat_sys_modelling_up, pr.h_rat_sys_modelling_down);

      }
    }

  }

}

//////////////////////////////////////////
void fill_stackedHistos(vector<process> &v_processes){
  for(int i = v_processes.size()-1; i>=0; i--){
    for(int j = 0; j <= i; j++){
      if(j==0){
	if(v_processes.at(j).hist == NULL){
	  cout << "The first process you add has to be non-empty!" << endl;
	  exit(1);
	}
	v_processes.at(i).hist_stack = (TH1D*) v_processes.at(j).hist->Clone();
      }else{
	if(v_processes.at(j).hist != NULL){
	  v_processes.at(i).hist_stack->Add(v_processes.at(j).hist);
	}
      }
    }
    v_processes.at(i).hist_stack->SetLineColor(kBlack);
    v_processes.at(i).hist_stack->SetLineWidth(0);
    v_processes.at(i).hist_stack->SetFillColor(v_processes.at(i).color);
    v_processes.at(i).hist_stack->SetLabelOffset(999);
    v_processes.at(i).hist_stack->SetLabelSize(0);
  }
}


//////////////////////////////////////////
void save_yields(TH1* hdat, vector<process> &v_processes, string name){
  //a) replace the .pdf with .txt
  size_t f = name.find(".pdf");
  name.replace(f,string(".pdf").length(),".txt");

  //b)create a file and save the yields into it
  ofstream ofile((name).c_str(), ios::out);
  if(!ofile) {  
    cout << "Error: could not open file: " << name << "for writing" << endl;
    return;
  }

  int nwid=14;
  ofile << setw(nwid) << "process" << setw(nwid) << "entries" << setw(nwid) << "integral" << setw(nwid) << "error" << setw(nwid) << "error/integ."  << endl;

  int entries = hdat->GetEntries();
  double error = 0.0;
  double integral = hdat->IntegralAndError(0, hdat->GetNbinsX()+1, error);
  double relError = TMath::IsNaN(error/integral) ? 0.0 : error/integral;
  
  ofile << endl;
  ofile << setw(nwid) << "data" << setw(nwid) << entries << setw(nwid) << integral << setw(nwid) << error << setw(nwid) << relError << endl;
  ofile << endl;
  
  int ensb(0);
  double er2sb(0.0), insb(0.0), resb(0.0);
  
  for(unsigned int i = 0; i< v_processes.size(); i++){
    entries = v_processes.at(i).hist->GetEntries();
    error = 0.0;
    integral = v_processes.at(i).hist->IntegralAndError(0, v_processes.at(i).hist->GetNbinsX()+1, error);
    relError = TMath::IsNaN(error/integral) ? 0.0 : error/integral;
    ofile << setw(nwid) << v_processes.at(i).tex_name << setw(nwid) << entries << setw(nwid) << integral << setw(nwid) << error << setw(nwid) << relError << endl;

    ensb += entries;
    insb += integral;
    er2sb += error*error;

  }//for all processes
 
  resb = sqrt(er2sb)/insb;

  ofile << endl;
  ofile << setw(nwid) << "sum of bkgs" << setw(nwid) << ensb << setw(nwid) << insb << setw(nwid) << sqrt(er2sb) << setw(nwid) << resb << endl;

  //c)close file
  ofile.close();
}

//////////////////////////////////////////
void addSystematics(TH1* hnom, TH1* halt, TH1* hup, TH1* hdown){

  int nBins = hnom->GetNbinsX();

  for (int i = 1; i< nBins+1; i++){
    if(halt->GetBinContent(i) < 1e-10) continue;
    if(hnom->GetBinContent(i) < 1e-10) continue;
    if(TMath::IsNaN(halt->GetBinContent(i))) continue; //filter out NaN
    float rat = (halt->GetBinContent(i)/hnom->GetBinContent(i)) - 1.0;
    if (fabs(rat) > 0.5){
      cout << "WARNING WARNING WARNING THERE IS A HUGE UNC" << fabs(rat) << endl;
      cout << halt->GetName() << endl;
    }

    if( rat < 0.0 ){
      float unc_old = hdown->GetBinContent(i);
      float unc_new = sqrt(unc_old*unc_old + rat*rat);
      hdown->SetBinContent(i, unc_new);
    }else{
      float unc_old = hup->GetBinContent(i);
      float unc_new = sqrt(unc_old*unc_old + rat*rat);
      hup->SetBinContent(i, unc_new);
    }

  }
  
}

//BM: TODO this is very ugly, make it nicer
//////////////////////////////////////////
void fillSystematicHistograms(TFile* f, TString name, TH1* hnom, TH1* hup, TH1* hdown, TH1* hMup, TH1* hMdown){

  if(!g_doSystematics) return;
  if(name.Contains("bkg")) return;

  name.Append("_Sys");

  TDirectory *dsys = f->GetDirectory("Systematics");
  TH1 *htmp;
  TKey *key;
  TIter next( dsys->GetListOfKeys() );
  while( (key = (TKey *) next()) ){
    string cname = key->GetClassName();
    TClass *cl = gROOT->GetClass(cname.c_str());
    if(! cl->InheritsFrom(TH1::Class()) ) continue;    
    htmp = (TH1*) key->ReadObj();
    TString hname = htmp->GetName();
    if (!hname.Contains(name)) continue;

    binsHist(htmp);

    if(hname.Contains("SysModelling")){
      addSystematics(hnom, htmp, hMup, hMdown);
    }else{
      addSystematics(hnom, htmp, hup, hdown);
    }
   
  }

}


//////////////////////////////////////////
TH1D* getTotalSystematics(vector<process> &v_processes, int ud, TH1* h_add){
  TH1D* htmp;
  if(ud == 0){
    htmp = (TH1D*) v_processes.at(0).h_rat_sys_up->Clone("h_tot_sys_up");
    htmp->Reset();
    for(int i = 1; i < htmp->GetNbinsX()+1; i++){
      float unc_tot = 0.0;
      float val_tot = 0.0;
      for(int j = 0; j < v_processes.size(); j++){
	float unc_tmp = v_processes.at(j).h_rat_sys_up->GetBinContent(i);
	float val_tmp = v_processes.at(j).hist->GetBinContent(i);
	unc_tot += unc_tmp*unc_tmp*val_tmp*val_tmp;
	val_tot += val_tmp;
      }
      htmp->SetBinContent(i, 1.0 + sqrt(unc_tot)/val_tot);
    }
  }else if(ud == 1){
    htmp = (TH1D*) v_processes.at(0).h_rat_sys_down->Clone("h_tot_sys_down");
    htmp->Reset();
    for(int i = 1; i < htmp->GetNbinsX()+1; i++){
      float unc_tot = 0.0;
      float val_tot = 0.0;
      for(int j = 0; j < v_processes.size(); j++){
	float unc_tmp = v_processes.at(j).h_rat_sys_down->GetBinContent(i);
	float val_tmp = v_processes.at(j).hist->GetBinContent(i);
	unc_tot += unc_tmp*unc_tmp*val_tmp*val_tmp;
	val_tot += val_tmp;
      }
      htmp->SetBinContent(i, 1.0 - 1.0*sqrt(unc_tot)/val_tot);
    }

  }else if(ud == 2){
    htmp = (TH1D*) v_processes.at(0).h_rat_sys_modelling_up->Clone("h_tot_sys_modelling_up");
    htmp->Reset();
    for(int i = 1; i < htmp->GetNbinsX()+1; i++){
      float unc_tot = 0.0;
      float val_tot = 0.0;
      for(int j = 0; j < v_processes.size(); j++){
	float unc_tmp = v_processes.at(j).h_rat_sys_modelling_up->GetBinContent(i);
	float val_tmp = v_processes.at(j).hist->GetBinContent(i);
	unc_tot += unc_tmp*unc_tmp*val_tmp*val_tmp;
	val_tot += val_tmp;
      }
      htmp->SetBinContent(i, 1.0 + sqrt(unc_tot)/val_tot);
    }

  }else if(ud == 3){
    htmp = (TH1D*) v_processes.at(0).h_rat_sys_modelling_down->Clone("h_tot_sys_modelling_down");
    htmp->Reset();
    for(int i = 1; i < htmp->GetNbinsX()+1; i++){
      float unc_tot = 0.0;
      float val_tot = 0.0;
      for(int j = 0; j < v_processes.size(); j++){
	float unc_tmp = v_processes.at(j).h_rat_sys_modelling_down->GetBinContent(i);
	float val_tmp = v_processes.at(j).hist->GetBinContent(i);
	unc_tot += unc_tmp*unc_tmp*val_tmp*val_tmp;
	val_tot += val_tmp;
      }
      htmp->SetBinContent(i, 1.0 - 1.0*sqrt(unc_tot)/val_tot);
    }
  }

  for(int i = 1; i < htmp->GetNbinsX()+1; i++){
    float nbc = htmp->GetBinContent(i);
    if (nbc != nbc) htmp->SetBinContent(i, 1.0);
  }

  if(h_add == nullptr){
    return htmp;
  }else{
    for(int i = 1; i < htmp->GetNbinsX()+1; i++){
      float nbc = (htmp->GetBinContent(i)-1.0)*(htmp->GetBinContent(i)-1.0) + (h_add->GetBinContent(i)-1.0)* (h_add->GetBinContent(i)-1.0);
      if(nbc != nbc) nbc = 0.0;
      if(ud == 1 || ud == 3){
	htmp->SetBinContent(i, 1.0 - sqrt(nbc));
      }else{
	htmp->SetBinContent(i, 1.0 + sqrt(nbc));
      }
    }
    return htmp;
  }

}


//////////////////////////////////////////
TGraphAsymmErrors* getErrorBand(TH1* hup, TH1* hdown){

  int nBins = hup->GetNbinsX();

  TGraphAsymmErrors* tgae = new TGraphAsymmErrors(nBins);
  TString name = hup->GetName();
  name.Append("_tgae");
  tgae->SetName(name);
  for(int i = 0; i < nBins; i++){
    float exl = hup->GetBinWidth(i+1)/2.0;
    float x = hup->GetBinCenter(i+1);
    float exh = hup->GetBinWidth(i+1)/2.0;
    float eyl = 1.0 - hdown->GetBinContent(i+1);
    float eyh = hup->GetBinContent(i+1) - 1.0;
    tgae->SetPoint(i, x, 1.0);
    tgae->SetPointError(i, exl, exh, eyl, eyh);
  }

  return tgae;

}


//////////////////////////////////////////
TGraphAsymmErrors* getTotalUncertaintyBand(TH1D* htot, TH1D* hup, TH1D* hdown){

  int nBins = htot->GetNbinsX();

  TGraphAsymmErrors* tgaetot = new TGraphAsymmErrors(nBins);
  TString name = htot->GetName();
  name.Append("_totalUncBand");
  tgaetot->SetName(name);
  for(int i = 0; i < nBins; i++){
    float exl = hup->GetBinWidth(i+1)/2.0;
    float x = hup->GetBinCenter(i+1);
    float exh = hup->GetBinWidth(i+1)/2.0;
    float eyl = 1.0 - hdown->GetBinContent(i+1);
    float eyh = hup->GetBinContent(i+1) - 1.0;
    float y = htot->GetBinContent(i+1);
    tgaetot->SetPoint(i, x, y);
    tgaetot->SetPointError(i, exl, exh, eyl*y, eyh*y);
  }

  return tgaetot;

}

scaleFactorProvider::scaleFactorProvider(std::string fileName) {
  std::ifstream file;
  file.open(fileName.c_str());
  if (!file.good()) {
    Error("scaleFactorProvider()", "Can't open file '%s'.", fileName.c_str());
    exit(EXIT_FAILURE);
  }

  Info("scaleFactorProvider()", "Reading file '%s'.", fileName.c_str());

  // process file
  while (!file.eof()) {
    // read line
    std::string lineString;
    getline(file, lineString);
    //std::cout << lineString << std::endl;
    
    // skip empty lines
    if (lineString.length() == 0) {
      continue;
    }

    // skip lines starting with #
    if (lineString.find("#") == 0) {
      continue;
    }

    // store in map
    std::stringstream line(lineString);
    double sf;
    std::string sampleName;
    line >> sampleName >> sf;
    m_scaleFactors[sampleName] = sf;
  }
  file.close();
}


double scaleFactorProvider::getScaleFactor(std::string sampleName){
  
  if (m_scaleFactors.count(sampleName) == 0) {
    Error("scaleFactorProvider::getScaleFactor()", "Unknown sampleName %s", sampleName.c_str());
    exit(EXIT_FAILURE);
  }

  return m_scaleFactors[sampleName];
}
