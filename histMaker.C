
bool DEBUG = kTRUE;

void histMaker(const char* fileName="test")
{
  //Get TFile
  TFile* f = new TFile(fileName, "READ");
  if(!f->IsOpen())
  { std::cout << "!!! File Not Found !!!" << std::endl;
    exit(1); }  
  else
    std::cout<< "Opened " << fileName << std::endl;

  //For now, just define binning here
  const int numPtBins = 8;
  float lowpt[numPtBins]  = {1.3, 2.0, 3.0, 4.0, 
                             5.0, 6.0, 8.0, 10.0};
  float highpt[numPtBins] = {2.0, 3.0, 4.0, 5.0,
                             6.0, 8.0, 10., 20. };


  // Get Histos
  TH2F* eHadDelPhiPt    = (TH2F*)f->Get("hHadEDelPhiPt");
  TH2F* eHadDelPhiPt_US = (TH2F*)f->Get("hHadEEDelPhiPt_US");
  TH2F* eHadDelPhiPt_LS = (TH2F*)f->Get("hHadEEDelPhiPt_LS");
  TH2F* hadPtEPt = (TH2F*)f->Get("hHadPtEPt");
  TH1F* ePt = (TH1F*)f->Get("hEPt");
  TH1F* eePt_US = (TH1F*)f->Get("hEEPt_US");
  TH1F* eePt_LS = (TH1F*)f->Get("hEEPt_LS");
  float trigCount[numPtBins],trigCount_US[numPtBins],trigCount_LS[numPtBins];
  if(DEBUG) cout << "Get Hist." << endl;

  //Do Projection
  TH1D* eHadDelPhi[numPtBins];
  TH1D* eHadDelPhi_US[numPtBins];
  TH1D* eHadDelPhi_LS[numPtBins];
  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    eHadDelPhi[ptbin] = (TH1D*)eHadDelPhiPt->ProjectionY(Form("eHadDelPhi_%i",ptbin),eHadDelPhiPt->GetXaxis()->FindBin(lowpt[ptbin]),eHadDelPhiPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    eHadDelPhi_US[ptbin] = (TH1D*)eHadDelPhiPt_US->ProjectionY(Form("eHadDelPhi_US_%i",ptbin),eHadDelPhiPt_US->GetXaxis()->FindBin(lowpt[ptbin]),eHadDelPhiPt_US->GetXaxis()->FindBin(highpt[ptbin])-1);
    eHadDelPhi_LS[ptbin] = (TH1D*)eHadDelPhiPt_LS->ProjectionY(Form("eHadDelPhi_LS_%i",ptbin),eHadDelPhiPt_LS->GetXaxis()->FindBin(lowpt[ptbin]),eHadDelPhiPt_LS->GetXaxis()->FindBin(highpt[ptbin])-1);
    trigCount[ptbin] = ePt->Integral(ePt->GetXaxis()->FindBin(lowpt[ptbin]),ePt->GetXaxis()->FindBin(highpt[ptbin])-1);
    trigCount_US[ptbin] = eePt_US->Integral(eePt_US->GetXaxis()->FindBin(lowpt[ptbin]),eePt_US->GetXaxis()->FindBin(highpt[ptbin])-1);
    trigCount_LS[ptbin] = eePt_LS->Integral(eePt_LS->GetXaxis()->FindBin(lowpt[ptbin]),eePt_LS->GetXaxis()->FindBin(highpt[ptbin])-1);
  }
  if(DEBUG) cout << "Proj done." << endl;

  //create canvas
  const int numCanvas = numPtBins/9 + 1;
  TCanvas* dPhiPt[numCanvas];
  TCanvas* dPhiPt_US[numCanvas];
  TCanvas* dPhiPt_LS[numCanvas];
  for(int q=0; q<numCanvas; q++)
  {
    dPhiPt[q] = new TCanvas(Form("dPhiPt_%i",q),"pT Dependence of DelPhi",50,50,1050,1050);
    dPhiPt[q]->Divide(3,3);
    dPhiPt_US[q] = new TCanvas(Form("dPhiPt_US_%i",q),"pT Dependence of DelPhi US Pairs",50,50,1050,1050);
    dPhiPt_US[q]->Divide(3,3);
    dPhiPt_LS[q] = new TCanvas(Form("dPhiPt_LS_%i",q),"pT Dependence of DelPhi LS Pairs",50,50,1050,1050);
    dPhiPt_LS[q]->Divide(3,3);
  }
  TCanvas* ptCompare = new TCanvas("ptCompare","pT Comparison",50,50,1050,1050);
  if(DEBUG) cout << "Canvas made." << endl;

  // Draw on each pad
  TPaveText* lbl[numPtBins];
  char textLabel[100];
  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    if(DEBUG) cout << "pt loop: " << ptbin <<endl;
    lbl[ptbin] = new TPaveText(.67,.25,.85,.3,Form("NB NDC%i",ptbin));
    sprintf(textLabel,"%.2f < P_{T,e} < %.2f",lowpt[ptbin],highpt[ptbin]);
    lbl[ptbin]->AddText(textLabel);
    lbl[ptbin]->SetFillColor(kWhite);

    int activeCanvas = (int) ptbin/9;
    int activeBin = ptbin - activeCanvas*9; 
    dPhiPt[activeCanvas]->cd(activeBin+1);
    eHadDelPhi[ptbin]->Rebin(4);
    float binWidth = eHadDelPhi[ptbin]->GetXaxis()->GetBinWidth(10);
    eHadDelPhi[ptbin]->Scale(1./trigCount[ptbin]/binWidth);
    eHadDelPhi[ptbin]->GetXaxis()->SetRangeUser(-1.5,4.4);
    eHadDelPhi[ptbin]->SetMarkerStyle(20);
    eHadDelPhi[ptbin]->SetMarkerColor(kBlack);
    eHadDelPhi[ptbin]->SetLineColor(kBlack);
    eHadDelPhi[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    dPhiPt_US[activeCanvas]->cd(activeBin+1);
    eHadDelPhi_US[ptbin]->Rebin(4);
    binWidth = eHadDelPhi_US[ptbin]->GetXaxis()->GetBinWidth(10);
    eHadDelPhi_US[ptbin]->Scale(1./trigCount_US[ptbin]/binWidth);
    eHadDelPhi_US[ptbin]->GetXaxis()->SetRangeUser(-1.5,4.4);
    eHadDelPhi_US[ptbin]->SetMarkerStyle(20);
    eHadDelPhi_US[ptbin]->SetMarkerColor(kBlack);
    eHadDelPhi_US[ptbin]->SetLineColor(kBlack);
    eHadDelPhi_US[ptbin]->Draw();
    lbl[ptbin]->Draw("same");

    dPhiPt_LS[activeCanvas]->cd(activeBin+1);
    eHadDelPhi_LS[ptbin]->Rebin(4);
    binWidth = eHadDelPhi_LS[ptbin]->GetXaxis()->GetBinWidth(10);
    eHadDelPhi_LS[ptbin]->Scale(1./trigCount_LS[ptbin]/binWidth);
    eHadDelPhi_LS[ptbin]->GetXaxis()->SetRangeUser(-1.5,4.4);
    eHadDelPhi_LS[ptbin]->SetMarkerStyle(20);
    eHadDelPhi_LS[ptbin]->SetMarkerColor(kBlack);
    eHadDelPhi_LS[ptbin]->SetLineColor(kBlack);
    eHadDelPhi_LS[ptbin]->Draw();
    lbl[ptbin]->Draw("same");
  }  

  ptCompare->cd();
  gPad->SetLogz(1);
  hadPtEPt->Draw("colz");

  //Set front page
  TCanvas* fp = new TCanvas("fp","Front Page",100,0,1000,900);
  fp->cd();
  TBox *bLabel = new TBox(0.01, 0.88, 0.99, 0.99);
  bLabel->SetFillColor(38);
  bLabel->Draw();
  TLatex tl;
  tl.SetNDC();
  tl.SetTextColor(kWhite);
  tl.SetTextSize(0.033);
  char tlName[100];
  char tlName2[100];

  TString titlename = fileName;
  int found = titlename.Last('/');
  if(found >= 0){
    titlename.Replace(0, found+1, "");
  } 
  sprintf(tlName, "RUN 14 AuAu 200 GeV NPE-Hadron");
  tl.SetTextSize(0.05);
  tl.SetTextColor(kWhite);
  tl.DrawLatex(0.05, 0.92,tlName);

  TBox *bFoot = new TBox(0.01, 0.01, 0.99, 0.12);
  bFoot->SetFillColor(38);
  bFoot->Draw();
  tl.SetTextColor(kWhite);
  tl.SetTextSize(0.05);
  tl.DrawLatex(0.05, 0.05, (new TDatime())->AsString());
  tl.SetTextColor(kBlack);
  tl.SetTextSize(0.03);
  tl.DrawLatex(0.1, 0.14, titlename);

  // Place canvases in order
  TCanvas* temp = new TCanvas();
  char name[100];
  sprintf(name, "%s.pdf[", fileName);
  temp->Print(name);
  sprintf(name, "%s.pdf", fileName);
  temp = fp; // print front page
  temp->Print(name);
  for(int q=0; q<numCanvas; q++)
  {
    temp = dPhiPt[q]; // print data canvases
    temp->Print(name);
    temp = dPhiPt_US[q]; // print data canvases
    temp->Print(name);
    temp = dPhiPt_LS[q]; // print data canvases
    temp->Print(name);
  }
  temp = ptCompare;
  temp->Print(name);

  sprintf(name, "%s.pdf]", fileName);
  temp->Print(name);
}

