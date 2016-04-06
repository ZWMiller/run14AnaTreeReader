
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
  float lowpt[numPtBins]  = {0.2, 0.6, 1.0, 1.5, 
    2.0, 2.5, 3.5, 6.5};
  float highpt[numPtBins] = {0.6, 1.0, 1.5, 2.0,
    2.5, 3.5, 6.5, 10.};


  // Get Histos
  TH2F* eHadDelPhiPt = (TH2F*)f->Get("hHadgEDelPhiPt");
  TH2F* hadPtEPt = (TH2F*)f->Get("hHadgPtEPt");
  TH1F* ePt = (TH1F*)f->Get("hEPt");
  float trigCount[numPtBins];
  if(DEBUG) cout << "Get Hist." << endl;

  //Do Projection
  TH1D* eHadDelPhi[numPtBins];
  for(int ptbin=0; ptbin<numPtBins; ptbin++)
  {
    eHadDelPhi[ptbin] = (TH1D*)eHadDelPhiPt->ProjectionY(Form("eHadDelPhi_%i",ptbin),eHadDelPhiPt->GetXaxis()->FindBin(lowpt[ptbin]),eHadDelPhiPt->GetXaxis()->FindBin(highpt[ptbin])-1);
    trigCount[ptbin] = ePt->Integral(ePt->GetXaxis()->FindBin(lowpt[ptbin]),ePt->GetXaxis()->FindBin(highpt[ptbin])-1);
  }
  if(DEBUG) cout << "Proj done." << endl;

  //create canvas
  const int numCanvas = numPtBins/9 + 1;
  TCanvas* dPhiPt[numCanvas];
  for(int q=0; q<numCanvas; q++)
  {
    dPhiPt[q] = new TCanvas(Form("dPhiPt_%i",q),"pT Dependence of DelPhi",50,50,1050,1050);
    dPhiPt[q]->Divide(3,3);
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
  }
  temp = ptCompare;
  temp->Print(name);

  sprintf(name, "%s.pdf]", fileName);
  temp->Print(name);
}

