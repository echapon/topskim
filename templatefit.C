#include "TFractionFitter.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"

const float dphicut = 2.;

TH1F *cat(TH1F *h1, TH1F *h2); // concatenate 2 histos

void templatefit() {
   TFile *fdata = TFile::Open("emuSkim_data.root");
   TFile *fttbar = TFile::Open("emuSkim_ttbarMC.root");
   TTree *tdata = (TTree*) fdata->Get("skimTree");
   TTree *tmc = (TTree*) fttbar->Get("skimTree");
   TCanvas *c1 = new TCanvas();

   // make the histos
   tdata->Draw("discr_csvV1[0]>>hOS_hdphi(10,0,1)",Form("nJt>0&&emuPt>0&&lepPt[0]+lepPt[1]>50&&emuSgn==0&&emuDphi>%f",dphicut));
   TH1F *hOS_hdphi = (TH1F*) gDirectory->Get("hOS_hdphi");
   hOS_hdphi->GetXaxis()->SetTitle("csvV1 (leading jet)");
   hOS_hdphi->GetYaxis()->SetTitle("Entries");
   tdata->Draw("discr_csvV1[0]>>hSS_hdphi(10,0,1)",Form("nJt>0&&emuPt>0&&lepPt[0]+lepPt[1]>50&&emuSgn!=0&&emuDphi>%f",dphicut));
   TH1F *hSS_hdphi = (TH1F*) gDirectory->Get("hSS_hdphi");
   hSS_hdphi->GetXaxis()->SetTitle("csvV1 (leading jet)");
   hSS_hdphi->GetYaxis()->SetTitle("Entries");
   tdata->Draw("discr_csvV1[0]>>hOS_ldphi(10,0,1)",Form("nJt>0&&emuPt>0&&lepPt[0]+lepPt[1]>50&&emuSgn==0&&emuDphi<%f",dphicut));
   TH1F *hOS_ldphi = (TH1F*) gDirectory->Get("hOS_ldphi");
   hOS_ldphi->GetXaxis()->SetTitle("csvV1 (leading jet)");
   hOS_ldphi->GetYaxis()->SetTitle("Entries");
   tdata->Draw("discr_csvV1[0]>>hSS_ldphi(10,0,1)",Form("nJt>0&&emuPt>0&&lepPt[0]+lepPt[1]>50&&emuSgn!=0&&emuDphi<%f",dphicut));
   TH1F *hSS_ldphi = (TH1F*) gDirectory->Get("hSS_ldphi");
   hSS_ldphi->GetXaxis()->SetTitle("csvV1 (leading jet)");
   hSS_ldphi->GetYaxis()->SetTitle("Entries");
   tmc->Draw("discr_csvV1[0]>>hOS_ldphi_mc(10,0,1)",Form("nJt>0&&emuPt>0&&lepPt[0]+lepPt[1]>50&&emuSgn==0&&emuDphi<%f",dphicut));
   TH1F *hOS_ldphi_mc = (TH1F*) gDirectory->Get("hOS_ldphi_mc");
   hOS_ldphi_mc->GetXaxis()->SetTitle("csvV1 (leading jet)");
   hOS_ldphi_mc->GetYaxis()->SetTitle("Entries");
   tmc->Draw("discr_csvV1[0]>>hOS_hdphi_mc(10,0,1)",Form("nJt>0&&emuPt>0&&lepPt[0]+lepPt[1]>50&&emuSgn==0&&emuDphi>%f",dphicut));
   TH1F *hOS_hdphi_mc = (TH1F*) gDirectory->Get("hOS_hdphi_mc");
   hOS_hdphi_mc->GetXaxis()->SetTitle("csvV1 (leading jet)");
   hOS_hdphi_mc->GetYaxis()->SetTitle("Entries");
   TH1F *hOS_alldphi = cat(hOS_ldphi,hOS_hdphi);
   TH1F *hSS_alldphi = cat(hSS_ldphi,hOS_hdphi);
   TH1F *hOS_alldphi_mc = cat(hOS_ldphi_mc,hOS_hdphi_mc);

   // // get the normalisation from high dphi
   // float nOS_hdphi = hOS_hdphi->GetEntries();
   // float nSS_hdphi = hSS_hdphi->GetEntries();
   // float norm = nOS_hdphi/nSS_hdphi;
   // cout << "norm: " << nOS_hdphi << "/" << nSS_hdphi << " = " << norm << endl;
   // float dnorm = norm*sqrt(1./nOS_hdphi+1./nSS_hdphi);

   // hSS_hdphi->Scale(norm);
   // hSS_hdphi->SetLineColor(kRed);
   // hSS_hdphi->SetMarkerColor(kRed);
   // hSS_hdphi->Draw();
   // hOS_hdphi->Draw("same");
   // c1->SaveAs("plot_hdphi.pdf");
   // cout << "Delta phi > " << dphicut << ": OS vs SS" << endl;
   // cout << "Chi2/ndf: " << hSS_hdphi->Chi2Test(hOS_hdphi,"NORM UU CHI2") << "/"  << hSS_hdphi->GetNbinsX() << " (pvalue: " << hSS_hdphi->Chi2Test(hOS_hdphi,"NORM UU") << ")" << endl;
   // cout << "KS: " << hSS_hdphi->KolmogorovTest(hOS_hdphi) << endl;

   // // go the signal region (low dphi)
   // hSS_ldphi->Scale(norm);
   // hSS_ldphi->SetLineColor(kRed);
   // hSS_ldphi->SetMarkerColor(kRed);
   // hSS_ldphi->Draw();
   // hOS_ldphi->Draw("same");
   // c1->SaveAs("plot_ldphi.pdf");
   // cout << "Delta phi < " << dphicut << ": OS vs SS" << endl;
   // cout << "Chi2/ndf: " << hSS_ldphi->Chi2Test(hOS_ldphi,"NORM UU CHI2") << "/"  << hSS_ldphi->GetNbinsX() << " (pvalue: " << hSS_ldphi->Chi2Test(hOS_ldphi,"NORM UU") << ")" << endl;
   // cout << "KS: " << hSS_ldphi->KolmogorovTest(hOS_ldphi) << endl;
   // float nSS_ldphi_10 = hSS_ldphi->GetBinContent(10);
   // float nSS_ldphi_10_err = nSS_ldphi_10*sqrt(pow(hSS_ldphi->GetBinError(10)/nSS_ldphi_10,2)+pow(dnorm/norm,2));
   // float nOS_ldphi_10 = hOS_ldphi->GetBinContent(10);
   // float nOS_ldphi_10_err = hOS_ldphi->GetBinError(10);
   // cout << "last bin: scaled SS = " << nSS_ldphi_10 << " +/- " << nSS_ldphi_10_err << ", OS = " << nOS_ldphi_10 << " +/- " << nOS_ldphi_10_err << endl;
   // float chi2_10 = pow(nSS_ldphi_10-nOS_ldphi_10,2) / (pow(nSS_ldphi_10_err,2)+pow(nOS_ldphi_10_err,2));
   // cout << "chi2: " << chi2_10 << " (pvalue: " << TMath::Prob(chi2_10,1) << ", zvalue = " << sqrt(2)*TMath::ErfInverse(1.-TMath::Prob(chi2_10,1)) << ")" << endl;

   // TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
   // mc->Add(hSS_ldphi);
   // mc->Add(hOS_ldphi_mc);
   // TFractionFitter* fit = new TFractionFitter(hOS_ldphi, mc); // initialise
   // fit->Constrain(0,0.0,1.0); 
   // fit->Constrain(1,0.0,1.0); 
   // Int_t status = fit->Fit();               // perform the fit
   // std::cout << "fit status: " << status << std::endl;
   // if (status == 0) {                       // check on fit status
   //    TH1F* result = (TH1F*) fit->GetPlot();
   //    hOS_ldphi->Draw("Ep");
   //    THStack *hs = new THStack("hs","");
   //    double frac0, dfrac0, frac1, dfrac1;
   //    fit->GetResult(0,frac0,dfrac0);
   //    fit->GetResult(1,frac1,dfrac1);
   //    hSS_ldphi->Scale(frac0*hOS_ldphi->Integral(1,10)/hSS_ldphi->Integral(1,10));
   //    hSS_ldphi->SetFillColor(kRed);
   //    hs->Add(hSS_ldphi);
   //    hOS_ldphi_mc->Scale(frac1*hOS_ldphi->Integral(1,10)/hOS_ldphi_mc->Integral(1,10));
   //    hOS_ldphi_mc->SetFillColor(kBlue);
   //    hs->Add(hOS_ldphi_mc);
   //    hs->Draw("hist same");

   //    hOS_ldphi->Draw("Ep same");
   //    result->Draw("same");
   //    c1->RedrawAxis();
   //    c1->SaveAs("plot_OS_ldphi_template.pdf");
   //    // cout << hOS_ldphi->GetBinContent(10) << " " << result->GetBinContent(10) << endl;

   //    cout << "Delta phi < " << dphicut << ": OS vs template fit" << endl;
   //    cout << "Chi2/ndf: " << hOS_ldphi->Chi2Test(result,"UW CHI2") << "/"  << hOS_ldphi->GetNbinsX() << " (pvalue: " << hOS_ldphi->Chi2Test(result,"UW") << ")" << endl;
   //    cout << "KS: " << hOS_ldphi->KolmogorovTest(result,"UW") << endl;
   //    double nsig, dnsig;
   //    nsig = hOS_ldphi_mc->IntegralAndError(1,10,dnsig);
   //    cout << "I found " << nsig << " +/- " << nsig*sqrt(pow(dnsig/nsig,2)+pow(dfrac1/frac1,2)) << " signal events" << endl;
   // }

   // and now, fit simultaneously the low and high dphi

   TObjArray *mc2 = new TObjArray(2);        // MC histograms are put in this array
   mc2->Add(hSS_alldphi);
   mc2->Add(hOS_alldphi_mc);
   TFractionFitter* fit2 = new TFractionFitter(hOS_alldphi, mc2); // initialise
   fit2->Constrain(0,0.0,1.0); 
   fit2->Constrain(1,0.0,1.0); 
   Int_t status2 = fit2->Fit();               // perform the fit
   std::cout << "fit status: " << status2 << std::endl;
   if (status2 == 0) {                       // check on fit status
      TH1F* result = (TH1F*) fit2->GetPlot();
      hOS_alldphi->Draw("Ep");
      THStack *hs = new THStack("hs2","");
      double frac0, dfrac0, frac1, dfrac1;
      fit2->GetResult(0,frac0,dfrac0);
      fit2->GetResult(1,frac1,dfrac1);
      hSS_alldphi->Scale(frac0*hOS_alldphi->Integral(1,20)/hSS_alldphi->Integral(1,20));
      hSS_alldphi->SetFillColor(kRed);
      hs->Add(hSS_alldphi);
      hOS_alldphi_mc->Scale(frac1*hOS_alldphi->Integral(1,20)/hOS_alldphi_mc->Integral(1,20));
      hOS_alldphi_mc->SetFillColor(kBlue);
      hs->Add(hOS_alldphi_mc);
      hs->Draw("hist same");

      hOS_alldphi->Draw("Ep same");
      result->Draw("same");
      c1->RedrawAxis();
      c1->SaveAs("plot_OS_alldphi_template.pdf");
      // cout << hOS_ldphi->GetBinContent(10) << " " << result->GetBinContent(10) << endl;

      cout << "Delta phi </> " << dphicut << ": OS vs template fit" << endl;
      cout << "Chi2/ndf: " << hOS_alldphi->Chi2Test(result,"UW CHI2") << "/"  << hOS_alldphi->GetNbinsX() << " (pvalue: " << hOS_alldphi->Chi2Test(result,"UW") << ")" << endl;
      cout << "KS: " << hOS_alldphi->KolmogorovTest(result,"UW") << endl;
      double nsig, dnsig;
      nsig = hOS_alldphi_mc->IntegralAndError(1,20,dnsig);
      cout << "I found " << nsig << " +/- " << nsig*sqrt(pow(dnsig/nsig,2)+pow(dfrac1/frac1,2)) << " signal events" << endl;
   }
}

TH1F *cat(TH1F *h1, TH1F *h2) {
   TH1F *ans = new TH1F(Form("%s_%s",h1->GetName(),h2->GetName()),"Concatenated histo",h1->GetNbinsX()+h2->GetNbinsX(),0,2);
   for (int i=1; i<=h1->GetNbinsX(); i++) {
      ans->SetBinContent(i,h1->GetBinContent(i));
      ans->SetBinError(i,h1->GetBinError(i));
   }
   for (int i=1; i<=h2->GetNbinsX(); i++) {
      ans->SetBinContent(i+h1->GetNbinsX(),h2->GetBinContent(i));
      ans->SetBinError(i+h1->GetNbinsX(),h2->GetBinError(i));
   }
   return ans;
}
