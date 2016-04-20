#include "TMinuit.h"
#include "TVirtualFitter.h"

#define NPLOT 6
#define EID 2
#define NFRAME 2
#define TRIG 4
#define NFILE 41
#define NPT 6
#define NUNCERTAINTY 10

TString TrigName[4] = {"MB","HT0","HT1","HT2"};
TString VarName[2] = {"theta","phi"};
TString FrameName[2] = {"Helicity","CS"};
TString Frame[2] ={"HX","CS"};
Float_t z1[10],z2[10],errorz1[10],errorz2[10],x[10],y[10];
Double_t chisquare,chisquare1,chisquare2;

Double_t func1(Double_t *x,Double_t *par)
{
    Double_t value1=par[0]*(1+par[1]*x[0]*x[0]);
    return value1;
}

Double_t func2(Double_t *y,Double_t *par)
{
    Double_t value2=par[3]*(1+2*par[2]*TMath::Cos(2*y[0])/(3+par[1]));
    return value2;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    const Int_t nbins = 10;
    Int_t i;
    Double_t xx,yy;
    Double_t chisq=0,chisq1=0,chisq2=0;
    Double_t delta[2];
    Int_t nbin=0;
    for(i=0;i<nbins;i++){
        xx = x[i];
        yy = y[i];
        if(fabs(errorz1[i]>0.000001)) delta[0] = (z1[i]-func1(&xx,par))/errorz1[i];
        else delta[0]=0;
        if(fabs(errorz2[i]>0.000001)) delta[1] = (z2[i]-func2(&yy,par))/errorz2[i];
        else delta[1]=0;
        chisq += (delta[0]*delta[0])+(delta[1]*delta[1]);
        chisq1 += delta[0]*delta[0];
        chisq2 += delta[1]*delta[1];
    }
    f=chisq;
    chisquare1=chisq1;
    chisquare2=chisq2;
    chisquare=chisq;
}


void Systematic()
{
    static Double_t vstart[4][NPT][5];
    static Double_t step[4][NPT][5];
    
    double PtEdge[7] = {0,2,3,4,6,8,14};
    
    TH1F *hist[NFILE][4][NPT][NFRAME][2];
    TH1F *hist1[NFILE][4][NPT][NFRAME][2];
    TH1F *defaulthist[4][NPT][2][2];
    
    TFile *inputfile[NFILE];
    TFile *defaultinputfile;
    
    int nbin;
    
    Double_t arglist[10];
    Int_t ierflg = 0;
    bool flag[4][NPT][2][2][4];
    int  empty[TRIG][NPT][NFRAME][EID][4];
    Int_t  fitflag[4][NPT][2][2][4];
    Double_t Chisquare[4][NPT][2][2];
    double fParamVal[NFILE][4][NPT][2][8];
    double fParamErr[NFILE][4][NPT][2][8];
    
    Double_t maximum[10][TRIG][NPT][2][2];
    Double_t maximumphi[10][TRIG][NPT][2][2];
    Double_t uncertainty[6][TRIG][NPT][2][EID];
    Double_t uncertaintyphi[6][TRIG][NPT][2][EID];
    
    Double_t matrix[4][4];
    Double_t covariant[TRIG][NPT][2][EID];
    
    for(int trig=0;trig<4;trig++){
        for(int npt=0;npt<NPT;npt++){
            for(int nframe=0;nframe<2;nframe++){
                for(int eID=0;eID<2;eID++){
                    for(int kk=0;kk<10;kk++){
                        maximum[kk][trig][npt][nframe][eID]=0;
                        maximumphi[kk][trig][npt][nframe][eID]=0;
                    }
                }
            }
        }
    }
    
    TGraphErrors *Uncertainties[6][TRIG][NFRAME][EID];
    
    Double_t LAMBDAERR[NPT],PT[NPT],LAMBDA[NPT],PTERR[NPT];
    Double_t LAMBDAPHI[NPT],LAMBDAPHIERR[NPT];
    Double_t LAMBDAINV[NPT],LAMBDAINVERR[NPT];
    
    TF1 *fit1[NPT][2][2]; //npt nframe default/revised
    TF1 *fit2[NPT][2][2];
    for(int nfile=0;nfile<1;nfile++){
        if(nfile==38) continue;
        for(int eID =1;eID<3;eID++){
            inputfile[nfile] = new TFile(Form("~/jpsi/test20160210_%d/corrected%deID.root",nfile,eID),"read");
            for(int trig=0;trig<TRIG;trig++){
                for(int npt=0;npt<NPT;npt++){
                    for(int nframe=0;nframe<2;nframe++){
                        for(int var=0;var<4;var++){
                            flag[trig][npt][nframe][eID-1][var]=false;
                            empty[trig][npt][nframe][eID-1][var]=0;
                            fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]=0;
                            fParamErr[nfile][trig][npt][nframe][1+4*(eID-1)]=0;
                        }
                    }
                }
            }
            
            for(int trig=0;trig<TRIG;trig++){
                for(int npt=0;npt<NPT;npt++){
                    for(int nframe=0;nframe<NFRAME;nframe++){
                        
                        vstart[trig][npt][0]=0.1;
                        vstart[trig][npt][1]=0.1;
                        vstart[trig][npt][2]=0.1;
                        vstart[trig][npt][3]=-0.1;
                        vstart[trig][npt][4]=0.001;
                        
                        step[trig][npt][0]=0.001;
                        step[trig][npt][1]=0.001;
                        step[trig][npt][2]=0.001;
                        step[trig][npt][3]=0.001;
                        step[trig][npt][4]=0.001;
                        
                        if(inputfile[nfile]->Get(Form("HistJpsiPhase_%d_%d_%d",trig,npt,2*nframe))!=0x0) {
                            flag[trig][npt][nframe][eID-1][0]=true;
                            flag[trig][npt][nframe][eID-1][1]=true;
                            hist[nfile][trig][npt][nframe][0] = (TH1F*)inputfile[nfile]->Get(Form("HistJpsiPhase_%d_%d_%d",trig,npt,2*nframe));
                            hist[nfile][trig][npt][nframe][0]->SetName(Form("HistPhase_%dnfile_%dtrig_%dnpt_%dnframe_%deID_theta",nfile,trig,npt,nframe,eID));
                            if(eID==1) {
                                hist1[nfile][trig][npt][nframe][0]=(TH1F*)hist[nfile][trig][npt][nframe][0]->Clone();// save the 1eID default corrected distribution
                                hist1[nfile][trig][npt][nframe][0]->SetName(Form("1eID_%dnfile_%dtrig_%dnpt_%dnframe_theta",nfile,trig,npt,nframe,eID));
                            }
                            for(int nbin=1;nbin<11;nbin++){
                                z1[nbin-1] = hist[nfile][trig][npt][nframe][0]->GetBinContent(nbin);
                                errorz1[nbin-1] = hist[nfile][trig][npt][nframe][0]->GetBinError(nbin);
                                x[nbin-1]=hist[nfile][trig][npt][nframe][0]->GetBinCenter(nbin);
                                if(fabs(z1[nbin-1]<0.0001)) empty[trig][npt][nframe][eID-1][0]++;
                            }
                        }
                        if(inputfile[nfile]->Get(Form("HistJpsiPhase_%d_%d_%d",trig,npt,2*nframe+1))!=0x0) {
                            flag[trig][npt][nframe][eID-1][2]=true;
                            flag[trig][npt][nframe][eID-1][3]=true;
                            hist[nfile][trig][npt][nframe][1] = (TH1F*)inputfile[nfile]->Get(Form("HistJpsiPhase_%d_%d_%d",trig,npt,2*nframe+1));
                            hist[nfile][trig][npt][nframe][1]->SetName(Form("HistPhase_%d_%d_%d_%d_phi",trig,npt,nframe,eID));
                            if(eID==1) {
                                hist1[nfile][trig][npt][nframe][1]=(TH1F*)hist[nfile][trig][npt][nframe][1]->Clone();// save the 1eID default corrected distribution
                                hist1[nfile][trig][npt][nframe][1]->SetName(Form("1eID_%dnfile_%dtrig_%dnpt_%dnframe_phi",nfile,trig,npt,nframe,eID));
                            }
                            for(int nbin=1;nbin<11;nbin++){
                                z2[nbin-1] = hist[nfile][trig][npt][nframe][1]->GetBinContent(nbin);
                                errorz2[nbin-1] = hist[nfile][trig][npt][nframe][1]->GetBinError(nbin);
                                y[nbin-1]=hist[nfile][trig][npt][nframe][1]->GetBinCenter(nbin);
                                if(fabs(z2[nbin-1])<0.0001)empty[trig][npt][nframe][eID-1][1]++;
                            }
                        }
                        int icon;
                        TMinuit *gMinuit = new TMinuit(5);
                        gMinuit->SetFCN(fcn);
                        
                        arglist[0] = 1;
                        gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
                        
                        gMinuit->mnparm(0,"norm1",vstart[trig][npt][0],step[trig][npt][0],0,0,ierflg);
                        gMinuit->mnparm(1,"lambda1",vstart[trig][npt][1],step[trig][npt][1],0,0,ierflg);
                        fitflag[trig][npt][nframe][eID-1][1]=ierflg;
                        gMinuit->mnparm(2,"lambda2",vstart[trig][npt][2],step[trig][npt][2],0,0,ierflg);
                        fitflag[trig][npt][nframe][eID-1][2]=ierflg;
                        gMinuit->mnparm(3,"norm2",vstart[trig][npt][3],step[trig][npt][3],0,0,ierflg);
                        
                        arglist[0] = 500.;
                        arglist[1] = 1.;
                        gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
                        Double_t amin,edm,errdef;
                        Int_t nvpar,nparx,icstat;
                        gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
                        Chisquare[trig][npt][nframe][0] = chisquare1;
                        Chisquare[trig][npt][nframe][1] = chisquare2;
                        for(int var=0;var<4;var++){
                            if(flag[trig][npt][nframe][eID-1][var]==true){
                                gMinuit->GetParameter(var,fParamVal[nfile][trig][npt][nframe][var+4*(eID-1)],fParamErr[nfile][trig][npt][nframe][var+4*(eID-1)]);
                                //			cout<<"fParamVal============"<<"nfile"<<nfile<<"   trig"<<trig<<"   npt"<<npt<<"  nframe"<<nframe<<"   eID"<<eID<<"    "<<fParamVal[nfile][trig][npt][nframe][var+4*(eID-1)]<<endl;
                            }
                            else {
                                fParamVal[nfile][trig][npt][nframe][var+4*(eID-1)]=0;
                                fParamErr[nfile][trig][npt][nframe][var+4*(eID-1)]=0;
                            }
                        }
                        gMinuit->mnemat(&matrix[0][0],4);
                        covariant[trig][npt][nframe][eID-1] = matrix[1][2];
                    }
                }
            }
            
            gStyle->SetOptStat(0);
            gStyle->SetFillStyle(0);
            
            TCanvas *canvas[5][3];
            TLegend *legend;
            TLatex text;
            TLatex CHI2;
            TLatex deltalambda;
            gStyle->SetLegendBorderSize(0);
            
            canvas[0][0] = new TCanvas(Form("lambda%deID",eID),Form("lambda%deID",eID),1200,800);
            canvas[0][0]->Divide(3,2);
            TGraphErrors *plot[4][6][NFRAME][EID];
            
            double lambda[NPT],lambda_err[NPT],lambdaphi[NPT],lambdaphierr[NPT],lambdainv[NPT],lambdainverr[NPT];
            double minimum_err=1.0;
            double pt[NPT],pterr[NPT];
            for(int nplot=0;nplot<6;nplot++){
                if(nplot==0) var=1,nframe=0;
                if(nplot==1) var=2,nframe=0;
                //		if(nplot==2) continue;
                if(nplot==3) var=1,nframe=1;
                if(nplot==4) var=2,nframe=1;
                //		if(nplot==5) continue;
                canvas[0][0]->cd(nplot+1);
                legend = new TLegend(0.75,0.4,0.89,0.89);
                
                for(int trig=0;trig<4;trig++){
                    for(int npt=0;npt<NPT;npt++){
                        pt[npt]=0.5*(PtEdge[npt]+PtEdge[npt+1]);
                        pterr[npt]=0.5*(PtEdge[npt+1]-PtEdge[npt]);
                        //			if(flag[trig][npt][nframe][eID-1][var]==kTRUE && fitflag[trig][npt][nframe][eID-1][var]==0 && empty[trig][npt][nframe][eID-1][var]<=4 && !(trig==3 && npt==2) &&!(trig==1 && npt==0 &&(nframe==0 || nframe==3) &&eID==2 )) {
                        if(flag[trig][npt][nframe][eID-1][var]==kTRUE && fitflag[trig][npt][nframe][eID-1][var]==0 && empty[trig][npt][nframe][eID-1][var]<=4 && !(trig==3 && npt==2) ) {
                            lambda[npt] = fParamVal[0][trig][npt][nframe][var+4*(eID-1)];
                            lambda_err[npt] = fParamErr[0][trig][npt][nframe][var];
                            minimum_err=fParamErr[0][trig][npt][nframe][var]<minimum_err?fParamErr[0][trig][npt][nframe][var]:minimum_err;
                            //       if((nplot==2 || nplot==3) && npt>=3){
                            //         lambda[npt] = -100.;
                            //           lambda_err[npt] =0.;
                            //     }
                            //   else{
                            //     lambda[npt] = -100.;
                            //   lambda_err[npt] =0.;
                            // }
                        }
                        else{
                            lambda[npt] = -100.;
                            lambda_err[npt] =0.;
                        }
                        //	if(nplot==1 || nplot==4)lambdaphi[npt]= fParamVal[0][trig][npt][nframe][2+4*(eID-1)];
                        if(nplot==2 || nplot==5)lambda[npt] = (fParamVal[0][trig][npt][nframe][1+4*(eID-1)]+3*fParamVal[0][trig][npt][nframe][2+4*(eID-1)])/(1.-fParamVal[0][trig][npt][nframe][2+4*(eID-1)]);
                        //cout<<"fParamVal1========="<<"   trig"<<trig<<"   npt"<<npt<<"   nframe"<<nframe<<"   eID"<<eID<<"    "<<fParamVal[0][trig][npt][nframe][1+4*eID]<<"  lambda[npt]===="<<lambda[npt]<<"   theta"<<fParamVal[0][trig][npt][nframe][1+4*(eID-1)]<<"  phi"<<fParamVal[0][trig][npt][nframe][2+4*(eID-1)]<<endl;
                        if(nplot==2 || nplot==5)lambda_err[npt] = TMath::Sqrt(fParamErr[nfile][trig][npt][nframe][1+4*(eID-1)]*fParamErr[nfile][trig][npt][nframe][1+4*(eID-1)]/((1-fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)])*(1-fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]))+fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]*fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]*TMath::Power((3+fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)])/((1-fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)])*(1-fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)])),2)+2*(3+fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)])/TMath::Power((1-fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]),3)*covariant[trig][npt][nframe][eID-1]);
                    }
                    plot[trig][nplot][nframe][eID-1]=new TGraphErrors(NPT,pt,lambda,pterr,lambda_err);
                    plot[trig][nplot][nframe][eID-1]->SetMarkerStyle(20+trig+4*(eID-1));
                    plot[trig][nplot][nframe][eID-1]->SetMarkerColor(1+trig+4*(eID-1));
                    plot[trig][nplot][nframe][eID-1]->SetLineColor(1+trig+4*(eID-1));
                    
                    plot[trig][nplot][nframe][eID-1]->SetName(Form("plot_%s_%s_%d_%deID",TrigName[trig].Data(),VarName[var-1].Data(),nframe,eID));
                    plot[trig][nplot][nframe][eID-1]->SetTitle(Form("Parameter #lambda_{#%s} in %s frame",VarName[var-1].Data(),FrameName[nframe].Data()));
                    plot[trig][nplot][nframe][eID-1]->GetYaxis()->SetRangeUser(-2,2);
                    plot[trig][nplot][nframe][eID-1]->GetXaxis()->SetRangeUser(0,14);
                    plot[trig][nplot][nframe][eID-1]->GetXaxis()->SetTitle("Jpsi Pt");
                    plot[trig][nplot][nframe][eID-1]->GetXaxis()->SetLabelSize();
                    plot[trig][nplot][nframe][eID-1]->GetYaxis()->SetTitle(Form("#lambda_{#%s} in %s frame",VarName[var-1].Data(),FrameName[nframe].Data()));
                    if(trig==2)continue;   // do not plot HT1 right now
                    if (nplot==2 || nplot==5) {
                        continue;
                    }
                    if(trig==1)plot[trig][nplot][nframe][eID-1]->Draw("ap");
                    else plot[trig][nplot][nframe][eID-1]->Draw("psame");
                    if(trig!=0)	{
                        legend->AddEntry(plot[trig][nplot][nframe][eID-1],Form("2D HT%d",trig-1),"p");
                    }
                    else {
                        //legend->AddEntry(plot[trig][nplot][nframe][eID-1],"2D MB","p");
                    }
                    legend->Draw("same");
                }
            }
            
            canvas[0][0]->SaveAs(Form("lambda%deID.pdf",eID));
            gStyle->SetLabelSize(0.3);
            for(int trig=0;trig<4;trig++){
                canvas[trig][0] = new TCanvas(Form("%s_%d_canvas",TrigName[trig].Data(),eID),Form("%s_%d_canvas",TrigName[trig].Data(),eID),1500,1200);
                canvas[trig][0]->Divide(NPT,4);
                for(int npt=0;npt<NPT;npt++){
                    for(int nframe=0;nframe<NFRAME;nframe++){
                        fit1[npt][nframe][0] = new TF1(Form("fit1_%d_%d",npt,0),func1,-1,1,2);
                        fit2[npt][nframe][0] = new TF1(Form("fit2_%d_%d",npt,0),func2,-TMath::Pi(),TMath::Pi(),4);
                        fit1[npt][nframe][1] = new TF1(Form("fit1_%d_%d",npt,1),func1,-1,1,2);
                        fit2[npt][nframe][1] = new TF1(Form("fit2_%d_%d",npt,1),func2,-TMath::Pi(),TMath::Pi(),4);
                        
                        for(int var=0;var<2;var++){
                            canvas[trig][0]->cd(npt+NPT*var+2*NPT*nframe+1);
                            if(flag[trig][npt][nframe][eID-1][var+1]==true){
                                if(eID==1){
                                    if(hist1[0][trig][npt][nframe][var]!=0x0 && hist1[nfile][trig][npt][nframe][var]!=0x0){hist1[0][trig][npt][nframe][var]->Draw();
                                        hist1[nfile][trig][npt][nframe][var]->Draw("same");
                                        hist1[0][trig][npt][nframe][var]->SetMarkerStyle(22);
                                        hist1[nfile][trig][npt][nframe][var]->SetMarkerStyle(26);
                                    }
                                }
                                else{
                                    if(hist[nfile][trig][npt][nframe][var]!=0x0){
                                        hist[0][trig][npt][nframe][var]->Draw();
                                        hist[nfile][trig][npt][nframe][var]->Draw("same");
                                        hist[0][trig][npt][nframe][var]->SetMarkerStyle(22);
                                        hist[nfile][trig][npt][nframe][var]->SetMarkerStyle(26);
                                    }
                                }
                                hist[nfile][trig][npt][nframe][var]->SetTitle(Form("%s %1.f GeV/c < Pt <%1.f GeV/c %deID",TrigName[trig].Data(),PtEdge[npt],PtEdge[npt+1],eID));
                                hist[nfile][trig][npt][nframe][var]->GetYaxis()->SetTitle(Form("%s %s",FrameName[nframe].Data(),VarName[var].Data()));
                                fit1[npt][nframe][0]->SetParameters(fParamVal[0][trig][npt][nframe][0+4*(eID-1)],fParamVal[0][trig][npt][nframe][1+4*(eID-1)]);
                                fit2[npt][nframe][0]->SetParameters(fParamVal[0][trig][npt][nframe][0+4*(eID-1)],fParamVal[0][trig][npt][nframe][1+4*(eID-1)],fParamVal[0][trig][npt][nframe][2+4*(eID-1)],fParamVal[0][trig][npt][nframe][3+4*(eID-1)]);
                                fit1[npt][nframe][1]->SetParameters(fParamVal[nfile][trig][npt][nframe][0+4*(eID-1)],fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]);
                                fit2[npt][nframe][1]->SetParameters(fParamVal[nfile][trig][npt][nframe][0+4*(eID-1)],fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)],fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)],fParamVal[nfile][trig][npt][nframe][3+4*(eID-1)]);
                                
                                cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
                                text.SetTextFont(43);
                                text.SetTextSize(14);
                                text.SetNDC();
                                CHI2.SetTextFont(43);
                                CHI2.SetTextSize(14);
                                CHI2.SetNDC();
                                deltalambda.SetTextFont(43);
                                deltalambda.SetTextSize(14);
                                deltalambda.SetNDC();
                                
                                if(var==0){
                                    fit1[npt][nframe][0]->Draw("same");
                                    fit1[npt][nframe][1]->Draw("same");
                                    deltalambda.DrawLatex(0.2,0.8,Form("#Delta#lambda_{#theta} = %1.2f-%1.2f = %1.2f",fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)],fParamVal[0][trig][npt][nframe][1+4*(eID-1)],fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)]));
                                }
                                else {
                                    fit2[npt][nframe][0]->Draw("same");
                                    fit2[npt][nframe][1]->Draw("same");
                                    deltalambda.DrawLatex(0.2,0.8,Form("#Delta#lambda_{#phi} =%1.2f-%1.2f= %1.2f",fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)],fParamVal[0][trig][npt][nframe][2+4*(eID-1)],fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)]));
                                }
                            }
                        }
                        if(nfile>=1 && nfile<=11){
                            cout<<"00000000000000000==========="<<maximum[0][trig][npt][nframe][eID-1]<<endl;
                            maximum[0][trig][npt][nframe][eID-1]+=(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)])/11;
                            maximumphi[0][trig][npt][nframe][eID-1] +=(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])/11;
                            cout<<"flag============"<<flag[trig][npt][nframe][eID-1][1]<<endl;
                            cout<<"irflag==========="<<ierflg<<endl;
                            cout<<"1111111111111111111===="<<nfile<<" trig"<<trig<<"  npt"<<npt<<"   eID"<<eID<<"        "<<(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])<<"     "<<maximumphi[0][trig][npt][nframe][eID-1]<<endl;
                        }
                        if(nfile>=12 && nfile<=13){
                            maximum[1][trig][npt][nframe][eID-1]+=(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)])/2;
                            maximumphi[1][trig][npt][nframe][eID-1] +=(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])/2;
                            cout<<"222222222222222222==========="<<nfile<<"   trig"<<trig<<"   npt"<<npt<<"  eID"<<eID<<"    "<<(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])<<"           "<<maximumphi[1][trig][npt][nframe][eID-1]<<endl;
                        }
                        if(nfile>=14 && nfile<=21){//weight
                            maximum[2][trig][npt][nframe][eID-1]+=(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)])/8;
                            maximumphi[2][trig][npt][nframe][eID-1]+=(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])/8;
                            cout<<"333333333333333============== "<<nfile<<"  trig"<<trig<<"   npt"<<npt<<"   eID"<<eID<<"        "<<(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])<<"       "<<maximumphi[2][trig][npt][nframe][eID-1]<<endl;
                        }
                        if(nfile>=22 && nfile<=25){ //nhistfit
                            cout<<"444444444444444444444444444444"<<maximum[3][trig][npt][nframe][eID-1]<<endl;
                            if((fabs(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)]))>=maximum[3][trig][npt][nframe][eID-1]) maximum[3][trig][npt][nframe][eID-1] = fabs(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)]);
                            if((fabs(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)]))>=maximumphi[3][trig][npt][nframe][eID-1]) maximumphi[3][trig][npt][nframe][eID-1] = fabs(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)]);
                            cout<<"44444444444444444444=========="<<nfile<<"    trig"<<trig<<"    npt"<<npt<<"   eID"<<eID<<"      "<<(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])<<"         "<<maximumphi[3][trig][npt][nframe][eID-1]<<endl;
                        }
                        if(nfile>=26 && nfile<=29){//nsigma
                            cout<<"55555555555555555555555555555555"<<maximum[4][trig][npt][nframe][eID-1]<<endl;
                            if((fabs(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)]))>=maximum[4][trig][npt][nframe][eID-1]) maximum[4][trig][npt][nframe][eID-1] = fabs(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)]);
                            
                            if((fabs(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)]))>=maximumphi[4][trig][npt][nframe][eID-1]) maximumphi[4][trig][npt][nframe][eID-1] = fabs(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)]);
                            cout<<"55555555555555555555555================="<<nfile<<"   trig"<<trig<<"    npt"<<npt<<"   eID"<<eID<<"       "<<(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])<<"          "<<maximumphi[4][trig][npt][nframe][eID-1]<<endl;
                        }
                        if(nfile>=30 && nfile<=31){  // trigger
                            maximum[5][trig][npt][nframe][eID-1]+=(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)])/2;
                            maximumphi[5][trig][npt][nframe][eID-1]+=(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])/2;
                            cout<<"6666666666666666666==============="<<nfile<<"   trig"<<trig<<"  npt"<<npt<<"   eID"<<eID<<"        "<<(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])<<"          "<<maximumphi[5][trig][npt][nframe][eID-1]<<endl;
                        }
                        if(nfile>=32 && nfile<=35){
                            if(fabs(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)])>=maximum[6][trig][npt][nframe][eID-1]) maximum[6][trig][npt][nframe][eID-1] = fabs(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)]);
                            
                            if(fabs(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])>=maximumphi[6][trig][npt][nframe][eID-1]) maximumphi[6][trig][npt][nframe][eID-1] = fabs(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)]);
                            cout<<"77777777777777==================="<<nfile<<"  trig"<<trig<<"   mpt"<<npt<<"   eID"<<eID<<"    "<<(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])<<"           "<<maximumphi[6][trig][npt][nframe][eID-1]<<endl;
                        }
                        if(nfile>=39 && nfile<=40){
                            maximum[7][trig][npt][nframe][eID-1]+=(fParamVal[nfile][trig][npt][nframe][1+4*(eID-1)]-fParamVal[0][trig][npt][nframe][1+4*(eID-1)])/2;
                            maximumphi[7][trig][npt][nframe][eID-1]+=(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])/2;
                            cout<<"88888888888888888888==============="<<nfile<<"          "<<"         "<<(fParamVal[nfile][trig][npt][nframe][2+4*(eID-1)]-fParamVal[0][trig][npt][nframe][2+4*(eID-1)])<<"           "<<maximumphi[7][trig][npt][nframe][eID-1]<<endl;
                        }
                    }
                    canvas[trig][0]->SaveAs(Form("Uncertainties/%d_%s_%deID.pdf",nfile,TrigName[trig].Data(),eID));
                }
            }
        }
    }
    canvas[0][3] = new TCanvas("uncertainties","uncertainties",1200,800);
    canvas[0][3]->Divide(3,2);
    
    for(int eID=0;eID<2;eID++){
        for(int nplot=0;nplot<NPLOT;nplot++){
            canvas[0][3]->cd(nplot+1);
            
            legend = new TLegend(0.7,0.7,0.89,0.89);
            if(nplot==0) nframe=0;
            if(nplot==1) nframe=0;
            if(nplot==2) nframe=0;
            if(nplot==3) nframe=1;
            if(nplot==4) nframe=1;
            if(nplot==5) nframe=1;
            
            for(int trig=1;trig<4;trig++){
                for(int npt=0;npt<NPT;npt++){
                    for(int i=0;i<8;i++){
                        //		cout<<"i======="<<i<<"trig"<<trig<<"  npt"<<npt<<"   eID"<<eID<<"maximumphi================"<<maximumphi[i][trig][npt][nframe][eID]*maximumphi[i][trig][npt][nframe][eID]<<endl;
                        if(maximum[i][trig][nframe][eID]>1) cout<<"i======="<<i<<"trig"<<trig<<"  npt"<<npt<<"   eID"<<eID<<"maximum========="<<maximum[i][trig][npt][nframe][eID]<<endl;
                        if(maximumphi[i][trig][nframe][eID]>1) cout<<"i======="<<i<<"trig"<<trig<<"  npt"<<npt<<"   eID"<<eID<<"maximumphi========="<<maximumphi[i][trig][npt][nframe][eID]<<endl;
                        uncertainty[nplot][trig][npt][nframe][eID] += maximum[i][trig][npt][nframe][eID]*maximum[i][trig][npt][nframe][eID];
                        uncertaintyphi[nplot][trig][npt][nframe][eID] += maximumphi[i][trig][npt][nframe][eID]*maximumphi[i][trig][npt][nframe][eID];
                    }
                    if(uncertainty[nplot][trig][npt][nframe][eID]>1) cout<<" sum uncertainty========"<<uncertainty[nplot][trig][npt][nframe][eID]<<endl;
                    if(uncertaintyphi[nplot][trig][npt][nframe][eID]>1)		cout<<"sum uncertaintyphi==============="<<uncertaintyphi[nplot][trig][npt][nframe][eID]<<endl;
                    uncertainty[nplot][trig][npt][nframe][eID]=TMath::Sqrt(uncertainty[nplot][trig][npt][nframe][eID]);
                    uncertaintyphi[nplot][trig][npt][nframe][eID]=TMath::Sqrt(uncertaintyphi[nplot][trig][npt][nframe][eID]);
                    cout<<"=======================uncertainty"<<uncertainty[nplot][trig][npt][nframe][eID]<<endl;
                    cout<<" ======================uncertaintyphi"<<uncertaintyphi[nplot][trig][npt][nframe][eID]<<endl;
                    
                    PT[npt] = 0.5*(PtEdge[npt]+PtEdge[npt+1]);
                    PTERR[npt] = 0.15;
                    if(flag[trig][npt][nframe][eID][1]==kTRUE && fitflag[trig][npt][nframe][eID][1]==0) {
                        LAMBDA[npt] = fParamVal[0][trig][npt][nframe][1+4*(eID)];
                        LAMBDAERR[npt] = uncertainty[nplot][trig][npt][nframe][eID];
                    }
                    else {
                        LAMBDA[npt] = -100;
                        LAMBDAERR[npt] = 0.;
                    }
                    
                    if(flag[trig][npt][nframe][eID][2]==kTRUE && fitflag[trig][npt][nframe][eID][2]==0){					   
                        LAMBDAPHI[npt] = fParamVal[0][trig][npt][nframe][2+4*(eID)];
                        LAMBDAPHIERR[npt] = uncertaintyphi[nplot][trig][npt][nframe][eID];	
                    }
                    else {
                        LAMBDAPHI[npt] = -100;
                        LAMBDAPHIERR[npt] = 0.;	
                    }
                    
                    if(flag[trig][npt][nframe][eID][1]==kTRUE && flag[trig][npt][nframe][eID][2]==kTRUE && (LAMBDA[npt]>-99 && LAMBDAPHI[npt]>-99 )){
                        LAMBDAINV[npt] = (fParamVal[0][trig][npt][nframe][1+4*(eID)] + 3.*fParamVal[0][trig][npt][nframe][2+4*(eID)])/(1.-fParamVal[0][trig][npt][nframe][2+4*(eID)]);
                        //						LAMBDAINVERR[npt] = TMath::Sqrt(uncertainty[nplot][trig][npt][nframe][eID]*uncertainty[nplot][trig][npt][nframe][eID]/((1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt]))+uncertaintyphi[nplot][trig][npt][nframe][eID]*uncertaintyphi[nplot][trig][npt][nframe][eID]*((3+LAMBDA[npt])*(3+LAMBDA[npt])/((1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt])))+2*(3+LAMBDA[npt])/((1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt]))*covariant[trig][npt][nframe][eID]*covariant[trig][npt][nframe][eID]);
                        LAMBDAINVERR[npt] = TMath::Sqrt(LAMBDAERR[npt]*LAMBDAERR[npt]/((1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt]))+LAMBDAPHIERR[npt]*LAMBDAPHIERR[npt]*((3+LAMBDA[npt])*(3+LAMBDA[npt])/((1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt])))+2*(3+LAMBDA[npt])/((1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt])*(1-LAMBDAPHI[npt]))*covariant[trig][npt][nframe][eID]*covariant[trig][npt][nframe][eID]);
                    }
                    else{
                        LAMBDAINV[npt] = -100;
                        LAMBDAINVERR[npt] = 0.;
                    }
                }
                if(nplot==0 || nplot==3){
                    Uncertainties[nplot][trig][nframe][eID] = new TGraphErrors(NPT,PT,LAMBDA,PTERR,LAMBDAERR);
                }
                if(nplot==1 || nplot==4) {
                    Uncertainties[nplot][trig][nframe][eID] = new TGraphErrors(NPT,PT,LAMBDAPHI,PTERR,LAMBDAPHIERR);
                }	
                if(nplot==2 || nplot==5){
                    Uncertainties[nplot][trig][nframe][eID] = new TGraphErrors(NPT,PT,LAMBDAINV,PTERR,LAMBDAINVERR);
                }
                Uncertainties[nplot][trig][nframe][eID]->SetName(Form("uncertainties_%dtrig_%dnframe_%deID_%dnplot",trig,nframe,eID,nplot));
                Uncertainties[nplot][trig][nframe][eID]->GetXaxis()->SetRangeUser(0,14);
                if(nplot==2 || nplot==5) Uncertainties[nplot][trig][nframe][eID]->GetYaxis()->SetRangeUser(-5,5);
                else Uncertainties[nplot][trig][nframe][eID]->GetYaxis()->SetRangeUser(-2,2);
                
                Uncertainties[nplot][trig][nframe][eID]->GetXaxis()->SetLabelSize();
                Uncertainties[nplot][trig][nframe][eID]->SetFillColor(1+trig+eID*4);
                Uncertainties[nplot][trig][nframe][eID]->SetMarkerSize();
                
                if(trig==1){
                    Uncertainties[nplot][trig][nframe][eID]->Draw("ap");
                }
                else Uncertainties[nplot][trig][nframe][eID]->Draw("psame");
                
                legend->AddEntry(Uncertainties[nplot][trig][nframe][eID],TrigName[trig].Data(),"p");				
                
                
                plot[trig][nplot][nframe][eID]->SetMarkerColor(1+trig+eID*4);
                plot[trig][nplot][nframe][eID]->SetLineColor(1+trig+eID*4);
                /*				cout<<"   nplot"<<nplot<<"  trig"<<trig<<"  nframe"<<nframe<<"   eID"<<eID<<endl;
                 cout<<"////////////////======================//////////////////////////////////////////////////"<<endl;
                 plot[trig][nplot][nframe][eID]->Print("all");
                 cout<<"////////////////======================//////////////////////////////////////////////////"<<endl;
                 Uncertainties[nplot][trig][nframe][eID]->Print("all");
                 */
                plot[trig][nplot][nframe][eID]->SetMarkerStyle(20+trig+eID*4);
                Uncertainties[nplot][trig][nframe][eID]->SetMarkerStyle(20+trig+eID*4);
                Uncertainties[nplot][trig][nframe][eID]->SetMarkerColor(1+trig);
                Uncertainties[nplot][trig][nframe][eID]->SetLineColor(1+trig+eID*4);
            }
            
            legend->Draw("same");
        }
        canvas[0][3]->SaveAs(Form("test%deID.pdf",eID+1));
    }
}
