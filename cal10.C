#define NTRIG 4
#define NPHASE 4
#define NPT 6
#define NDRAW 4
#define NPHI 10
#define NFRAME 2
#define NVAR 3
#define EID 2

Double_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};
TString TrigName[NTRIG]={"MB","HT0","HT1","HT2"};
TString PhaseName[NPHASE]={"HXTheta","HXPhi","CSTheta","CSPhi"};
TString FrameName[NFRAME]={"HX","CS"};
//TString VarName[NVAR]={"#lambda_{#theta}","#lambda_{#phi}","#lambda_{inv}"};
TString VarName[NVAR]={"lambda_theta","lambda_phi","lambda_inv"};
void addtext(int n, int k, int l int eID, TH1F *hist);

void cal10(){
	TFile *datafile[NTRIG];
	TFile *mcfile = new TFile("rootfile/OutFile_cent_0_9.root");
	Double_t xx[NPT],xxerr[NPT];
	Double_t LambdaTheta[NPT];
	Double_t LambdaThetaErr[NPT];
	Double_t LambdaPhi[NPT];
	Double_t LambdaPhiErr[NPT];
	Double_t LambdaInv[NPT];
	Double_t LambdaInvErr[NPT];

	Double_t datapt[3] = {2.48,3.52,4.74};
	Double_t datapterr[3]={0,0,0};
	Double_t datalambda[3] = {0.145,-0.476,-0.617};
	Double_t datasyst[3] = {0.297,0.160,0.260};

	TGraphErrors *run09=new TGraphErrors(3,datapt,datalambda,datapterr,datasyst);
	run09->SetMarkerStyle(24);
	run09->SetMarkerColor(2);
	run09->SetLineColor(2);

	TGraphAsymmErrors *Eff[NTRIG][NPHASE][NPT];
	TGraphErrors *Lambda_graph[NTRIG][NPHASE][3]; // MB HT0 HT1 HT2; theta_hx phi_hx theta_cs phi_cs; theta phi invariant;

	TF1 *func1;
	TF1 *func2;

	TH1F *HistEff[NTRIG][NPHASE][NPT];
	TH2F *HistPhasePtUL[NTRIG][NPHASE];
	TH2F *HistPhasePtLI[NTRIG][NPHASE];
	TH1F *HistPhase0[NTRIG][NPHASE][NPT];
	TH1F *HistPhase1[NTRIG][NPHASE][NPT];
	TH1F *HistPhase2[NTRIG][NPHASE][NPT];
	TH1F *HistPhase[NTRIG][NPHASE][NPT];
	TH2F *DefaultTH2F[NTRIG][NPHASE][2];
	TH1F *DefaultTH1F[NTRIG][NPHASE][NPT][3];

	TFile *DefaultFile[NTRIG][EID];

	TCanvas *canvas[NTRIG][NPHASE][NDRAW];
	TLegend *rawdatalegend;
    TLine *line;

	TFile *raw;
	TFile *efficiency;
	TFile *corrected;
	TFile *Lambda;
	TCanvas *ncanvas;

	TString datafilename[4]={"rootfile/mb_1TrkPid_0.ana.root","rootfile/ht0_trg3_1TrkPid_0.ana.root","rootfile/ht1_trg4_1TrkPid_0.ana.root","rootfile/ht2_trg5_1TrkPid_0.ana.root"};

	TString embhistname2[]={"hMcJpsiCosThetaPt","hMBJpsiCosThetaPt","hHt0JpsiCosThetaPt","hHt1JpsiCosThetaPt","hHt2JpsiCosThetaPt","hMcJpsiPhiPt","hMBJpsiPhiPt","hHtJpsiPhiPt","hHt1JpsiPhiPt","hHt2JpsiPhiPt","hMcJpsiCosThetaPtCS","hMBJpsiCosThetaPtCS","hHt0JpsiCosThetaPtCS","hHt1JpsiCosThetaPtCS","hHt2JpsiCosThetaPtCS","hMcJpsiPhiPtCS","hMBJpsiPhiPtCS","hHt0JpsiPhiPtCS","hHt1JpsiPhiPtCS","hHt2JpsiPhiPtCS"};
	TString embhistname1[NPHASE*(NTRIG+1)]={"hMcJpsiCosThetaPt","hMBJpsiCosThetaPt1","hHt0JpsiCosThetaPt1","hHt1JpsiCosThetaPt1","hHt2JpsiCosThetaPt1","hMcJpsiPhiPt","hMBJpsiPhiPt1","hHt0JpsiPhiPt1","hHt1JpsiPhiPt1","hHt2JpsiPhiPt1","hMcJpsiCosThetaPtCS","hMBJpsiCosThetaPtCS1","hHt0JpsiCosThetaPtCS1","hHt1JpsiCosThetaPtCS1","hHt2JpsiCosThetaPtCS1","hMcJpsiPhiPtCS","hMBJpsiPhiPtCS1","hHt0JpsiPhiPtCS1","hHt1JpsiPhiPtCS1","hHt2JpsiPhiPtCS1"};


	TString datahistname[NPHASE*2] = {"hJpsiCosThetaPt","hJpsiCosThetaPtBG","hJpsiPhiPtHX","hJpsiPhiPtHXBG","hJpsiCosThetaPtCS","hJpsiCosThetaPtCSBG","hJpsiPhiPtCS","hJpsiPhiPtCSBG"};


	TString datafilename1[NTRIG]={"rootfile/mb_1TrkPid_0.ana.root","rootfile/ht0_trg3_1TrkPid_0.ana.root","rootfile/ht1_trg4_1TrkPid_0.ana.root","rootfile/ht2_trg5_1TrkPid_0.ana.root"};

	TString datafilename2[NTRIG]={"rootfile/mb_2TrkPid_0.ana.root","rootfile/ht0_trg3_2TrkPid_0.ana.root","rootfile/ht1_trg4_2TrkPid_0.ana.root","rootfile/ht2_trg5_2TrkPid_0.ana.root"};

	TLatex STARpreliminary;
	TLatex ptbintext;
    
    ptbintext.SetTextSize(0.05);

	gStyle->SetFillStyle(0);	
	gStyle->SetLegendBorderSize(0);


	for(int eIDmode=1;eIDmode<3;eIDmode++){
		for(int trig=0;trig<4;trig++){
			if(trig==0)DefaultFile[trig][eIDmode-1]=new TFile(Form("../test20160210_0/rootfile/mb_%dTrkPid_0.ana.root",eIDmode));
			else DefaultFile[trig][eIDmode-1] = new TFile(Form("../test20160210_0/rootfile/ht%d_trg%d_%dTrkPid_0.ana.root",trig-1,trig+2,eIDmode));
		}

		//	TString datafilename[4]={"mb_2TrkPid_0.ana.root","ht0_trg3_2TrkPid_0.ana.root","ht1_trg4_2TrkPid_0.ana.root","ht2_trg5_2TrkPid_0.ana.root"};
		if(eIDmode==1){
			raw = new TFile("raw1eID.root","recreate");
			efficiency = new TFile("efficiency1eID.root","recreate");
			corrected = new TFile("corrected1eID.root","recreate");
		}
		if(eIDmode==2){
			raw = new TFile("raw2eID.root","recreate");
			efficiency = new TFile("efficiency2eID.root","recreate");
			corrected = new TFile("corrected2eID.root","recreate");
		}
		for (int trig=0; trig<NTRIG; trig++) {
			for (int nphase =0; nphase<NPHASE; nphase++) {
				int nframe=0;
				if(nphase<2) nframe=0;
				else nframe=1;
				for (int ndraw=0; ndraw<NDRAW; ndraw++) {
					canvas[trig][nphase][ndraw] = new TCanvas(Form("%s_%s_Draw%d",TrigName[trig].Data(),PhaseName[nphase].Data(),ndraw),"",1200,800);
					canvas[trig][nphase][ndraw]->Divide(3,2);
				}
				for (int npt=0;npt<NPT;npt++){
					int flag =1;
					int min = ((TH2F*)mcfile->Get(embhistname1[0]))->GetYaxis()->FindBin(PtEdge[npt]);
					int max = ((TH2F*)mcfile->Get(embhistname1[0]))->GetYaxis()->FindBin(PtEdge[npt+1])-1;
					xx[npt] = (PtEdge[npt]+PtEdge[npt+1])/2-0.1*trig;
					xxerr[npt] = (PtEdge[npt+1]-PtEdge[npt])/2;

					if(eIDmode==1){Eff[trig][nphase][npt] = new TGraphAsymmErrors((TH1F*)((TH2F*)mcfile->Get(embhistname1[nphase*5+trig+1]))->ProjectionX(Form("hRcJpsiCosThetaIn%dPt%dphase%d",trig,npt,nphase),min,max),(TH1F*)((TH2F*)mcfile->Get(embhistname1[nphase*5]))->ProjectionX(Form("hMcJpsiCosThetaIn%dPt%dphase%d",trig,npt,nphase),min,max),"N");
						if(nphase%2==0) HistEff[trig][nphase][npt]=new TH1F(Form("hEff%s%sPt%d",TrigName[trig].Data(),PhaseName[nphase].Data(),npt),"",10,-1,1);
						else HistEff[trig][nphase][npt]=new TH1F(Form("hEff%s%sPt%d",TrigName[trig].Data(),PhaseName[nphase].Data(),npt),"",10,-TMath::Pi(),TMath::Pi());
					}
					else{
						Eff[trig][nphase][npt] = new TGraphAsymmErrors((TH1F*)((TH2F*)mcfile->Get(embhistname2[nphase*5+trig+1]))->ProjectionX(Form("hRcJpsiCosThetaIn%dPt%dphase%d",trig,npt,nphase),min,max),(TH1F*)((TH2F*)mcfile->Get(embhistname2[nphase*5]))->ProjectionX(Form("hMcJpsiCosThetaIn%dPt%dphase%d",trig,npt,nphase),min,max),"N");
						if(nphase%2==0) HistEff[trig][nphase][npt]=new TH1F(Form("hEff%s%sPt%d",TrigName[trig].Data(),PhaseName[nphase].Data(),npt),"",10,-1,1);
						else HistEff[trig][nphase][npt]=new TH1F(Form("hEff%s%sPt%d",TrigName[trig].Data(),PhaseName[nphase].Data(),npt),"",10,-TMath::Pi(),TMath::Pi());
					}
					for (int nphi=0; nphi<NPHI; nphi++) {
						Double_t x,y;
						Eff[trig][nphase][npt]->GetPoint(nphi,x,y);
						Eff[trig][nphase][npt]->Print("all");
						HistEff[trig][nphase][npt]->SetBinContent(nphi+1,y);
						HistEff[trig][nphase][npt]->SetBinError(nphi+1,Eff[trig][nphase][npt]->GetErrorY(nphi));
					}

					canvas[trig][nphase][0]->cd(npt+1);
					HistEff[trig][nphase][npt]->SetMinimum(0);
					HistEff[trig][nphase][npt]->SetMaximum(0.5);
					HistEff[trig][nphase][npt]->SetMarkerStyle(20);

                   
                    
					HistEff[trig][nphase][npt]->Draw();
					HistEff[trig][nphase][npt]->GetYaxis()->SetTitle("efficiency");
					HistEff[trig][nphase][npt]->GetYaxis()->SetTitleOffset(0.9);
					HistEff[trig][nphase][npt]->GetYaxis()->SetTitleSize(.05);
				
                    ptbintext.SetNDC();
                    ptbintext.DrawLatex(0.2,0.8,Form("%1.fGeV/c < J/#psi Pt < %1.fGeV/c",PtEdge[npt],PtEdge[npt+1]));
                    

					efficiency->cd();
					HistEff[trig][nphase][npt]->Write();
					addtext(npt,trig,nphase,eIDmode,HistEff[trig][nphase][npt]);

					if(eIDmode==1)datafile[trig] = new TFile(datafilename1[trig]);
					else datafile[trig] = new TFile(datafilename2[trig]);
					HistPhasePtUL[trig][nphase] = (TH2F*)datafile[trig]->Get(datahistname[2*nphase]);
					HistPhasePtLI[trig][nphase] = (TH2F*)datafile[trig]->Get(datahistname[2*nphase+1]);
					DefaultTH2F[trig][nphase][0] = (TH2F*)DefaultFile[trig][eIDmode-1]->Get(datahistname[2*nphase]);
					DefaultTH2F[trig][nphase][1] = (TH2F*)DefaultFile[trig][eIDmode-1]->Get(datahistname[2*nphase+1]);

					HistPhase1[trig][nphase][npt] = (TH1F*)HistPhasePtUL[trig][nphase]->ProjectionX(Form("1hRawJpsiCosThetaIn%dPt%d",trig,npt),min,max);
					HistPhase2[trig][nphase][npt] = (TH1F*)HistPhasePtLI[trig][nphase]->ProjectionX(Form("2hRawJpsiCosThetaIn%dPt%d",trig,npt),min,max);
					HistPhase0[trig][nphase][npt]=(TH1F*)HistPhase1[trig][nphase][npt]->Clone();
					HistPhase0[trig][nphase][npt]->SetName(Form("0hRawJpsiIn%dPt%dPhase%d%d",trig,npt,nphase,eIDmode));
					HistPhase0[trig][nphase][npt]->Add(HistPhase2[trig][nphase][npt],-1);
					DefaultTH1F[trig][nphase][npt][1] = (TH1F*)DefaultTH2F[trig][nphase][0]->ProjectionX(Form("default1hRawJpsiCosThetaIn%dPt%d",trig,npt),min,max);
					DefaultTH1F[trig][nphase][npt][2] = (TH1F*)DefaultTH2F[trig][nphase][1]->ProjectionX(Form("default2hRawJpsiCosThetaIn%dPt%d",trig,npt),min,max);
					DefaultTH1F[trig][nphase][npt][0]=(TH1F*)DefaultTH1F[trig][nphase][npt][1]->Clone();
					DefaultTH1F[trig][nphase][npt][0]->SetName(Form("default0hRawJpsiIn%dPt%dPhase%d%d",trig,npt,nphase,eIDmode));
					DefaultTH1F[trig][nphase][npt][0]->Add(DefaultTH1F[trig][nphase][npt][2],-1);

					HistPhase1[trig][nphase][npt]->SetMarkerStyle(29);
					HistPhase1[trig][nphase][npt]->SetMarkerColor(4);
					HistPhase1[trig][nphase][npt]->SetLineColor(4);
					HistPhase1[trig][nphase][npt]->SetMarkerSize(2.3);
					DefaultTH1F[trig][nphase][npt][1]->SetMarkerStyle(20);
					DefaultTH1F[trig][nphase][npt][1]->SetMarkerColor(4);
					DefaultTH1F[trig][nphase][npt][1]->SetLineColor(4);
					DefaultTH1F[trig][nphase][npt][1]->SetMarkerSize(0.9);
					HistPhase2[trig][nphase][npt]->SetMarkerStyle(22);
					HistPhase2[trig][nphase][npt]->SetMarkerColor(1);
					HistPhase2[trig][nphase][npt]->SetLineColor(1);
					HistPhase2[trig][nphase][npt]->SetMarkerSize(1.8);
					DefaultTH1F[trig][nphase][npt][2]->SetMarkerStyle(21);
					DefaultTH1F[trig][nphase][npt][2]->SetMarkerColor(1);
					DefaultTH1F[trig][nphase][npt][2]->SetLineColor(1);
					DefaultTH1F[trig][nphase][npt][2]->SetMarkerSize(0.9);
					HistPhase0[trig][nphase][npt]->SetMarkerStyle(20);
					HistPhase0[trig][nphase][npt]->SetMarkerColor(2);
					HistPhase0[trig][nphase][npt]->SetLineColor(2);
					HistPhase0[trig][nphase][npt]->SetMarkerSize(1.1);
					DefaultTH1F[trig][nphase][npt][0]->SetMarkerStyle(22);
					DefaultTH1F[trig][nphase][npt][0]->SetMarkerColor(2);
					DefaultTH1F[trig][nphase][npt][0]->SetLineColor(2);
					DefaultTH1F[trig][nphase][npt][0]->SetMarkerSize(0.9);	

					canvas[trig][nphase][1]->cd(npt+1);
					if(npt!=5 )	HistPhase1[trig][nphase][npt]->SetMinimum(HistPhase2[trig][nphase][npt]->GetMinimum()-25);
					elseif(npt==5 )	HistPhase1[trig][nphase][npt]->SetMinimum(HistPhase2[trig][nphase][npt]->GetMinimum()-10);
                    HistPhase1[trig][nphase][npt]->Draw();
					HistPhase2[trig][nphase][npt]->Draw("same");
					HistPhase0[trig][nphase][npt]->Draw("same");
					addtext(npt,trig,nphase,eIDmode,HistPhase1[trig][nphase][npt]);
					if(npt==3){
                        if(trig==1 && eIDmode==2){
                            HistPhase1[trig][nphase][npt]->SetMinimum(HistPhase2[trig][nphase][npt]->GetMinimum()-3);
                            HistPhase1[trig][nphase][npt]->SetMaximum(35);
                        }
                        if(trig==3 && eIDmode==1 && nphase==1){
                            HistPhase1[trig][nphase][npt]->SetMaximum(400);
                            HistPhase1[trig][nphase][npt]->SetMinimum(0);
                        }
						STARpreliminary.SetNDC();
						STARpreliminary.SetTextFont(2);
						STARpreliminary.SetTextSize(0.08);
						STARpreliminary.DrawLatex(0.2,0.7,"STAR preliminary");
					}
					ptbintext.SetNDC();
					ptbintext.DrawLatex(0.2,0.8,Form("%1.fGeV/c < J/#psi Pt < %1.fGeV/c",PtEdge[npt],PtEdge[npt+1]));

                    
                    if (npt==4) {
                        if(trig==1 && eIDmode==2){
                            HistPhase1[trig][nphase][npt]->SetMinimum(0);
                            HistPhase1[trig][nphase][npt]->SetMaximum(16);
                        }
                        if (trig==3 && eIDmode==2) {
                            HistPhase1[trig][nphase][npt]->SetMinimum(HistPhase2[trig][nphase][npt]->GetMinimum()-5);
                        }
                        if (trig==3 && eIDmode==1 && nphase==1) {
                            HistPhase1[trig][nphase][npt]->SetMaximum(180);
                            HistPhase1[trig][nphase][npt]->SetMinimum(0);
                        }
                    }
                    
					if(npt==5){
                        if (trig==3 && eIDmode==1 && nphase==1) {
                            HistPhase1[trig][nphase][npt]->SetMaximum(85);
                            HistPhase1[trig][nphase][npt]->SetMinimum(0);
						}
                        if (trig==1 && eIDmode==2 ){
                            if(nphase==0) {
                                HistPhase1[trig][nphase][npt]->SetMinimum(0);
                                HistPhase1[trig][nphase][npt]->SetMaximum(10);
                            }
                            if(nphase==1){
                                HistPhase1[trig][nphase][npt]->SetMinimum(0);
                                HistPhase1[trig][nphase][npt]->SetMaximum(10);
                            }
                        }
                        if (trig==3 && eIDmode==2) {
                            HistPhase1[trig][nphase][npt]->SetMinimum(HistPhase2[trig][nphase][npt]->GetMinimum()-5);
                            HistPhase1[trig][nphase][npt]->SetMaximum(35);
                        }
						//if(nphase==1) HistPhase1[trig][nphase][npt]->SetMaximum(100);
						rawdatalegend = new TLegend(0.45,0.5,0.8,0.8);
						rawdatalegend->SetTextSize(0.06);
				//		rawdatalegend->AddEntry(,Form("%1.fGeV/c < Jpsi Pt < %1.fGeV/c",PtEdge[npt],PtEdge[npt+1]));
						rawdatalegend->AddEntry(HistPhase1[trig][nphase][npt],"Unlike-Sign","p");
						rawdatalegend->AddEntry(HistPhase2[trig][nphase][npt],"Like-Sign","p");
						rawdatalegend->AddEntry(HistPhase0[trig][nphase][npt],"J/#psi Signal","p"); 
					//	rawdatalegend->Draw("same");
					}
					canvas[trig][nphase][3]->cd(npt+1);
					//		HistPhase1[trig][nphase][npt]->SetMinimum(HistPhase2[trig][nphase][npt]->GetMinimum()-25);
					HistPhase1[trig][nphase][npt]->Draw();
					HistPhase2[trig][nphase][npt]->Draw("same");
					HistPhase0[trig][nphase][npt]->Draw("same");

					DefaultTH1F[trig][nphase][npt][0]->Draw("same");
					DefaultTH1F[trig][nphase][npt][1]->Draw("same");
					DefaultTH1F[trig][nphase][npt][2]->Draw("same");
					//			HistPhase1[trig][nphase][npt]->SetTitle(Form("%s in %s %deID",TrigName[trig].Data(),FrameName[nframe].Data(),eIDmode));		
					rawdatalegend = new TLegend(0.7,0.6,0.9,0.9);
					rawdatalegend->AddEntry(HistPhase0[trig][nphase][npt],"revised signal","p");
					rawdatalegend->AddEntry(HistPhase1[trig][nphase][npt],"revised unlike","p");
					rawdatalegend->AddEntry(HistPhase2[trig][nphase][npt],"revised like","p");
					rawdatalegend->AddEntry(DefaultTH1F[trig][nphase][npt][0],"default signal","p");
					rawdatalegend->AddEntry(DefaultTH1F[trig][nphase][npt][1],"default unlike","p");
					rawdatalegend->AddEntry(DefaultTH1F[trig][nphase][npt][2],"default like","p");
					rawdatalegend->Draw("same");



					addtext(npt,trig,nphase,eIDmode,HistPhase1[trig][nphase][npt]);


					HistPhase[trig][nphase][npt]=(TH1F*)HistPhase0[trig][nphase][npt]->Clone();
					HistPhase[trig][nphase][npt]->SetName(Form("HistJpsiPhase_%d_%d_%d",trig,npt,nphase));
					HistPhase[trig][nphase][npt]->Divide(HistEff[trig][nphase][npt]);
					if(HistPhase[trig][nphase][npt]->Integral()>0){
						HistPhase[trig][nphase][npt]->Scale(1/(HistPhase[trig][nphase][npt]->Integral()));
						corrected->cd();
						HistPhase[trig][nphase][npt]->Write();
					}
					else flag=0;

					canvas[trig][nphase][2]->cd(npt+1);
					HistPhase[trig][nphase][npt]->Draw();
					addtext(npt,trig,nphase,eIDmode,HistPhase[trig][nphase][npt]);
					gStyle->SetOptFit(1);
					gStyle->SetOptStat(0);

					func1 = new TF1("func1","[0]*(1+[1]*x*x)",-1,1);
					func2 = new TF1("func2","[2]*(1+2*[1]*TMath::Cos(2*x)/(3+[0]))",-TMath::Pi(),TMath::Pi());

					if (nphase%2==0) {
						func1->SetParameters(0.1,-0.5);
						int empty1=0;
						for(int nbin=1;nbin<=10;nbin++) {
							if(HistPhase[trig][nphase][npt]->GetBinContent(nbin)<=0) empty1++;
						}
						if (empty1>=6) flag=0;
						HistPhase[trig][nphase][npt]->Fit("func1","R");
						//                   if(func1->GetChisquare()/func1->GetNDF()>3) flag=0;
						if(flag==1) {
							LambdaTheta[npt] = func1->GetParameter(1);
							LambdaThetaErr[npt] = func1->GetParError(1);
						} else {
							LambdaTheta[npt] = -100;
							LambdaThetaErr[npt] = 0;
						}
					} else {
						func2->SetParameters(LambdaTheta[npt],0.05,0.1);
						func2->FixParameter(0,LambdaTheta[npt]);

						int empty2 = 0;
						for(int nbin=1;nbin<=10;nbin++) {
							if(HistPhase[trig][nphase][npt]->GetBinContent(nbin)<=0) empty2++;
						}
						if (empty2>=6) flag=0;
						HistPhase[trig][nphase][npt]->Fit("func2","R");
						//                 if(func2->GetChisquare()/func2->GetNDF()>3) flag=0;
						if(flag==1) {
							LambdaPhi[npt] = func2->GetParameter(1);
							LambdaPhiErr[npt] = func2->GetParError(1);
						} else {
							LambdaPhi[npt] = -100;
							LambdaPhiErr[npt] = 0;
						}
						if(LambdaPhi[npt]>-99&&LambdaPhi[npt]>-99) {
							LambdaInv[npt] = (LambdaTheta[npt]+3*LambdaPhi[npt])/(1-LambdaPhi[npt]);
							LambdaInvErr[npt] = TMath::Sqrt(TMath::Power(LambdaThetaErr[npt]/(1-LambdaPhi[npt]),2)+TMath::Power(LambdaPhiErr[npt]*(3+LambdaTheta[npt])/((1-LambdaPhi[npt])*(1-LambdaPhi[npt])),2));
						} else {
							LambdaInv[npt] = -100;
							LambdaInvErr[npt] = 0;
						}
					}
				}

				if(eIDmode==1){
					canvas[trig][nphase][0]->SaveAs(Form("figures/%s_%s_eff1eID.pdf",TrigName[trig].Data(),PhaseName[nphase].Data()));
					canvas[trig][nphase][1]->SaveAs(Form("figures/%s_%s_raw1eID.pdf",TrigName[trig].Data(),PhaseName[nphase].Data()));
					canvas[trig][nphase][2]->SaveAs(Form("figures/%s_%s_fit1eID.pdf",TrigName[trig].Data(),PhaseName[nphase].Data()));
					canvas[trig][nphase][3]->SaveAs(Form("default_%s_%s_raw1eID.pdf",TrigName[trig].Data(),PhaseName[nphase].Data()));
				}

				if(eIDmode==2){
					canvas[trig][nphase][0]->SaveAs(Form("figures/%s_%s_eff2eID.pdf",TrigName[trig].Data(),PhaseName[nphase].Data()));
					canvas[trig][nphase][1]->SaveAs(Form("figures/%s_%s_raw2eID.pdf",TrigName[trig].Data(),PhaseName[nphase].Data()));
					canvas[trig][nphase][2]->SaveAs(Form("figures/%s_%s_fit2eID.pdf",TrigName[trig].Data(),PhaseName[nphase].Data()));
					canvas[trig][nphase][3]->SaveAs(Form("default_%s_%s_raw2eID.pdf",TrigName[trig].Data(),PhaseName[nphase].Data()));
				}

				if(nphase==1||nphase==3) {
					Lambda_graph[trig][nframe][0]=new TGraphErrors(NPT,xx,LambdaTheta,xxerr,LambdaThetaErr);
					Lambda_graph[trig][nframe][0]->GetYaxis()->SetRangeUser(-2,2);
					Lambda_graph[trig][nframe][0]->GetXaxis()->SetRangeUser(0,14);
					Lambda_graph[trig][nframe][1] = new TGraphErrors(NPT,xx,LambdaPhi,xxerr,LambdaPhiErr);
					Lambda_graph[trig][nframe][1]->GetYaxis()->SetRangeUser(-2,2);
					Lambda_graph[trig][nframe][1]->GetXaxis()->SetRangeUser(0,14);
					Lambda_graph[trig][nframe][2] = new TGraphErrors(NPT,xx,LambdaInv,xxerr,LambdaInvErr);
					Lambda_graph[trig][nframe][2]->GetYaxis()->SetRangeUser(-5,5);
					Lambda_graph[trig][nframe][2]->GetXaxis()->SetRangeUser(0,14);
				}
			}
		}

		TLegend *legend=new TLegend(0.7,0.7,0.9,0.9);
		legend->AddEntry(Lambda_graph[1][0][0],"HT0","lep");
		legend->AddEntry(Lambda_graph[2][0][0],"HT1","lep");
		legend->AddEntry(Lambda_graph[3][0][0],"HT2","lep");
		legend->AddEntry(Lambda_graph[0][0][0],"MB","lep");

		if(eIDmode==1)Lambda = new TFile("Lambda1eID.root","RECREATE"); 
		if(eIDmode==2)Lambda = new TFile("Lambda2eID.root","RECREATE");
		if(eIDmode==1){
			ncanvas= new TCanvas("Lambda1eID","",1200,800);
		}
		if(eIDmode==2){
			ncanvas= new TCanvas("Lambda2eID","",1200,800);
		}
		ncanvas->Divide(3,2);
		for(int trig=0;trig<NTRIG;trig++) {
			for(int nframe=0;nframe<2;nframe++) {
				for(int nvar=0;nvar<3;nvar++) {
					ncanvas->cd(nframe*3+nvar+1);
					Lambda_graph[trig][nframe][nvar]->SetTitle(Form("%s in %s frame %deID",VarName[nvar].Data(),FrameName[nframe].Data(),eIDmode));
					Lambda_graph[trig][nframe][nvar]->SetName(Form("%s_%s_%d",VarName[nvar].Data(),FrameName[nframe].Data(),trig));
					Lambda_graph[trig][nframe][nvar]->GetYaxis()->SetTitle(VarName[nvar].Data());
					Lambda_graph[trig][nframe][nvar]->GetXaxis()->SetTitle("J/#psi p_{T} (GeV/c)");

					Lambda_graph[trig][nframe][nvar]->GetYaxis()->SetTitleSize(0.045);
					Lambda_graph[trig][nframe][nvar]->GetXaxis()->SetTitleSize(0.045);
					Lambda_graph[trig][nframe][nvar]->Write();

					if(trig==0) {
						Lambda_graph[trig][nframe][nvar]->Draw("AP");
						legend->Draw();
						if(nframe==0&&nvar==0) run09->Draw("SAMEPE");
					}
					Lambda_graph[trig][nframe][nvar]->SetMarkerStyle(20);
					Lambda_graph[trig][nframe][nvar]->SetMarkerSize(1);
					Lambda_graph[trig][nframe][nvar]->SetMarkerColor(1+trig);
					Lambda_graph[trig][nframe][nvar]->SetLineColor(1+trig);
					Lambda_graph[trig][nframe][nvar]->Draw("SAMEPE");
				}
			}
		}
		Lambda->Close();
		ncanvas->SaveAs(Form("figures/Lambda%deID.pdf",eIDmode));

		//		TCanvas *Run2009 = new TCanvas("Run2009","Run2009",1200,800);
		//		TLegend * legend2009 = new TLegend(0.7,0.7,0.9,0.9);
		//		Run2009->cd();
		//		TLatex tex2009;

		//		Lambda_graph[1][0][0]->SetMarkerStyle(29);
		//		Lambda_graph[1][0][0]->SetMarkerSize(3);
		//		Lambda_graph[1][0][0]->SetMarkerColor(1);
		//		Lambda_graph[1][0][0]->SetLineColor(1);
		//		Lambda_graph[1][0][0]->SetTitle("#lambda_{#theta} as a function of J/#psi P_{t} in p+p collision at #surds = 200GeV");
		//		Lambda_graph[1][0][0]->Draw("ap");

		//		run09->SetMarkerStyle(30);
		//		run09->SetMarkerSize(3);
		//		run09->Draw("samepe");
		//		tex2009.SetTextAlign(11);
		//		tex2009.SetNDC(kTRUE);
		//		tex2009.DrawLatex(0.2,0.8,"p+p collision #surds = 200GeV");


		//		legend2009->AddEntry(run09,"STAR RUN9","p");
		//		legend2009->AddEntry(Lambda_graph[1][0][0],"STAR RUN12","p");


		//		legend2009->Draw();

		//		Run2009->SaveAs("figures/Run2009.pdf");
		/*		if(eIDmode==1){		
				TFile *lambdatheta1eID = new TFile("lambdatheta1eID.root","recreate");
				lambdatheta1eID->cd();
				Lambda_graph[1][0][0]->Write();
				}
				if(eIDmode==2){
				TFile *lambdatheta2eID = new TFile("lambdatheta2eID.root","recreate");
				lambdatheta2eID->cd();
				Lambda_graph[1][0][0]->Write();
				}
		 */
	}

	effmacro();
}

void addtext(int n, int k, int l, int eID, TH1F *hist){
	TString PtRange[] = {"0<p_{T}<2 GeV/c","2<p_{T}<3 GeV/c","3<p_{T}<4 GeV/c","4<p_{T}<6 GeV/c","6<p_{T}<8 GeV/c","8<p_{T}<14 GeV/c"};
	TString triggername[] = {"MB ","HT0","HT1","HT2"};
	TLatex  tex;
	tex.SetTextAlign(11);
	tex.SetNDC(kTRUE);
	//tex.DrawText(0.3,0.8,PtRange[n]);
	//tex.DrawClone("SAME");

	hist->GetXaxis()->SetTitleOffset(0.6);
	hist->GetXaxis()->SetTitleSize(0.06);
	if (l==0) {
		hist->SetTitle(triggername[k].Append(" in HX frame ")+PtRange[n]+Form("  %deID",eID));
		hist->GetXaxis()->SetTitle("cos#theta");
	}
	else if(l==1){
		hist->SetTitle(triggername[k].Append(" in HX frame")+PtRange[n]+Form("  %deID",eID));
		hist->GetXaxis()->SetTitle("#varphi");
	}
	else if(l==2){
		hist->SetTitle(triggername[k].Append(" in CS frame")+PtRange[n]+Form("  %deID",eID));
		hist->GetXaxis()->SetTitle("cos#theta");
	}
	else{
		hist->SetTitle(triggername[k].Append(" in CS frame")+PtRange[n]+Form("  %deID",eID));
		hist->GetXaxis()->SetTitle("#varphi");
	}
}

void effmacro(){
	TFile *mcfile = new TFile("rootfile/OutFile_Cent_0_9.root","read");
	TFile *mcfile1 = new TFile("~/jpsi/test20160210_Barbara/rootfile/OutFile_cent_0_9.root","read");

	TString embhistname2[]={"hMcJpsiCosThetaPt","hMBJpsiCosThetaPt","hHt0JpsiCosThetaPt","hHt1JpsiCosThetaPt","hHt2JpsiCosThetaPt","hMcJpsiPhiPt","hMBJpsiPhiPt","hHtJpsiPhiPt","hHt1JpsiPhiPt","hHt2JpsiPhiPt","hMcJpsiCosThetaPtCS","hMBJpsiCosThetaPtCS","hHt0JpsiCosThetaPtCS","hHt1JpsiCosThetaPtCS","hHt2JpsiCosThetaPtCS","hMcJpsiPhiPtCS","hMBJpsiPhiPtCS","hHt0JpsiPhiPtCS","hHt1JpsiPhiPtCS","hHt2JpsiPhiPtCS"};
	TString embhistname1[NPHASE*(NTRIG+1)]={"hMcJpsiCosThetaPt","hMBJpsiCosThetaPt1","hHt0JpsiCosThetaPt1","hHt1JpsiCosThetaPt1","hHt2JpsiCosThetaPt1","hMcJpsiPhiPt","hMBJpsiPhiPt1","hHt0JpsiPhiPt1","hHt1JpsiPhiPt1","hHt2JpsiPhiPt1","hMcJpsiCosThetaPtCS","hMBJpsiCosThetaPtCS1","hHt0JpsiCosThetaPtCS1","hHt1JpsiCosThetaPtCS1","hHt2JpsiCosThetaPtCS1","hMcJpsiPhiPtCS","hMBJpsiPhiPtCS1","hHt0JpsiPhiPtCS1","hHt1JpsiPhiPtCS1","hHt2JpsiPhiPtCS1"};
	TH1F *hN[NTRIG][NPT][EID][2];
	TH1F *hD[NTRIG][NPT][EID][2];
	TCanvas *canvas[4][EID];

	TH1F *ratio[NPT];
	TLegend *legend[5];
	gStyle->SetOptStat(0);
	for(int eID=1;eID<3;eID++){
		for(int trig=0;trig<NTRIG;trig++){
			canvas[trig][eID-1] = new TCanvas(Form("canvaseff_%d_eID%d",trig,eID),Form("canvaseff_%d",trig),1200,800);
			canvas[trig][eID-1]->Divide(3,2);
			int nphase = 0;
			for(int npt=0;npt<NPT;npt++){
				canvas[trig][eID-1]->cd(npt+1);
				int min = ((TH2F*)mcfile->Get(embhistname1[0]))->GetYaxis()->FindBin(PtEdge[npt]);
				int max = ((TH2F*)mcfile->Get(embhistname1[0]))->GetYaxis()->FindBin(PtEdge[npt+1])-1;
				if(eID==1){				
					hN[trig][npt][eID-1][0]=(TH1F*)((TH2F*)mcfile->Get(embhistname1[nphase*5+trig+1]))->ProjectionX(Form("hRcJpsiCosThetaIn%dPt%dphase%deID%d",trig,npt,nphase,eID),min,max);
					hN[trig][npt][eID-1][0]->SetName(Form("hN0_npt%d_%d",npt,eID));
					hD[trig][npt][eID-1][0]=(TH1F*)((TH2F*)mcfile->Get(embhistname1[nphase*5]))->ProjectionX(Form("hMcJpsiCosThetaIn%dPt%dphase%deID%d",trig,npt,nphase,eID),min,max);
					hD[trig][npt][eID-1][0]->SetName(Form("hD0_npt%d_%d",npt,eID));
					hN[trig][npt][eID-1][1]=(TH1F*)((TH2F*)mcfile1->Get(embhistname1[nphase*5+trig+1]))->ProjectionX(Form("hRcJpsiCosThetaIn%dPt%dphase%deID%d",trig,npt,nphase,eID),min,max);
					hN[trig][npt][eID-1][1]->SetName(Form("hN1_npt%d_d",npt,eID));
					hD[trig][npt][eID-1][1]=(TH1F*)((TH2F*)mcfile1->Get(embhistname1[nphase*5]))->ProjectionX(Form("hMcJpsiCosThetaIn%dPt%dphase%deID%d",trig,npt,nphase,eID),min,max);
					hD[trig][npt][eID-1][1]->SetName(Form("hD1_npt%d_d",npt,eID));
				}			
				else{
					hN[trig][npt][eID-1][0]=(TH1F*)((TH2F*)mcfile->Get(embhistname2[nphase*5+trig+1]))->ProjectionX(Form("hRcJpsiCosThetaIn%dPt%dphase%deID%d",trig,npt,nphase,eID),min,max);
					hN[trig][npt][eID-1][0]->SetName(Form("hN0_npt%d_%d",npt,eID));
					hD[trig][npt][eID-1][0]=(TH1F*)((TH2F*)mcfile->Get(embhistname2[nphase*5]))->ProjectionX(Form("hMcJpsiCosThetaIn%dPt%dphase%deID%d",trig,npt,nphase,eID),min,max);
					hD[trig][npt][eID-1][0]->SetName(Form("hD0_npt%d_%d",npt,eID));
					hN[trig][npt][eID-1][1]=(TH1F*)((TH2F*)mcfile1->Get(embhistname2[nphase*5+trig+1]))->ProjectionX(Form("hRcJpsiCosThetaIn%dPt%dphase%deID%d",trig,npt,nphase,eID),min,max);
					hN[trig][npt][eID-1][1]->SetName(Form("hN1_npt%d_d",npt,eID));
					hD[trig][npt][eID-1][1]=(TH1F*)((TH2F*)mcfile1->Get(embhistname2[nphase*5]))->ProjectionX(Form("hMcJpsiCosThetaIn%dPt%dphase%deID%d",trig,npt,nphase,eID),min,max);
					hD[trig][npt][eID-1][1]->SetName(Form("hD1_npt%d_d",npt,eID));
				}
				hD[trig][npt][eID-1][1]->SetMinimum(0);
				hD[trig][npt][eID-1][0]->SetMarkerStyle(24);
				hD[trig][npt][eID-1][0]->SetMarkerSize(1.3);
				hN[trig][npt][eID-1][0]->SetMarkerStyle(24);
				hN[trig][npt][eID-1][0]->SetMarkerSize(1.3);
				hD[trig][npt][eID-1][1]->SetMarkerColor(kRed);
				hD[trig][npt][eID-1][1]->SetLineColor(kRed);
				hN[trig][npt][eID-1][1]->SetMarkerColor(kBlack);
				hN[trig][npt][eID-1][1]->SetLineColor(kBlack);

				hD[trig][npt][eID-1][1]->Draw();
				if(eID==1)	hD[trig][npt][eID-1][1]->SetTitle(Form("%s %1.fGeV/c < Jpsi Pt < %1.fGeV/c %s",TrigName[trig].Data(),PtEdge[npt],PtEdge[npt+1],"1eID"));
				else 		hD[trig][npt][eID-1][1]->SetTitle(Form("%s %1.fGeV/c < Jpsi Pt < %1.fGeV/c %s",TrigName[trig].Data(),PtEdge[npt],PtEdge[npt+1],"2eID"));
				hN[trig][npt][eID-1][1]->Draw("same");
				hD[trig][npt][eID-1][0]->Draw("same");
				hN[trig][npt][eID-1][0]->Draw("same");
				hD[trig][npt][eID-1][1]->SetMarkerStyle(20);
				hN[trig][npt][eID-1][1]->SetMarkerStyle(20);
				hD[trig][npt][eID-1][0]->SetMarkerColor(kRed);
				hD[trig][npt][eID-1][0]->SetLineColor(kRed);
				hN[trig][npt][eID-1][0]->SetMarkerColor(kBlack);
				hN[trig][npt][eID-1][0]->SetLineColor(kBlack);
				legend[npt] = new TLegend(0.7,0.5,0.9,0.7);
				legend[npt]->AddEntry(hN[trig][npt][eID-1][0],"revised","lep");
				legend[npt]->AddEntry(hN[trig][npt][eID-1][1],"default","lep");
				legend[npt]->Draw();	
			}
			canvas[trig][eID-1]->SaveAs(Form("canvaseff_%d_eID%d.pdf",trig,eID));
		}
	}
}
