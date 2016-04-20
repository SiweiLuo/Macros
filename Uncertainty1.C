#define NTRIG 4
#define NPT 6
#define NVAR 3
#define NPHASE 4 
#define NFRAME 2
#define EID 2
#define PLOT 6

Double_t PtEdge[NPT+1] = {0,2,3,4,6,8,14};

void Uncertainty1(){
	TFile *lambda1;
	TFile *lambda2;
	TFile *efficiency1;
	TFile *efficiency2;
	TFile *correcteddata1;
	TFile *correcteddata2;
	TFile *uncertainty;
	TFile *rawdata[NTRIG][EID][2];

	TH1F *histeff[NTRIG][NPT][2];

	TGraphErrors *graph[4][2][3][3];
	TString TrigName[NTRIG]={"MB","HT0","HT1","HT2"};
	TString PhaseName[NPHASE]={"HXTheta","HXPhi","CSTheta","CSPhi"};
	TString FrameName[NFRAME]={"HX","CS"};
	TString VarName[NVAR]={"lambda_theta","lambda_phi","lambda_inv"};

	TH1F *ratio[NTRIG][NPT];
	TCanvas *canvaseff[NTRIG];
	TLegend *legendeff;

	TH1F *histfit[NTRIG][NPT][2];
	TLegend *legendfit;

	TCanvas *canvas[PLOT][NTRIG][EID];
	TLegend *legend;

	TCanvas *canvasfit = new TCanvas("canvasfit","canvasfit",1200,800);
	TLatex text;
	//				TF1 *fun1 = new TF1("fun1","[0]*(1+[1]*x*x)",-1,1);
	//				fun1->SetParameters(0.1,0.5);
	//				fun1->SetLineColor(1);
	//				TF1 *fun2 = new TF1("fun2","[0]*(1+[1]*x*x)",-1,1);
	//				fun1->SetParameters(0.1,0.5);
	//				fun2->SetLineColor(2);

	for(int eID=1;eID<3;eID++){
		for(int trig=0;trig<4;trig++){
			if(trig==0)rawdata[trig][eID-1][0] = new TFile(Form("rootfile/mb_%dTrkPid_0.ana.root",eID));
			else rawdata[trig][eID-1][0] = new TFile(Form("rootfile/ht%d_trg%d_%dTrkPid_0.ana.root",trig-1,trig+2,eID));
			if(trig==0)rawdata[trig][eID-1][1] = new TFile(Form("../test20160210_0/rootfile/mb_%dTrkPid_0.ana.root",eID));
			else rawdata[trig][eID-1][1] = new TFile(Form("../test20160210_0/rootfile/ht%d_trg%d_%dTrkPid_0.ana.root",trig-1,trig+2,eID));
		}

		if(eID==1){
			lambda1 = new TFile("Lambda1eID.root");
			lambda2 = new TFile("~/jpsi/test20160210_0/Lambda1eID.root");
			efficiency1 = new TFile("efficiency1eID.root","read");
			efficiency2 = new TFile("~/jpsi/test20160210_0/efficiency1eID.root","read");
			correcteddata1 = new TFile("corrected1eID.root","read");
			correcteddata2 = new TFile("~/jpsi/test20160210_0/corrected1eID.root","read");
		}
		if(eID==2){
			lambda1 = new TFile("Lambda2eID.root");
			lambda2 = new TFile("~/jpsi/test20160210_0/Lambda2eID.root");
			efficiency1 = new TFile("efficiency2eID.root","read");
			efficiency2 = new TFile("~/jpsi/test20160210_0/efficiency2eID.root","read");
			correcteddata1 = new TFile("corrected2eID.root","read");
			correcteddata2 = new TFile("~/jpsi/test20160210_0/corrected2eID.root","read");
		}
		uncertainty = new TFile(Form("uncertainty%deID.root",eID),"recreate");

		Double_t xx[NPT],yy[NPT],xerr[NPT],yerr[NPT];
		canvas[0][0][eID] = new TCanvas(Form("uncertainty%deID",eID),Form("uncertainty%deID",eID),1200,800);
		canvas[0][0][eID]->Divide(3,2);		
		for(int ntrig=0;ntrig<4;ntrig++){
			//canvas[0][trig][eID] = new TCanvas(Form("raw_%s_%deID",TrigName[trig].Data(),eID),Form("raw_%s_%deID",TrigName[trig].Data(),eID),1200,800);
			//canvas[0][trig][eID]->Divide(3,2);
			for(int nframe=0;nframe<2;nframe++){
				for(int nvar=0;nvar<3;nvar++){
					graph[ntrig][nframe][nvar][0] = (TGraphErrors*)lambda1->Get(Form("%s_%s_%d",VarName[nvar].Data(),FrameName[nframe].Data(),ntrig));
					graph[ntrig][nframe][nvar][1] = (TGraphErrors*)lambda2->Get(Form("%s_%s_%d",VarName[nvar].Data(),FrameName[nframe].Data(),ntrig));

					graph[ntrig][nframe][nvar][0]->SetName(Form("%s_%s_%d_0",VarName[nvar].Data(),FrameName[nframe].Data(),ntrig));
					graph[ntrig][nframe][nvar][1]->SetName(Form("%s_%s_%d_1",VarName[nvar].Data(),FrameName[nframe].Data(),ntrig));

					for(int npoint=0;npoint<5;npoint++){
						graph[ntrig][nframe][nvar][1]->GetPoint(npoint,xx[npoint],yy[npoint]);
						xx[npoint] = xx[npoint]+0.3;
						graph[ntrig][nframe][nvar][1]->SetPoint(npoint,xx[npoint],yy[npoint]);
					}
					graph[ntrig][nframe][nvar][0]->Write();
					graph[ntrig][nframe][nvar][1]->Write();
				}
			}
		}

		//		TLegend *legend1;
		for(int i=1;i<7;i++){
			int nframe =0;
			if(i>3) nframe=1;
			canvas[0][0][eID]->cd(i);
			legend = new TLegend(0.7,0.5,1,0.9);
			graph[0][nframe][(i-1)%3][0]->Draw("ap");
			graph[0][nframe][(i-1)%3][0]->SetMarkerColor(6);
			graph[0][nframe][(i-1)%3][0]->SetMarkerStyle(30);		
			graph[0][nframe][(i-1)%3][0]->SetLineColor(6);	
			graph[0][nframe][(i-1)%3][0]->SetMarkerSize(1.3);
			graph[0][nframe][(i-1)%3][1]->Draw("samepe");
			graph[0][nframe][(i-1)%3][1]->SetMarkerColor(6);
			graph[0][nframe][(i-1)%3][1]->SetLineColor(6);
			graph[0][nframe][(i-1)%3][1]->SetMarkerStyle(29);
			graph[0][nframe][(i-1)%3][1]->SetMarkerSize(1.2);
			for(int ntrig=1;ntrig<4;ntrig++) {
				graph[ntrig][nframe][(i-1)%3][0]->Draw("samepe");
				graph[ntrig][nframe][(i-1)%3][0]->SetMarkerColor(ntrig);
				graph[ntrig][nframe][(i-1)%3][0]->SetLineColor(ntrig);
				graph[ntrig][nframe][(i-1)%3][0]->SetMarkerStyle(23+ntrig);
				graph[ntrig][nframe][(i-1)%3][0]->SetMarkerSize(1.3);
				graph[ntrig][nframe][(i-1)%3][1]->Draw("samepe");
				graph[ntrig][nframe][(i-1)%3][1]->SetMarkerColor(ntrig);
				graph[ntrig][nframe][(i-1)%3][1]->SetLineColor(ntrig);
				graph[ntrig][nframe][(i-1)%3][1]->SetMarkerStyle(19+ntrig);
				graph[ntrig][nframe][(i-1)%3][1]->SetMarkerSize(1.2);
			}

			for(int trig=0;trig<4;trig++){
				if(trig==0){
					legend->AddEntry(graph[trig][0][0][0],"MB revised","p");
					legend->AddEntry(graph[trig][0][0][1],"MB default","p");
				}
				else{
					legend->AddEntry(graph[trig][0][0][0],Form("HT%d revised",trig-1),"p");
					legend->AddEntry(graph[trig][0][0][1],Form("HT%d default",trig-1),"p");
				}
			}
			legend->Draw();
		}
		canvas[0][0][eID]->SaveAs(Form("uncertainty%deID.pdf",eID));

		for(int trig=0;trig<4;trig++){
			canvas[1][trig][eID] = new TCanvas(Form("eff_%s_%deID",TrigName[trig].Data(),eID),Form("eff_%s_%deID",TrigName[trig].Data(),eID),1200,800);
			canvas[1][trig][eID]->Divide(3,2);

			canvas[2][trig][eID] = new TCanvas(Form("ratio_%s_%deID",TrigName[trig].Data(),eID),Form("ratio_%s_%deID",TrigName[trig].Data(),eID),1200,800);
			canvas[2][trig][eID]->Divide(3,2);
			for(int npt=0;npt<NPT;npt++){	
				canvas[1][trig][eID]->cd(npt+1);				
				cout<<"npt ==========="<<npt<<endl;
					cout<<"string==========="<<Form("hEff%sHXThetaPt%d",TrigName[trig].Data(),npt)<<"   pointer===="<<(TH1F*)efficiency1->Get(Form("hEff%sHXThetaPt%d",TrigName[trig].Data(),npt))<<endl;	
				histeff[trig][npt][0] = (TH1F*)efficiency1->Get(Form("hEff%sHXThetaPt%d",TrigName[trig].Data(),npt));
				histeff[trig][npt][1] = (TH1F*)efficiency2->Get(Form("hEff%sHXThetaPt%d",TrigName[trig].Data(),npt));
				histeff[trig][npt][0]->SetTitle(Form("%s %1.f GeV/c <Pt < %1.f GeV/c %deID",TrigName[trig].Data(),PtEdge[npt],PtEdge[npt+1],eID));
				histeff[trig][npt][0]->GetXaxis()->SetTitle("Cos(#theta)");
				histeff[trig][npt][0]->SetMarkerStyle(24);
				histeff[trig][npt][1]->SetMarkerStyle(20);
				histeff[trig][npt][0]->Draw();
				histeff[trig][npt][1]->Draw("same");
				ratio[trig][npt] = (TH1F*)histeff[trig][npt][0]->Clone();
				ratio[trig][npt]->SetName(Form("ratio%sPt%d",TrigName[trig].Data(),npt));
				ratio[trig][npt]->SetTitle(Form("%s %1.fGeV/c<Pt<%1.fGeV/c %deID",TrigName[trig].Data(),PtEdge[npt],PtEdge[npt+1],eID));
				ratio[trig][npt]->GetYaxis()->SetRangeUser(0.,2.);
				ratio[trig][npt]->GetXaxis()->SetTitle("Cos(#theta)");
				ratio[trig][npt]->SetEntries(histeff[trig][npt][0]->GetEntries());
				ratio[trig][npt]->Divide(histeff[trig][npt][1]);
				canvas[2][trig][eID]->cd(npt+1);
				ratio[trig][npt]->Draw();
			}
			canvas[1][trig][eID]->SaveAs(Form("eff_%s_%deID.pdf",TrigName[trig].Data(),eID));
			canvas[2][trig][eID]->SaveAs(Form("ratio_%s_%deID.pdf",TrigName[trig].Data(),eID));
		}

		gStyle->SetOptFit();
		TLatex text;
		text.SetNDC();
		for(trig=0;trig<4;trig++){
			canvas[4][trig][eID] = new TCanvas(Form("fit_%s_%deID",TrigName[trig].Data(),eID),Form("fit_%s_%deID",TrigName[trig].Data(),eID),1200,800);
			canvas[4][trig][eID]->Divide(3,2);
			for( npt=0;npt<NPT;npt++){
				canvas[4][trig][eID]->cd(npt+1);

				TF1 *fun1 = new TF1("fun1","[0]*(1+[1]*x*x)",-1,1);
				fun1->SetParameters(0.1,0.5);
				fun1->SetLineColor(1);
				TF1 *fun2 = new TF1("fun2","[0]*(1+[1]*x*x)",-1,1);
				fun2->SetParameters(0.1,0.5);
				fun2->SetLineColor(2);
			if((TH1F*)correcteddata1->Get(Form("HistJpsiPhase_%d_%d_0",trig,npt))==0x0)	cout<<"========================"<<(TH1F*)correcteddata1->Get(Form("HistJpsiPhase_%d_%d_0",trig,npt))<<"     trig"<<trig<<"    npt"<<npt<<"    "<<endl;
				if((TH1F*)correcteddata1->Get(Form("HistJpsiPhase_%d_%d_0",trig,npt))!=0x0){
					histfit[trig][npt][0] = (TH1F*)correcteddata1->Get(Form("HistJpsiPhase_%d_%d_0",trig,npt));
					histfit[trig][npt][0]->SetMarkerStyle(24);
					histfit[trig][npt][0]->SetMarkerSize(1.5);
					histfit[trig][npt][0]->SetMarkerColor(1);
					histfit[trig][npt][0]->SetLineColor(1);
					histfit[trig][npt][0]->Draw();
					histfit[trig][npt][0]->Fit("fun1","R");
					histfit[trig][npt][0]->SetTitle(Form("%s %1.f GeV/c <J/#psi Pt <%1.f GeV/c %deID",TrigName[trig].Data(),PtEdge[npt],PtEdge[npt+1],eID));

				}
				if((TH1F*)correcteddata2->Get(Form("HistJpsiPhase_%d_%d_0",trig,npt))!=0x0)	{
					histfit[trig][npt][1] = (TH1F*)correcteddata2->Get(Form("HistJpsiPhase_%d_%d_0",trig,npt));
					histfit[trig][npt][1]->SetMarkerStyle(20);
					histfit[trig][npt][1]->SetMarkerSize(1);
					histfit[trig][npt][1]->SetMarkerColor(2);
					histfit[trig][npt][1]->SetLineColor(2);
					histfit[trig][npt][1]->Draw("same");
					histfit[trig][npt][1]->Fit("fun2","R");
					text.DrawLatex(0.2,0.2,Form("#Delta#lambda = %f",fun1->GetParameter(1)-fun2->GetParameter(1)));
				}
			}
			canvas[4][trig][eID]->SaveAs(Form("fit%seID%d.pdf",TrigName[trig].Data(),eID));
		}
	}
}
