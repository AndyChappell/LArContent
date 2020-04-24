void read_avail_clusters(const char* filename)
{
    TFile* file = TFile::Open(filename);

    TTreeReader reader("cluster_tree", file);

    TTreeReaderValue<Int_t> eventId(reader, "eventId");
    TTreeReaderValue<Int_t> viewId(reader, "viewId");
    TTreeReaderValue<Int_t> clusterId(reader, "clusterId");
    TTreeReaderValue<Int_t> isTrack(reader, "isTrack");
    TTreeReaderValue<std::vector<float>> clusterProb(reader, "clusterProb");

    const int Nbins{18};
    const float edges[Nbins + 1]{1.f,2.f,3.f,5.f,10.f,15.f,20.f,25.f,30.f,35.f,40.f,50.f,75.f,100.f,150.f,200.f,250.f,500.f,3000.f};
    TH1F nhitsHist("nhits", "NHits Distribution", Nbins, edges);
    TH1F recoCorrHist("nreco", "Correct fraction", Nbins, edges);
    recoCorrHist.SetXTitle("N hits");
    recoCorrHist.SetYTitle("Correct fraction");
    TH1F netCorrHist("nnet", "Correct fraction", Nbins, edges);
    TH1F totHist("tot", "N", Nbins, edges);
    int muRmsPerf[2][3]{0}, muRms90Perf[2][3]{0}, tailPerf[2][3]{0};
    std::cout << "Clusters: " << reader.GetEntries() << std::endl;
    TCanvas canvas("c1", "c1", 1024, 768);
    canvas.cd();
    const int low{100}, high{1000};
    canvas.Print(std::string("avail_clusters_" + std::to_string(high) + ".pdf[").c_str());
    while (reader.Next())
    {
        int nhits = clusterProb->size();
        nhitsHist.Fill(nhits);
        //if (nhits <= low || nhits > high)
        //    continue;
        std::string idStr = "_E" + std::to_string(*eventId) + "_V" + std::to_string(*viewId) + "_C" + std::to_string(*clusterId);
        std::string clusterHistname = "cluster_hist" + idStr;
        TH1F* hist = new TH1F(clusterHistname.c_str(), "P(Track) Distribution", 25, 0.f, 1.f);
        for (int i = 0; i < nhits; ++i)
        {
            hist->Fill(clusterProb->at(i));
        }
        float mean = hist->GetMean();
        float rms = hist->GetRMS();
        delete hist;

        std::sort(clusterProb->begin(), clusterProb->end());
        // Calculated left and right fractions
        int nShower{0}, nTrack{0};
        for (float p : *clusterProb)
        {
            if (p < 0.5)
                nShower++;
            else
                nTrack++;
        }
        float fTrack{static_cast<float>(nTrack) / (nShower + nTrack)};

        std::vector<float> sample(clusterProb->begin(), clusterProb->end());
        const int N = static_cast<int>(std::floor(0.1 * nhits));
        for (int i = 0; i < N; ++i)
        {
            std::vector<float> lsample(sample.begin(), sample.end() - 1);
            TH1F lhist("lhist", "P(Track) Distribution", 25, 0.f, 1.f);
            for (float p : lsample)
                lhist.Fill(p);
            std::vector<float> rsample(sample.begin() + 1, sample.end());
            TH1F rhist("rhist", "P(Track) Distribution", 25, 0.f, 1.f);
            for (float p : rsample)
                rhist.Fill(p);
            float lrms = lhist.GetRMS();
            float rrms = lhist.GetRMS();
            if (lrms < rrms)
                sample = lsample;
            else
                sample = rsample;
        }
        // Build 90% sample of hits with lowest RMS
        hist = new TH1F(clusterHistname.c_str(), "P(Track) Distribution", 25, 0.f, 1.f);
        for (float p : sample)
            hist->Fill(p);
        float mean90 = hist->GetMean();
        float rms90 = hist->GetRMS();

        hist->SetDirectory(0);
        hist->Scale(1.f / hist->Integral(), "nosw2");
        hist->Draw();
        canvas.Update();

        TPaveStats* stats = dynamic_cast<TPaveStats*>(canvas.GetPrimitive("stats"));
        stats->SetName("clusterstats");
        TList* listOfLines = stats->GetListOfLines();

        Font_t font = stats->GetLineWith("Mean")->GetTextFont();
        Float_t fontSize = stats->GetLineWith("Mean")->GetTextSize();
        // Remove lines
        listOfLines->Remove(stats->GetLineWith("Mean"));
        listOfLines->Remove(stats->GetLineWith("Std Dev"));

        // Add lines
        std::string meanStr = "Mean = " + std::to_string(mean);
        TLatex meanTxt(0, 0, meanStr.c_str());
        meanTxt.SetTextFont(font); meanTxt.SetTextSize(fontSize);
        listOfLines->Add(&meanTxt);

        std::string rmsStr = "RMS = " + std::to_string(rms);
        TLatex rmsTxt(0, 0, rmsStr.c_str());
        rmsTxt.SetTextFont(font); rmsTxt.SetTextSize(fontSize);
        listOfLines->Add(&rmsTxt);

        std::string mean90Str = "Mean90 = " + std::to_string(mean90);
        TLatex mean90Txt(0, 0, mean90Str.c_str());
        mean90Txt.SetTextFont(font); mean90Txt.SetTextSize(fontSize);
        listOfLines->Add(&mean90Txt);

        std::string rms90Str = "RMS90 = " + std::to_string(rms90);
        TLatex rms90Txt(0, 0, rms90Str.c_str());
        rms90Txt.SetTextFont(font); rms90Txt.SetTextSize(fontSize);
        listOfLines->Add(&rms90Txt);

        std::string tfracStr = "Tfrac = " + std::to_string(fTrack);
        TLatex tfracTxt(0, 0, tfracStr.c_str());
        tfracTxt.SetTextFont(font); tfracTxt.SetTextSize(fontSize);
        listOfLines->Add(&tfracTxt);

        std::string mcStr = "MC Track = " + std::to_string(*isTrack);
        TLatex mcTxt(0, 0, mcStr.c_str());
        mcTxt.SetTextFont(font); mcTxt.SetTextSize(fontSize);
        listOfLines->Add(&mcTxt);

        // Stop the stats being redrawn
        hist->SetStats(0);
        canvas.Modified();

        //canvas.Print(std::string("avail_clusters_" + std::to_string(high) + ".pdf").c_str());

        int row{*isTrack ? 1 : 0};  // Reco correct : Reco wrong
        int state{mean >= 0.5 ?
            (mean - rms > 0.5 ? 1 : 2) :        // Track  : Ambiguous
            (mean + rms < 0.5 ? 0 : 2) };       // Shower : Ambiguous
        int col{state != 2 && *isTrack == state ?
            0 :                                 // Network correct
            state != 2 ? 1 : 2};                // Network wrong : Ambiguous
        muRmsPerf[row][col]++;
        totHist.Fill(nhits);
        if (row == 0)
        {
            recoCorrHist.Fill(nhits);
            if (col != 1)
                netCorrHist.Fill(nhits);
        }
        else
        {
            if (col == 0)
                netCorrHist.Fill(nhits);
        }

        state = mean90 >= 0.5 ?
            (mean90 - rms90 > 0.5 ? 1 : 2) :    // Track  : Ambiguous
            (mean90 + rms90 < 0.5 ? 0 : 2);     // Shower : Ambiguous
        col = state != 2 && *isTrack == state ?
            0 :                                 // Network correct
            state != 2 ? 1 : 2;                 // Network wrong : Ambiguous
        muRms90Perf[row][col]++;

        float threshold = 2.f / 3.f;
        state = fTrack >= threshold ?
            1 :                                 // Track
            fTrack <= 1.f - threshold ? 0 : 2;  // Shower : Ambiguous
        col = state != 2 && *isTrack == state ?
            0 :                                 // Network correct
            state != 2 ? 1 : 2;                 // Network wrong : Ambiguous
        tailPerf[row][col]++;


        delete hist;
    }
    canvas.Print(std::string("avail_clusters_" + std::to_string(high) + ".pdf]").c_str());

    gPad->SetLogx(); gPad->SetLogy();
    nhitsHist.SetMinimum(1.f);
    nhitsHist.Draw();
    std::cout << "Underflow: " << nhitsHist.GetBinContent(0) << " Overflow: " << nhitsHist.GetBinContent(Nbins+1) << std::endl;
    canvas.Print("nhits.pdf");
    gPad->SetLogx(0); gPad->SetLogy(0);

    gStyle->SetOptStat(0);
    gPad->SetLogx(1);
    recoCorrHist.Divide(&totHist); recoCorrHist.SetMinimum(0.f); recoCorrHist.SetMaximum(1.f); recoCorrHist.SetLineColor(kBlack);
    netCorrHist.Divide(&totHist); netCorrHist.SetMinimum(0.f); netCorrHist.SetMaximum(1.f); netCorrHist.SetLineColor(kGreen + 1);
    recoCorrHist.Draw();
    netCorrHist.Draw("same");

    TLegend legend(0.6, 0.8, 0.9, 0.9);
    legend.AddEntry(&recoCorrHist, "Standard Reco", "l");
    legend.AddEntry(&netCorrHist, "Network", "l");
    legend.Draw();

    canvas.Print("correct_fraction.pdf");
    gPad->SetLogx(0);

    std::cout << "RMS:" << std::endl;
    for (int r = 0; r <= 1; ++r)
    {
        for (int c = 0; c <= 2; ++c)
        {
            std::cout << muRmsPerf[r][c] << "   ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nRMS90:" << std::endl;
    for (int r = 0; r <= 1; ++r)
    {
        for (int c = 0; c <= 2; ++c)
        {
            std::cout << muRms90Perf[r][c] << "   ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nTails:" << std::endl;
    for (int r = 0; r <= 1; ++r)
    {
        for (int c = 0; c <= 2; ++c)
        {
            std::cout << tailPerf[r][c] << "   ";
        }
        std::cout << std::endl;
    }

    file->Close();
}

