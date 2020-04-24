void read_tree(const char* filename)
{
    TFile* file = TFile::Open(filename);

    TTreeReader reader("pfo_tree", file);

    TTreeReaderValue<Int_t> eventId(reader, "eventId");
    TTreeReaderValue<Int_t> pfoId(reader, "pfoId");
    TTreeReaderValue<Int_t> clusterId(reader, "clusterId");
    TTreeReaderValue<Int_t> pfoIsTrack(reader, "pfoIsTrack");
    TTreeReaderValue<Float_t> clusterTrackEnergy(reader, "clusterTrackEnergy");
    TTreeReaderValue<Float_t> clusterEnergy(reader, "clusterEnergy");
    TTreeReaderValue<std::vector<float>> clusterProb(reader, "clusterProb");

    std::map<std::pair<int, int>, float> pfoTrackEnergy;
    std::map<std::pair<int, int>, float> pfoEnergy;
    std::map<std::pair<int, int>, bool> pfoRecoTrack;
    std::map<std::pair<int, int>, std::vector<float>*> pfoProb;

    TCanvas canvas("c1", "c1", 1024, 768);
    canvas.cd();
    canvas.Print("clusters.pdf[");
    while (reader.Next())
    {
        const std::pair key = std::make_pair(*eventId, *pfoId);
        if (pfoProb.find(key) == pfoProb.end())
        {
            std::string pfoHistname = "pfo_hist_E" + std::to_string(*eventId) + "_P" + std::to_string(*pfoId);
            pfoProb.insert(std::make_pair(key, new std::vector<float>()));
            pfoTrackEnergy.insert(std::make_pair(key, 0.f));
            pfoEnergy.insert(std::make_pair(key, 0.f));
            pfoRecoTrack.insert(std::make_pair(key, *pfoIsTrack == 1));
        }
        std::string idStr = "_E" + std::to_string(*eventId) + "_P" + std::to_string(*pfoId) + "_C" + std::to_string(*clusterId);
        std::string clusterHistname = "cluster_hist" + idStr;
        TH1F* hist = new TH1F(clusterHistname.c_str(), "P(Track) Distribution", 25, 0.f, 1.f);
        for (int i = 0; i < clusterProb->size(); ++i)
        {
            hist->Fill(clusterProb->at(i));
            pfoProb.at(key)->push_back(clusterProb->at(i));
            pfoTrackEnergy.at(key) += *clusterTrackEnergy;
            pfoEnergy.at(key) += *clusterEnergy;
        }
        float mean = hist->GetMean();
        float rms = hist->GetRMS();
        delete hist;

        std::sort(clusterProb->begin(), clusterProb->end());
        std::vector<float> sample(clusterProb->begin(), clusterProb->end());
        const int N = static_cast<int>(std::floor(0.1 * clusterProb->size()));
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

        std::string recoStr = "Reco Track = " + std::to_string(pfoRecoTrack.at(key));
        TLatex recoTxt(0, 0, recoStr.c_str());
        recoTxt.SetTextFont(font); recoTxt.SetTextSize(fontSize);
        listOfLines->Add(&recoTxt);

        bool isMC = *clusterTrackEnergy >= (0.5f * *clusterEnergy) ? true : false;
        std::string mcStr = "MC Track = " + std::to_string(isMC);
        TLatex mcTxt(0, 0, mcStr.c_str());
        mcTxt.SetTextFont(font); mcTxt.SetTextSize(fontSize);
        listOfLines->Add(&mcTxt);

        // Stop the stats being redrawn
        hist->SetStats(0);
        canvas.Modified();

        canvas.Print("clusters.pdf");

        delete hist;
    }
    canvas.Print("clusters.pdf]");

    TH1F* mcTracks = new TH1F("mcTracks", "P(Track) Distribution", 25, 0.f, 1.f);
    mcTracks->SetLineWidth(2);
    TH1F* mcShowers = new TH1F("mcShowers", "P(Track) Distribution", 25, 0.f, 1.f);
    mcShowers->SetLineColor(kRed);
    mcShowers->SetLineWidth(2);

    int muRmsPerf[2][3]{0}, muRms90Perf[2][3]{0}, tailPerf[2][3]{0};
    canvas.cd();
    canvas.Print("pfos.pdf[");
    for (const auto [key, probs] : pfoProb)
    {
        std::string pfoHistname = "pfo_hist_E" + std::to_string(key.first) + "_P" + std::to_string(key.second);
        TH1F* hist = new TH1F(pfoHistname.c_str(), "P(Track) Distribution", 25, 0.f, 1.f);
        for (float p : *probs)
            hist->Fill(p);
        float mean = hist->GetMean();
        float rms = hist->GetRMS();
        delete hist;

        std::sort(probs->begin(), probs->end());
        // Calculated left and right fractions
        int nShower{0}, nTrack{0};
        for (float p : *probs)
        {
            if (p < 0.5)
                nShower++;
            else
                nTrack++;
        }
        float fShower{static_cast<float>(nShower) / (nShower + nTrack)};
        float fTrack{1.f - fShower};

        // Cut sample size to 90% based on minimum RMS
        std::vector<float> sample(probs->begin(), probs->end());
        const int N = static_cast<int>(std::floor(0.1 * probs->size()));
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
        hist = new TH1F(pfoHistname.c_str(), "P(Track) Distribution", 25, 0.f, 1.f);
        for (float p : sample)
            hist->Fill(p);
        float mean90 = hist->GetMean();
        float rms90 = hist->GetRMS();

        bool isMC = pfoTrackEnergy.at(key) >= (0.5f * pfoEnergy.at(key)) ? true : false;
        if (isMC)
        {
            for (int i = 0; i < hist->GetNbinsX(); ++i)
                mcTracks->SetBinContent(i, mcTracks->GetBinContent(i) + hist->GetBinContent(i));
        }
        else
        {
            for (int i = 0; i < hist->GetNbinsX(); ++i)
                mcShowers->SetBinContent(i, mcShowers->GetBinContent(i) + hist->GetBinContent(i));
        }

        //hist->Scale(1.f / hist->Integral(), "nosw2 width");
        hist->Scale(1.f / hist->Integral(), "nosw2");
        hist->Draw();
        canvas.Update();

        TPaveStats* stats = dynamic_cast<TPaveStats*>(canvas.GetPrimitive("stats"));
        stats->SetName("pfostats");
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

        std::string sfracStr = "Sfrac = " + std::to_string(fShower);
        TLatex sfracTxt(0, 0, sfracStr.c_str());
        sfracTxt.SetTextFont(font); sfracTxt.SetTextSize(fontSize);
        listOfLines->Add(&sfracTxt);

        std::string tfracStr = "Tfrac = " + std::to_string(fTrack);
        TLatex tfracTxt(0, 0, tfracStr.c_str());
        tfracTxt.SetTextFont(font); tfracTxt.SetTextSize(fontSize);
        listOfLines->Add(&tfracTxt);

        std::string recoStr = "Reco Track = " + std::to_string(pfoRecoTrack.at(key));
        TLatex recoTxt(0, 0, recoStr.c_str());
        recoTxt.SetTextFont(font); recoTxt.SetTextSize(fontSize);
        listOfLines->Add(&recoTxt);

        std::string mcStr = "MC Track = " + std::to_string(isMC);
        TLatex mcTxt(0, 0, mcStr.c_str());
        mcTxt.SetTextFont(font); mcTxt.SetTextSize(fontSize);
        listOfLines->Add(&mcTxt);

        // Stop the stats being redrawn
        hist->SetStats(0);
        canvas.Modified();

        canvas.Print("pfos.pdf");

        int row{pfoRecoTrack.at(key) == isMC ? 0 : 1};  // Reco correct : Reco wrong
        int state{mean >= 0.5 ?
            (mean - rms > 0.5 ? 1 : 2) :        // Track  : Ambiguous
            (mean + rms < 0.5 ? 0 : 2) };       // Shower : Ambiguous
        int col{state != 2 && isMC == state ?
            0 :                                 // Network correct
            state != 2 ? 1 : 2};                // Network wrong : Ambiguous
        muRmsPerf[row][col]++;

        state = mean90 >= 0.5 ?
            (mean90 - rms90 > 0.5 ? 1 : 2) :    // Track  : Ambiguous
            (mean90 + rms90 < 0.5 ? 0 : 2);     // Shower : Ambiguous
        col = state != 2 && isMC == state ?
            0 :                                 // Network correct
            state != 2 ? 1 : 2;                 // Network wrong : Ambiguous
        muRms90Perf[row][col]++;

        float threshold = 2.f / 3.f;
        state = fTrack >= threshold ?
            1 :                                 // Track
            fTrack <= 1.f - threshold ? 0 : 2;  // Shower : Ambiguous
        col = state != 2 && isMC == state ?
            0 :                                 // Network correct
            state != 2 ? 1 : 2;                 // Network wrong : Ambiguous
        tailPerf[row][col]++;

        delete hist;
    }
    canvas.Print("pfos.pdf]");

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

    gStyle->SetOptStat(0);
    canvas.cd();
    mcTracks->Scale(1.f / mcTracks->Integral(), "nosw2");
    mcShowers->Scale(1.f / mcShowers->Integral(), "nosw2");
    mcTracks->SetMaximum(std::max(mcTracks->GetMaximum(), mcShowers->GetMaximum()));
    mcShowers->SetMaximum(std::max(mcTracks->GetMaximum(), mcShowers->GetMaximum()));
    mcTracks->Draw();
    mcShowers->Draw("same");

    canvas.Print("global.pdf");

    for (const auto [key, probs] : pfoProb)
        delete probs;

    delete mcTracks;
    delete mcShowers;

    file->Close();
}

