void read_tree()
{
    //gStyle->SetOptStat(0);
    TFile* file = TFile::Open("net.root");

    TTreeReader reader("net_tree", file);

    TTreeReaderValue<Int_t> eventId(reader, "eventId");
    TTreeReaderValue<Int_t> pfoId(reader, "pfoId");
    TTreeReaderValue<Int_t> clusterId(reader, "clusterId");
    TTreeReaderValue<Int_t> pfoIsTrack(reader, "pfoIsTrack");
    TTreeReaderValue<Float_t> clusterTrackEnergy(reader, "clusterTrackEnergy");
    TTreeReaderValue<Float_t> clusterEnergy(reader, "clusterEnergy");
    TTreeReaderValue<std::vector<float>> clusterProb(reader, "clusterProb");

    std::map<int, TH1F*> pfoHists;
    std::map<int, float> pfoTrackEnergy;
    std::map<int, float> pfoEnergy;
    std::map<int, bool> pfoRecoTrack;

    while (reader.Next())
    {
        if (pfoHists.find(*pfoId) == pfoHists.end())
        {
            std::string pfoHistname = "pfo_hist_E" + std::to_string(*eventId) + "_P" + std::to_string(*pfoId);
            pfoHists.insert(std::make_pair(*pfoId, new TH1F(pfoHistname.c_str(), "P(Track) Distribution", 25, 0.f, 1.f)));
            pfoHists.at(*pfoId)->SetDirectory(0);
            pfoTrackEnergy.insert(std::make_pair(*pfoId, 0.f));
            pfoEnergy.insert(std::make_pair(*pfoId, 0.f));
            pfoRecoTrack.insert(std::make_pair(*pfoId, *pfoIsTrack == 1));
        }
        std::string idStr = "_E" + std::to_string(*eventId) + "_P" + std::to_string(*pfoId) + "_C" + std::to_string(*clusterId);
        std::string clusterHistname = "cluster_hist" + idStr;
        TH1F frequencies(clusterHistname.c_str(), "Track probailities", 25, 0.f, 1.f);
        for (int i = 0; i < clusterProb->size(); ++i)
        {
            frequencies.Fill(clusterProb->at(i));
            pfoHists.at(*pfoId)->Fill(clusterProb->at(i));
            pfoTrackEnergy.at(*pfoId) += *clusterTrackEnergy;
            pfoEnergy.at(*pfoId) += *clusterEnergy;
        }
        frequencies.SetDirectory(0);
        TCanvas canvas("c1", "c1", 1024, 768);
        canvas.cd();
        frequencies.Draw();
        std::string filename = std::string(frequencies.GetName()) + ".pdf";
        canvas.SaveAs(filename.c_str());
    }

    for (const auto [key, hist] : pfoHists)
    {
        TCanvas canvas("c1", "c1", 1024, 768);
        canvas.cd();
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
        listOfLines->Remove(stats->GetLineWith("Std Dev"));

        // Add lines
        std::string rmsStr = "RMS = " + std::to_string(hist->GetRMS());
        TLatex rmsTxt(0, 0, rmsStr.c_str());
        rmsTxt.SetTextFont(font); rmsTxt.SetTextSize(fontSize);
        listOfLines->Add(&rmsTxt);

        std::string recoStr = "Reco Track = " + std::to_string(pfoRecoTrack.at(key));
        TLatex recoTxt(0, 0, recoStr.c_str());
        recoTxt.SetTextFont(font); recoTxt.SetTextSize(fontSize);
        listOfLines->Add(&recoTxt);

        bool isMC = pfoTrackEnergy.at(key) >= (0.5f * pfoEnergy.at(key)) ? true : false;
        std::string mcStr = "MC Track = " + std::to_string(isMC);
        TLatex mcTxt(0, 0, mcStr.c_str());
        mcTxt.SetTextFont(font); mcTxt.SetTextSize(fontSize);
        listOfLines->Add(&mcTxt);

        // Stop the stats being redrawn
        hist->SetStats(0);
        canvas.Modified();

        std::string filename = std::string(hist->GetName()) + ".pdf";
        canvas.SaveAs(filename.c_str());
    }

    file->Close();
}

