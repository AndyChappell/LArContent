void read_tree()
{
    TFile* file = TFile::Open("net.root");

    TTreeReader reader("net_tree", file);

    TTreeReaderValue<Int_t> eventId(reader, "eventId");
    TTreeReaderValue<Int_t> pfoId(reader, "pfoId");
    TTreeReaderValue<Int_t> clusterId(reader, "clusterId");
    TTreeReaderValue<Int_t> pfoIsTrack(reader, "pfoIsTrack");
    TTreeReaderValue<Float_t> mcTrackFraction(reader, "clusterMcTrackFraction");
    TTreeReaderValue<std::vector<float>> clusterProb(reader, "clusterProb");

    std::map<int, TH1F*> pfoHists;

    while (reader.Next())
    {
        if (pfoHists.find(*pfoId) == pfoHists.end())
        {
            std::string pfoHistname = "pfo_hist" + std::to_string(*eventId) + "." + std::to_string(*pfoId);
            pfoHists.insert(std::make_pair(*pfoId, new TH1F(pfoHistname.c_str(), "Track probabilities", 100, 0.f, 1.f)));
            pfoHists.at(*pfoId)->SetDirectory(0);
        }
        std::string idStr = std::to_string(*eventId) + "." + std::to_string(*pfoId) + "." + std::to_string(*clusterId);
        std::string clusterHistname = "cluster_hist" + idStr;
        TH1F frequencies(clusterHistname.c_str(), "Track probailities", 100, 0.f, 1.f);
        for (int i = 0; i < clusterProb->size(); ++i)
        {
            frequencies.Fill(clusterProb->at(i));
            pfoHists.at(*pfoId)->Fill(clusterProb->at(i));
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
        hist->Draw();
        std::string filename = std::string(hist->GetName()) + ".pdf";
        canvas.SaveAs(filename.c_str());
    }

    file->Close();
}

