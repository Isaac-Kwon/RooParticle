// Compare the "even" and "odd" integration results of the Rutherford inspect test.
// Usage: root -l 'even_odd.cpp'
//        root -l 'even_odd.cpp("even.root", "odd.root", "TT_0")'
// Top pad: imp vs SAngle for both files, bottom pad: imp vs (SAngle_even - SAngle_odd).
void even_odd(const char * fileEven = "./build/Data/test_rutherford_inspect_even_0.root",
              const char * fileOdd  = "./build/Data/test_rutherford_inspect_odd_0.root",
              const char * treeName = "TT_0")
{
    TFile * file1_ = TFile::Open(fileEven);
    TFile * file2_ = TFile::Open(fileOdd);
    if (!file1_ || file1_->IsZombie()) { Error("even_odd", "cannot open %s", fileEven); return; }
    if (!file2_ || file2_->IsZombie()) { Error("even_odd", "cannot open %s", fileOdd);  return; }

    TTree * tree1_ = (TTree*) file1_->Get(treeName);
    TTree * tree2_ = (TTree*) file2_->Get(treeName);
    if (!tree1_) { Error("even_odd", "no tree %s in %s", treeName, fileEven); return; }
    if (!tree2_) { Error("even_odd", "no tree %s in %s", treeName, fileOdd);  return; }
    tree1_->SetName("tree1_");
    tree2_->SetName("tree2_");

    tree1_->AddFriend("tree2_");

    TCanvas * c1 = new TCanvas("c1", "c1", 600, 400);
    TPad    * p1 = new TPad("p1", "p1", 0.0, 0.3, 1.0, 1.0);
    TPad    * p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 0.3);

    tree1_->Draw("imp:SAngle", "", "goff");
    TGraph  * g1 = new TGraph(tree1_->GetSelectedRows(), tree1_->GetV1(),  tree1_->GetV2() );
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.8);
    g1->SetName("g1");

    tree2_->Draw("imp:SAngle", "", "goff");
    TGraph  * g2 = new TGraph(tree2_->GetSelectedRows(), tree2_->GetV1(),  tree2_->GetV2() );
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.8);
    g2->SetName("g2");

    tree1_->Draw("imp:(SAngle-tree2_.SAngle)", "", "goff");
    TGraph  * g3 = new TGraph(tree1_->GetSelectedRows(), tree1_->GetV1(),  tree1_->GetV2() );
    g3->SetMarkerStyle(20);
    g3->SetMarkerSize(0.8);

    TMultiGraph * mg = new TMultiGraph();
    mg->Add(g1);
    mg->Add(g2);
    p1->Draw();
    p1->cd();
    mg->Draw("AP");
    c1->cd();
    p2->Draw();
    p2->cd();
    g3->Draw("AP");
}
