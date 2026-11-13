// Scattering angle vs impact parameter with the analytic Rutherford curve.
//
// Usage (batch):
//   root -l -b -q 'macro/plot_rutherford.C("Data/test_0.root")'
//   root -l -b -q 'macro/plot_rutherford.C("Data/test_0.root", "TT_0", "rutherford.png", 0.05, 4, 2, 79)'
//
// Arguments
//   file : ROOT file with the tree (branches "imp" [fm] and "SAngle" [rad])
//   tree : tree name (default TT_0)
//   out  : output image (png/pdf/...), nothing is saved when empty
//   v, m, q1, q2 : projectile velocity [c], mass [amu], charges [e] (defaults: alpha on Au)
//
// Analytic: theta = 2 atan( k q1 q2 / (m v^2 b) ),  k = 0.00154629 (fm, amu, c, e units).
// Residual pad: simulation - analytic. Note that "imp" and "v" are the values at the
// finite starting distance, not the asymptotic ones, so a few-percent systematic
// residual is expected.
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TStyle.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <vector>

void plot_rutherford(const char * file, const char * tree = "TT_0", const char * out = "",
                     double v = 0.05, double m = 4, double q1 = 2, double q2 = 79)
{
    const double k = 0.00154629;

    TFile * f = TFile::Open(file);
    if (!f || f->IsZombie()) { Error("plot_rutherford", "cannot open %s", file); return; }
    TTree * t = (TTree*) f->Get(tree);
    if (!t) { Error("plot_rutherford", "no tree %s in %s", tree, file); return; }
    if (!t->GetBranch("imp") || !t->GetBranch("SAngle")) {
        Error("plot_rutherford", "tree %s needs branches imp and SAngle", tree); return;
    }
    Long64_t n = t->GetEntries();
    if (n < 1) { Error("plot_rutherford", "tree %s is empty", tree); return; }

    double imp, ang;
    t->SetBranchAddress("imp", &imp);
    t->SetBranchAddress("SAngle", &ang);
    std::vector<double> b(n), a(n);
    for (Long64_t i = 0; i < n; i++) { t->GetEntry(i); b[i] = imp; a[i] = ang; }

    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int x, int y) { return b[x] < b[y]; });

    auto theory = [&](double bb) { return 2 * TMath::ATan(k * q1 * q2 / (m * v * v * bb)); };

    std::vector<double> xs(n), ys(n), res(n);
    double maxRes = 0, sum2 = 0;
    for (Long64_t j = 0; j < n; j++) {
        int i = idx[j];
        xs[j] = b[i]; ys[j] = a[i];
        res[j] = a[i] - theory(b[i]);
        if (std::fabs(res[j]) > std::fabs(maxRes)) maxRes = res[j];
        sum2 += res[j] * res[j];
    }
    double rms = std::sqrt(sum2 / n);
    printf("plot_rutherford: %lld events, v=%g m=%g q1=%g q2=%g\n", n, v, m, q1, q2);
    printf("  max |residual| = %.5g rad (signed %.5g), RMS residual = %.5g rad\n",
           std::fabs(maxRes), maxRes, rms);

    // smooth analytic curve
    const int nc = 400;
    double bmin = xs.front(), bmax = xs.back();
    std::vector<double> cx(nc), cy(nc);
    for (int i = 0; i < nc; i++) { cx[i] = bmin + (bmax - bmin) * i / (nc - 1); cy[i] = theory(cx[i]); }

    gStyle->SetOptStat(0);
    TCanvas * c = new TCanvas("c_rutherford", "Rutherford scattering", 700, 650);
    TPad * p1 = new TPad("p1", "p1", 0.0, 0.3, 1.0, 1.0);
    TPad * p2 = new TPad("p2", "p2", 0.0, 0.0, 1.0, 0.3);
    p1->SetBottomMargin(0.02); p1->SetGrid();
    p2->SetTopMargin(0.03); p2->SetBottomMargin(0.35); p2->SetGrid();
    p1->Draw(); p2->Draw();

    p1->cd();
    TGraph * gs = new TGraph(n, xs.data(), ys.data());
    gs->SetTitle(";;scattering angle [rad]");
    gs->SetMarkerStyle(20); gs->SetMarkerSize(0.7); gs->SetMarkerColor(kBlue + 1);
    gs->GetXaxis()->SetLimits(bmin - 0.02 * (bmax - bmin), bmax + 0.02 * (bmax - bmin));
    gs->GetXaxis()->SetLabelSize(0);
    gs->GetYaxis()->SetTitleSize(0.05); gs->GetYaxis()->SetTitleOffset(0.9);
    gs->Draw("AP");
    TGraph * gt = new TGraph(nc, cx.data(), cy.data());
    gt->SetLineColor(kRed + 1); gt->SetLineWidth(2);
    gt->Draw("L SAME");
    TLegend * lg = new TLegend(0.55, 0.72, 0.89, 0.88);
    lg->AddEntry(gs, "simulation", "p");
    lg->AddEntry(gt, "#theta = 2 atan(k q_{1}q_{2}/(m v^{2} b))", "l");
    lg->Draw();

    p2->cd();
    TGraph * gr = new TGraph(n, xs.data(), res.data());
    gr->SetTitle(";impact parameter b [fm];sim - analytic [rad]");
    gr->SetMarkerStyle(20); gr->SetMarkerSize(0.6); gr->SetMarkerColor(kBlue + 1);
    gr->GetXaxis()->SetLimits(bmin - 0.02 * (bmax - bmin), bmax + 0.02 * (bmax - bmin));
    gr->GetXaxis()->SetTitleSize(0.12); gr->GetXaxis()->SetLabelSize(0.10); gr->GetXaxis()->SetTitleOffset(1.0);
    gr->GetYaxis()->SetTitleSize(0.10); gr->GetYaxis()->SetLabelSize(0.09); gr->GetYaxis()->SetTitleOffset(0.45);
    gr->GetYaxis()->SetNdivisions(505);
    gr->Draw("AP");
    TLine * l0 = new TLine(bmin - 0.02 * (bmax - bmin), 0, bmax + 0.02 * (bmax - bmin), 0);
    l0->SetLineStyle(2); l0->Draw();

    c->cd();
    c->Update();
    if (out && out[0]) {
        c->SaveAs(out);
        printf("  saved %s\n", out);
    }
}
