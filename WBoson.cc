#include <iostream>
#include <cmath>
#include "Pythia8/Pythia.h"
#include <fstream>
#include <vector>


using namespace Pythia8;
using namespace std;

int main(){
    Pythia pythia;

    
    //This simulates an n + p collision.
    pythia.readString("Beams:idA = 2112");   // neutron
    pythia.readString("Beams:idB = 2212");   // proton

    // --- 2) Set a center-of-mass energy (keep as you used before; must be > ~80 GeV)
    pythia.readString("Beams:eCM = 200");   // GeV, chosen > m_W

    // --- 3) Enable weak single-boson W production.
    pythia.readString("WeakSingleBoson:ffbar2W = on");

    // --- 4) Restrict W decays to leptonic channels only.
    pythia.readString("24:onMode = off"); 
    // W+ -> e+ nu_e  :  24:onIfMatch = -11 12
    // W- -> e- anti-nu_e : 24:onIfMatch = 11 -12
    pythia.readString("24:onIfMatch = 11 -12");   // W- -> e- anti-nu_e
    pythia.readString("24:onIfMatch = -11 12");   // W+ -> e+ nu_e (kept for completeness)
    // Enable muon channels
    pythia.readString("24:onIfMatch = 13 -14");   // W- -> mu- anti-nu_mu
    pythia.readString("24:onIfMatch = -13 14");   // W+ -> mu+ nu_mu
    // --- 5) (Optional) Turn off QED FSR from final leptons to get a cleaner W mass peak
    pythia.readString("TimeShower:QEDshowerByL = off");

    Hist M("W Resonance Mass Distribution", 300, 50, 110);
    pythia.init();
    const int NBINS = 300;                // must match first arg of Hist
    const double XMIN = 50.0;             // must match lower edge in Hist
    const double XMAX = 110.0;            // must match upper edge in Hist
    const double BINWIDTH = (XMAX - XMIN) / double(NBINS);
    std::vector<double> manualHist(NBINS, 0.0);

    int Nevent = 5000;

    for (int i = 0; i < Nevent; ++i) {
        if (!pythia.next()) continue;
        int entries = pythia.event.size();
        for (int j = 0; j < entries; ++j) {
            int pid = pythia.event[j].id();
            if (pid != -24) continue;
            int d1 = pythia.event[j].daughter1();
            int d2 = pythia.event[j].daughter2();
            if (d1 <= 0 || d2 <= 0) continue;         
            // Use final-state daughters only
            if (!pythia.event[d1].isFinal() || !pythia.event[d2].isFinal()) continue;
            int pdg1 = pythia.event[d1].id();
            int pdg2 = pythia.event[d2].id();
            // For W- we expect (e- (=11) , anti-nu_e (=-12)) or (mu- (=13), anti-nu_mu (=-14))
            // Reconstruct invariant mass from the visible charged lepton and neutrino 4-vectors.
            // Note: neutrinos are present in the event record as particles with momentum,
            // so (p_lepton + p_nu).m() gives the full W invariant mass.
            Vec4 p1 = pythia.event[d1].p();
            Vec4 p2 = pythia.event[d2].p();
            double mrec = (p1 + p2).mCalc();
            M.fill(mrec);
            int ibin = int((mrec - XMIN) / BINWIDTH);
            if (ibin >= 0 && ibin < NBINS) manualHist[ibin] += 1.0;
        }
    }
    cout << M << endl;
    pythia.stat();
    // ---- write numeric bin contents to CSV for reliable external fits ----
{
    std::ofstream fout("WBins.csv");
    fout << "center,content\n";
    for (int i = 0; i < NBINS; ++i) {
        double center = XMIN + (i + 0.5) * BINWIDTH;
        fout << center << "," << manualHist[i] << "\n";
    }
    fout.close();
    std::cout << "W histogram exported to WBins.csv\n";
}

    HistPlot Mw("WBosonMasses");
    Mw.frame("WBosonMasses", "W Resonance Mass Distribution", "Masses (GeV)", "Entries");
    Mw.add(M);
    Mw.plot();

    
    return 0;
}


