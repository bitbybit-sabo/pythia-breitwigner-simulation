#include <iostream>
#include <cmath>
#include "Pythia8/Pythia.h"
#include <vector>
#include <fstream>

using namespace Pythia8;
using namespace std;

int main(){ 
    Pythia pythia;

    pythia.readString("Beams:idA=11");
    pythia.readString("Beams:idB=-11");
    pythia.readString("Beams:eCM=200");
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on"); //Z Boson-Photon continuum
    
    pythia.readString("23:onMode = off"); 
    pythia.readString("23:onIfMatch = 11 -11");  // allow e+ e-
    pythia.readString("23:onIfMatch = 13 -13");  // allow mu+ mu-
    pythia.readString("TimeShower:QEDshowerByL = off");


    Hist M("Z Resonance Mass Distribution", 400, 80, 110); //larger cut from 60 to 140GeV
    const int NBINS = 400;
    const double XMIN = 60.0;   // left edge (GeV)
    const double XMAX = 140.0;  // right edge (GeV)
    const double BINWIDTH = (XMAX - XMIN) / double(NBINS);
    std::vector<double> manualHist(NBINS, 0.0);

    pythia.init();

    int Nevent=5000;

    for (int i=0; i<Nevent; i++){
        cout<< "Event:" << i <<endl;

        if (!pythia.next()) continue;

        int entries= pythia.event.size();
        cout<<"Entries: "<< entries <<endl;

        for (int j=0; j<entries; j++){
            int id= pythia.event[j].id();
            if (id==23){
                if (!pythia.event[j].isFinal() && pythia.event[j].status() >= 0) continue;  // pick decaying Z only
                cout<<"Z Boson Detected "<<endl;
                int d1= pythia.event[j].daughter1();
                int d2= pythia.event[j].daughter2();
                if (!pythia.event[d1].isFinal() || !pythia.event[d2].isFinal()) continue; //pick final daughters only
                Vec4 p1= pythia.event[d1].p();
                Vec4 p2= pythia.event[d2].p();

                double mrec=(p1+p2).mCalc();
                M.fill(mrec);
                int ibin = int((mrec - XMIN) / BINWIDTH);
                if (ibin >= 0 && ibin < NBINS) manualHist[ibin] += 1.0;
            }
        }
    }
    cout<<M<<endl;
    pythia.stat();
{
    std::ofstream fout("ZBins.csv");
    fout << "center,content\n";
    for (int i = 0; i < NBINS; ++i) {
        double center = XMIN + (i + 0.5) * BINWIDTH;
        fout << center << "," << manualHist[i] << "\n";
    }
    fout.close();
    std::cout << "Z histogram exported to ZBins.csv\n"<<endl;
}

    HistPlot Mz("ZBosonMasses");
    Mz.frame("ZBosonMasses","Z Resonance Mass Distribution","Masses","Entries");
    Mz.add(M);
    Mz.plot();
    
    return 0;
}

