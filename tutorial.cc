#include <iostream>

#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main(){	

Pythia pythia;

pythia.readString("Beams:idA=11");
pythia.readString("Beams:idB=-11");
pythia.readString("Beams:eCM=91.2");

//set z-boson production to be on

pythia.readString("WeakSingleBoson:ffbar2gmZ=on");

//set only muon decay for Z bosons

pythia.readString("23:onMode = off");
pythia.readString("23:onIfAny = 13");

pythia.init();

int Nevent=100;

for (int i=0; i<Nevent; i++){

cout<< "Event:" << i <<endl;

if (!pythia.next()){
continue;}

int entries= pythia.event.size();
cout<<"Entries: "<< entries <<endl;

for (int j=0; j<entries; j++){

int id= pythia.event[j].id();

double mass= pythia.event[j].m();

std::string name= pythia.particleData.name(id);

double px= pythia.event[j].px();
double py= pythia.event[j].py();
double pz= pythia.event[j].pz();
double pabs=sqrt(pow(px,2)+pow(py,2)+pow(pz,2));

if (id==13 || id==-13){
cout<<"Muon Detected "<< id <<endl;
cout<<"Mass, Momentum and Energy are "<< mass <<","<< pabs <<","<< pythia.event[j].e() <<endl;
}
else{
cout <<"ID, Mass, Momentum and Name are "<< id <<","<< mass <<"," << pabs <<","<< name <<endl;
}
}

}

}
