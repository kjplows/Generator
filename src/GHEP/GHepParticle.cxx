//____________________________________________________________________________
/*!

\class   genie::GHepParticle

\brief   STDHEP-like event record entry that can fit a particle or a nucleus.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created November 18, 2004

*/
//____________________________________________________________________________

#include <cassert>
#include <iomanip>

#include <TMath.h>

#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"

using namespace genie;

using std::setw;
using std::setprecision;
using std::setfill;
using std::ios;
using std::cout;

const double kPCutOff    = 1e-15;
const double kOffShellDm = 0.002; // 2 MeV

ClassImp(GHepParticle)

//____________________________________________________________________________
namespace genie {
 ostream & operator<< (ostream& stream, const GHepParticle & particle)
 {
   particle.Print(stream);

   return stream;
 }
}
//___________________________________________________________________________
GHepParticle::GHepParticle()
{
  this->Init();
}
//___________________________________________________________________________
// TParticle-like constructor
GHepParticle::GHepParticle(int pdg, GHepStatus_t status,
            int mother1, int mother2, int daughter1, int daughter2,
                        const TLorentzVector & p, const TLorentzVector & v) :
fStatus(status),
fFirstMother(mother1),
fLastMother(mother2),
fFirstDaughter(daughter1),
fLastDaughter(daughter2)
{
  this->SetPdgCode(pdg);

  fP4 = new TLorentzVector(p);
  fX4 = new TLorentzVector(v);
}
//___________________________________________________________________________
// TParticle-like constructor
GHepParticle::GHepParticle(int pdg, GHepStatus_t status,
         int mother1, int mother2, int daughter1, int daughter2,
                           double px, double py, double pz, double E,
                                    double x, double y, double z, double t) :
fStatus(status),
fFirstMother(mother1),
fLastMother(mother2),
fFirstDaughter(daughter1),
fLastDaughter(daughter2)
{
  this->SetPdgCode(pdg);

  fP4 = new TLorentzVector(px,py,pz,E);
  fX4 = new TLorentzVector(x,y,z,t);
}
//___________________________________________________________________________
// Copy constructor
GHepParticle::GHepParticle(const GHepParticle & particle)
{
  this->Init();
  this->Copy(particle);
}
//___________________________________________________________________________
GHepParticle::~GHepParticle()
{
  this->CleanUp();
}
//___________________________________________________________________________
string GHepParticle::Name(void) const
{
  this->AssertIsKnownParticle();

  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);
  return p->GetName();
}
//___________________________________________________________________________
double GHepParticle::Mass(void) const
{
  this->AssertIsKnownParticle();

  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);
  return p->Mass();
}
//___________________________________________________________________________
double GHepParticle::Charge(void) const
{
  this->AssertIsKnownParticle();

  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);
  return p->Charge();
}
//___________________________________________________________________________
double GHepParticle::KinE(bool mass_from_pdg) const
{
  if(!fP4) {
    LOG("GHepParticle", pWARN) << "4-momentum not yet set!";
    return 0;
  }

  double E = fP4->Energy();
  double M = ( (mass_from_pdg) ? this->Mass() : fP4->M() );
  double K = E-M;

  K = TMath::Max(K,0.);
  return K;
}
//___________________________________________________________________________
int GHepParticle::Z(void) const
{
// Decoding Z from the PDG code

  if(!fIsNucleus) return -1;
  else
       return pdg::IonPdgCodeToZ(fPdgCode);
}
//___________________________________________________________________________
int GHepParticle::A(void) const
{
// Decoding A from the PDG code

  if(!fIsNucleus) return -1;
  else
       return pdg::IonPdgCodeToA(fPdgCode);
}
//___________________________________________________________________________
TLorentzVector * GHepParticle::GetP4(void) const 
{ 
// see GHepParticle::P4() for a method that does not create a new object and
// transfers its ownership 

  if(fP4) {
     TLorentzVector * p4 = new TLorentzVector(*fP4); 
     LOG("GHepParticle", pDEBUG) 
                    << "Return vp = " << utils::print::P4AsShortString(p4);
     return p4;
  } else {
    LOG("GHepParticle", pWARN) << "NULL 4-momentum TLorentzVector";
    return 0;
  }
}
//___________________________________________________________________________
TLorentzVector * GHepParticle::GetX4(void) const 
{ 
// see GHepParticle::X4() for a method that does not create a new object and
// transfers its ownership

  if(fX4) {
     TLorentzVector * x4 = new TLorentzVector(*fX4); 
     LOG("GHepParticle", pDEBUG) 
                          << "Return x4 = " << utils::print::X4AsString(x4);
     return x4;
  } else {
    LOG("GHepParticle", pWARN) << "NULL 4-position TLorentzVector";
    return 0;
  }
}
//___________________________________________________________________________
void GHepParticle::SetPdgCode(int code)
{
// Always set PDG code through this method to make sure that the fIsNucleus
// flag is set.

  fPdgCode = code;
  this->AssertIsKnownParticle();

  // GENIE uses the MINOS PDG extension for ions: PDG Code = [1AAAZZZ000]

  fIsNucleus = pdg::IsIon(code);

  // GENIE uses PDG codes in the range 1111111001 - 1111111999 for fake
  // particles + the rootino with pdg-code = 0

  fIsFake = ( (code == 0) || (code>1111111000 && code < 1111112000) );
}
//___________________________________________________________________________
void GHepParticle::SetMomentum(const TLorentzVector & p4)
{
  if(fP4)
      fP4->SetPxPyPzE( p4.Px(), p4.Py(), p4.Pz(), p4.Energy() );
  else
      fP4 = new TLorentzVector(p4);
}
//___________________________________________________________________________
void GHepParticle::SetMomentum(double px, double py, double pz, double E)
{
  if(fP4)
      fP4->SetPxPyPzE(px, py, pz, E);
  else
      fP4 = new TLorentzVector(px, py, pz, E);
}
//___________________________________________________________________________
void GHepParticle::SetPosition(const TLorentzVector & v4)
{
  this->SetPosition(v4.X(), v4.Y(), v4.Z(), v4.T());
}
//___________________________________________________________________________
void GHepParticle::SetPosition(double x, double y, double z, double t)
{
  LOG("GHepParticle", pDEBUG) 
            << "Setting position to (x = " << x << ", y = " 
                               << y << ", z = " << z << ", t = " << t << ")";

  if(fX4) fX4->SetXYZT(x,y,z,t);
  else    fX4 = new TLorentzVector(x,y,z,t);
}
//___________________________________________________________________________
void GHepParticle::SetEnergy(double E)
{
  this->SetMomentum(this->Px(), this->Py(), this->Pz(), E);
}
//___________________________________________________________________________
void GHepParticle::SetPx(double px)
{
  this->SetMomentum(px, this->Py(), this->Pz(), this->E());
}
//___________________________________________________________________________
void GHepParticle::SetPy(double py)
{
  this->SetMomentum(this->Px(), py, this->Pz(), this->E());
}
//___________________________________________________________________________
void GHepParticle::SetPz(double pz)
{
  this->SetMomentum(this->Px(), this->Py(), pz, this->E());
}
//___________________________________________________________________________
bool GHepParticle::IsOnMassShell(void) const
{
  this->AssertIsKnownParticle();

  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);

  double Mpdg = p->Mass();
  double M4p  = (fP4) ? fP4->M() : 0.;

//  return utils::math::AreEqual(Mpdg, M4p);

  return (TMath::Abs(M4p-Mpdg) < kOffShellDm);
}
//___________________________________________________________________________
bool GHepParticle::IsOffMassShell(void) const
{
  return (! this->IsOnMassShell());
}
//___________________________________________________________________________
void GHepParticle::Init(void)
{
  fIsNucleus = false;
  fIsFake    = true;

  fPdgCode = 0;
  fStatus  = kIStUndefined;

  fFirstMother   = -1;
  fLastMother    = -1;
  fFirstDaughter = -1;
  fLastDaughter  = -1;

  fP4 = new TLorentzVector(0,0,0,0);
  fX4 = new TLorentzVector(0,0,0,0);
}
//___________________________________________________________________________
void GHepParticle::CleanUp(void)
{
// deallocate memory

  if(fP4) delete fP4;
  if(fX4) delete fX4;
  fP4 = 0;
  fX4 = 0;
}
//___________________________________________________________________________
void GHepParticle::Reset(void)
{
// deallocate memory + initialize

  this->CleanUp();
  this->Init();
}
//___________________________________________________________________________
void GHepParticle::Clear(Option_t * option)
{
// implement the Clear(Option_t *) method so that the GHepParticle when is a
// member of a GHepRecord, gets deleted properly when calling TClonesArray's
// Clear("C")

  this->CleanUp();
}
//___________________________________________________________________________
void GHepParticle::Print(ostream & stream) const
{
  stream << "\n |";
  stream << setfill(' ') << setw(14) << this->Name()           << " | ";
  stream << setfill(' ') << setw(14) << this->PdgCode()        << " | ";
  stream << setfill(' ') << setw(6)  << this->Status()         << " | ";
  stream << setfill(' ') << setw(3)  << this->FirstMother()    << " | ";
  stream << setfill(' ') << setw(3)  << this->LastMother()     << " | ";
  stream << setfill(' ') << setw(3)  << this->FirstDaughter()  << " | ";
  stream << setfill(' ') << setw(3)  << this->LastDaughter()   << " | ";
  stream << setiosflags(ios::fixed)  << setprecision(3);
  stream << setfill(' ') << setw(6)  << this->Px()             << " | ";
  stream << setfill(' ') << setw(6)  << this->Py()             << " | ";
  stream << setfill(' ') << setw(6)  << this->Pz()             << " | ";
  stream << setfill(' ') << setw(6)  << this->E()              << " | ";
  stream << setfill(' ') << setw(6)  << this->Mass()           << " | ";
}
//___________________________________________________________________________
void GHepParticle::Print(Option_t * /*opt*/) const
{
// implement the TObject's Print(Option_t *) method

  this->Print(cout);
}
//___________________________________________________________________________
bool GHepParticle::Compare(const GHepParticle * p) const
{
// Do the comparisons in steps & put the ones that cost most
// in the inner-most {}

  bool same_particle = (this->fPdgCode == p->fPdgCode);
  bool same_status   = (this->fStatus  == p->fStatus );

  if( !same_particle || !same_status )  return false;
  else {
     if ( ! this->CompareFamily(p) )    return false;
     else {
       if( ! this->CompareMomentum(p) ) return false;
       else                             return true;
     }
  }
}
//___________________________________________________________________________
bool GHepParticle::CompareFamily(const GHepParticle * p) const
{
  bool same_family  = (
          this->fFirstMother == p->fFirstMother &&
                  this->fLastMother  == p->fLastMother &&
                          this->fFirstDaughter == p->fFirstDaughter &&
                                   this->fLastDaughter  == p->fLastDaughter
  );
  return same_family;
}
//___________________________________________________________________________
bool GHepParticle::CompareMomentum(const GHepParticle * p) const
{
  double dE  = TMath::Abs( this->E()  - p->E()  );
  double dPx = TMath::Abs( this->Px() - p->Px() );
  double dPy = TMath::Abs( this->Py() - p->Py() );
  double dPz = TMath::Abs( this->Pz() - p->Pz() );

  bool same_momentum =
       (dE < kPCutOff && dPx < kPCutOff && dPy < kPCutOff && dPz < kPCutOff);

  return same_momentum;
}
//___________________________________________________________________________
void GHepParticle::Copy(const GHepParticle & particle)
{
  this->SetStatus        (particle.Status());
  this->SetFirstMother   (particle.FirstMother());
  this->SetLastMother    (particle.LastMother());
  this->SetFirstDaughter (particle.FirstDaughter());
  this->SetLastDaughter  (particle.LastDaughter());
  this->SetMomentum      (*particle.P4());
  this->SetPosition      (*particle.X4());
  this->SetPdgCode       (particle.PdgCode());
}
//___________________________________________________________________________
void GHepParticle::AssertIsKnownParticle(void) const
{
  TParticlePDG * p = PDGLibrary::Instance()->Find(fPdgCode);

  if(!p) {
    LOG("GHepParticle", pFATAL)
       << "Found unknown particle [pdg-code = " << fPdgCode << "] !! "
              << "You might get a Nobel prize but now get a core dump. Bye!";
    assert(p);
  }
}
//___________________________________________________________________________
bool GHepParticle::operator == (const GHepParticle & p) const
{
  return (this->Compare(&p));
}
//___________________________________________________________________________
GHepParticle & GHepParticle::operator = (const GHepParticle & p)
{
  this->Copy(p);
  return (*this);
}
//___________________________________________________________________________

