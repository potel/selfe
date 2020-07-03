#include <fstream>
#include <iostream>
#include <armadillo>
#include <ctime>
using namespace arma;
using namespace std;
#include "selfe.h"
#include "structs.h"
#include "declarations.h"
ofstream misc1("misc1.txt");
ofstream misc2("misc2.txt");
ofstream misc3("misc3.txt");
ofstream misc4("misc4.txt");
ofstream misc5("misc5.txt");
ofstream misc6("misc6.txt");
ofstream misc7("misc7.txt");
ofstream misc8("misc8.txt");
ofstream output("output.txt");
ofstream phase("phase-shifts.txt");
ofstream cx("absorption.txt");

int main(int argc,char* argv[]){
	parameters *parm=new parameters;
	cout<<" Reading parameters from "<<argv[1]<<endl;
	const char* input=argv[1];
    ReadParameters(input,parm);
    TestBuild(parm);
}
void TestBuild(parameters *parm)
{
  cx_vec poles,c;
  int n,indx_mean_field,indx_coupling;
  double k,sigma_r,absorption,energy,energy_end,de;
  complejo phase_shift;
  cx_mat funpot=zeros<cx_mat>(parm->lagpoints,parm->lagpoints);
  nlpotential *V=new nlpotential(parm->points,parm->radius,parm->lagpoints);
  poles.zeros(parm->numpoles);
  c.zeros(parm->lagpoints);
  lagrange lag(parm->points,parm->lagpoints,parm->lag_box);
  lag.LagrangeBasis();
  vector <state> overlaps(parm->numpoles,state(parm->points,&lag));
  state *st=new state(parm->points,&lag);
  energy_end=parm->delta_e*parm->numpoles+5.*parm->gamma;
  cout<<"Upper energy limit: "<<energy_end<<"\n";
  de=parm->gamma/10;
  cout<<"Energy resolution: "<<de<<"\n";
  for (n=0;n<2;n++)
    {
      if(parm->mean_field==parm->potentials[n].id) indx_mean_field=n;
      if(parm->coupling==parm->potentials[n].id) indx_coupling=n;
      parm->potentials[n].Fill(&lag);
    }
  V->type="locnloc";
  V->pot_reduced=parm->potentials[indx_mean_field].pot_reduced;
  V->r_reduced=parm->potentials[indx_mean_field].r_reduced;
  lag.LagrangeBasis();
  for(n=0;n<parm->numpoles;n++)
    {
      poles(n)=parm->delta_e*(n+1)+parm->gamma*I;
      overlaps[n].energy=real(poles(n));
      overlaps[n].l=0;
      overlaps[n].j=0.5;
      overlaps[n].GenScatt(1.,0.,&parm->potentials[indx_mean_field],&lag,c);
      overlaps[n].wf_reduced=overlaps[n].wf_reduced%parm->potentials[indx_coupling].pot_reduced;
    }
  V->poles=poles;
  V->overlaps=overlaps;
  st->l=0;
  st->j=0.5;
  phase_shift=0.;
  st->energy=0.5;
  for(energy=de;energy<energy_end;energy+=de)
    {
      st->energy=energy;
      phase_shift=st->GenScatt(1.,0.,V,&lag,c);
      absorption=Absorption(st,V,1.,&lag);
      phase<<energy<<"  "<<(real(phase_shift))<<"  "<<(imag(phase_shift))<<"  "<<abs(phase_shift)<<"\n";
      cx<<energy<<"  "<<absorption<<"\n";
    }
}

///////////////////////////////////////////////////////////////////////
//  Generate scattering wavefunction, Numerov numerical integration  //
//      (only valid for local potential)                             // 
///////////////////////////////////////////////////////////////////////
complejo state::GenScatt(double q1q2,double mass,optical* pot){
  int i, i_1, i_2,status1,status2,points;
	double hbarx, dd, radio_1, radio_2, q, x, y, ex1, ex2, etac,delta_r,
      spinorbit,radio_match;
	complejo delta, derivada_log, fu1, fu2, factor;
    points=pot->r.n_elem;
	complejo *potential=new complejo[points];
	gsl_complex deltagsl;
	gsl_sf_result F1, G1, F2, G2, Fp, Gp;
	if(energy<0.) Error("Negative energy while trying to compute a scattering state\n");
    r=pot->r;
    radius=r(points-1);
	delta_r=r(points-1)-r(points-2);
	radio_match=radius-10.*delta_r;
	hbarx=HC*HC/(2.*AMU*mass);
	dd=delta_r*delta_r/hbarx;
	q=sqrt(energy/hbarx);
	etac=q1q2*mass*E2HC*AMU/(HC*q);
    for (i=0;i<points-1;i++) {
      potential[i]=pot->pot(i)+(l*(l+1.))*hbarx /(r(i)*r(i));
      //cout<<potential[i]<<"  "<<pot->pot[i]<<"\n";
      //exit(0);
	}
	wf(0)=1.e-10;
	wf(1)=(2.*(1.-0.416666667*dd*(-potential[0]+energy))*wf(0))/
			(1.+0.083333333*dd*(-potential[1]+energy));
	for (i=1;i<points-1;i++) {
      wf(i+1)=(2.*(1.-0.416666667*dd*(-potential[i]+energy))
               *wf(i)-(1.+0.083333333*dd*(-potential[i-1]+ energy))*wf(i-1))
				/(1.+0.083333333*dd*(-potential[i+1]+energy));
	}
	// Computing the phase shift delta
	radio_1=radio_match;
	i_1=int(ceil(radio_1/delta_r))-1;
	radio_1=delta_r*(i_1 + 1.);
	fu1=wf(i_1);
	radio_2=radio_1+2.*delta_r;
	i_2=int(ceil(radio_2/delta_r))-1;
	radio_2=delta_r*(i_2+1.);
	fu2=wf(i_2);
	derivada_log=fu1/fu2;
	status1=gsl_sf_coulomb_wave_FG_e(etac,q*radio_1,l,0,&F1,&Fp,&G1,&Gp,&ex1,&ex2);
	status2=gsl_sf_coulomb_wave_FG_e(etac,q*radio_2,l,0,&F2,&Fp,&G2,&Gp,&ex1,&ex2);
	x=real((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
	y=imag((F1.val-F2.val*derivada_log)/(-G1.val+G2.val*derivada_log));
	GSL_SET_COMPLEX(&deltagsl,x,y);
    //cout<<fu1<<"  "<<fu2<<"\n";
    //exit(0);
	delta=(gsl_complex_arctan(deltagsl).dat[0]+I*gsl_complex_arctan(deltagsl).dat[1]); // desfase
	factor=exp(I*(delta))*(cos(delta)*F1.val+sin(delta)*G1.val)/fu1;
	for (i=0;i<points;i++) {
      wf(i)=factor*wf(i);
      //misc1<<r(i)<<"  "<<real(wf(i))<<"\n";
      if(status1 || status2) wf(i)=0.;
	}
    //exit(0);
	delete[] potential;
	return delta;  
}


void Error(const string message)
{
  cout<<message<<"\n";
  output<<message<<"\n";
}
/*
 Generate scattering wavefunction computed with the Lagrange mesh method,
 and returns the phase shift. The resulting wavefunction is defined only in
 the mesh points. Valid for local and non local potentials.
      Potentials have to be defined at the mesh points, defined in lag.
      c returns the coefficients of the Lagrange expansion of the wavefunction.
*/
complejo state::GenScatt(double mass,double q1q2, nlpotential *potential, lagrange *lag,cx_vec &c)
{       
  double hbarx,central,part3,part4,ri,rj,q,etac,
    exp_F,exp_G,delta_a;
  complejo pot,Rmatrix,Hp,Hm,Hmp,Hpp,S,phase_shift,factor,vloc;
  int i,j;
  gsl_sf_result F,G,Fp,Gp;
  if(energy==0.) energy=0.000001;
  hbarx=HC*HC/(2.*AMU*mass);
  q=sqrt(energy/hbarx);
  etac=q1q2*mass*E2HC*AMU/(HC*q);
  cx_mat TLmatrix=zeros<cx_mat>(lag->N,lag->N);
  cx_mat Vmatrix=zeros<cx_mat>(lag->N,lag->N);
  cx_mat Gmatrix=zeros<cx_mat>(lag->N,lag->N);
  vec basis_a(lag->N);
  int numpoles,n;
  if(potential->type=="nloc" || potential->type=="locnloc") numpoles=potential->poles.n_elem;
  //cout<<"num poles: "<<numpoles<<"\n";
  basis_a=lag->basis.row(lag->basis.n_rows-1).t(); // vector of length N with the value of each Lagrange function at a
  // Initialization of Coulomb functions at a
  gsl_sf_coulomb_wave_FG_e(etac,q*lag->a,l,0,&F,&Fp,&G,&Gp,&exp_F,&exp_G);
  Hp=(G.val+I*F.val);
  Hm=(G.val-I*F.val);
  Hpp=q*(Gp.val+I*Fp.val);
  Hmp=q*(Gp.val-I*Fp.val);
  // Kinetic energy, Bloch term and energy eigenvalue
  for(i=0;i<lag->N;i++)
    {
      TLmatrix(i,i)=hbarx*((4.*lag->N*lag->N+4.*lag->N+3.0)*lag->x[i]*(1.-lag->x[i])-6.*lag->x[i]+1.)/
        (3.*lag->a*lag->a*(lag->x[i]*lag->x[i])*((1.-lag->x[i])*(1.-lag->x[i])))-energy;
      for(j=i+1;j<lag->N;j++)
        {
          part3=pow(-1.,i+j)/(lag->a*lag->a*sqrt(lag->x[i]*lag->x[j]*(1.-lag->x[i])*(1.-lag->x[j])));
          part4=(lag->N*lag->N*1.+lag->N*1.0+1.0+(lag->x[i]+lag->x[j]-2.*lag->x[i]*lag->x[j])/
                 ((lag->x[i]-lag->x[j])*(lag->x[i]-lag->x[j]))-1./(1.-lag->x[i])-1./(1.-lag->x[j]));
          TLmatrix(i,j)=hbarx*part3*part4;
          TLmatrix(j,i)=TLmatrix(i,j);
        }
    }
  // Potential (multiplied by rj*rj), including central potential and energy
  for(i=0;i<lag->N;i++)
    {
      ri=lag->a*lag->x[i];
      Vmatrix(i,i)=potential->pot_reduced(i)+hbarx*l*(l+1.)/(ri*ri);  // central potential and energy
      if(potential->type=="nloc" || potential->type=="locnloc")
        {
          for(j=0;j<lag->N;j++)
            {
              rj=lag->a*lag->x[j];
              potential->nlpot_reduced(i,j)=0.;
              for(n=0;n<numpoles;n++)
                {
                  potential->nlpot_reduced(i,j)-=1./(energy-potential->poles(n))*
                    conj(potential->overlaps[n].wf_reduced(i))*potential->overlaps[n].wf_reduced(j);
                }
              Vmatrix(i,j)+=lag->a*sqrt(lag->w[i]*lag->w[j]/4.)*potential->nlpot_reduced(i,j);
            }
        }
    }
  Gmatrix=inv(TLmatrix+Vmatrix); 
  Rmatrix=hbarx*dot(basis_a.t(),Gmatrix*basis_a)/lag->a;
  S=(Hm-lag->a*Rmatrix*Hmp)/(Hp-lag->a*Rmatrix*Hpp);
  phase_shift=-I*log(S)/2.;
  c=I*hbarx*(Hmp-S*Hpp)*Gmatrix*basis_a/2.;  // c coefficients (see eq. (3.13))
  for(i=0;i<lag->N;i++)
    {
      r_reduced(i)=lag->a*lag->x[i];
      wf_reduced(i)=c(i)/sqrt(lag->a*lag->w[i]/2.);
    }
  return phase_shift;
}

// Lagrange basis initializer, from 0 to a
void lagrange::LagrangeBasis()
{
  const double EPS=1.0e-10;
  int m,j,i;
  double z1,z,pp,p3,p2,p1,value,val,val_i;
  double F[N],dF[N];
  a1=0.;
  m=(N+1)/2;
  for (i=0;i<m;i++) {
    z=cos(PI*(i+0.75)/(N+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=0;j<N;j++) {
        p3=p2;
        p2=p1;
        p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
      }
      pp=N*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    x[i]=(-z+1.0)/2.0;
    x[N-1-i]=(z+1.0)/2.0;
    w[i]=2./((1.0-z*z)*pp*pp);
    w[N-1-i]=w[i];
  }
  for(m=0;m<r.size();m++)
    {
      for(i=1;i<=N;i++)
        {
          for(m=0;m<r.size();m++)
            {
              value=2.*r[m]/a-1.;
              if (value>1.) value=1.;
              if (value<-1.) value=-1.;
              gsl_sf_legendre_Pl_deriv_array(N,value,F,dF);
              for(i=1;i<=N;i++)
                {
                  val=(2.*r[m]/a)*dF[N]-a*x[i-1]*F[N]/(r[m]-a*x[i-1]);
                  val_i=(2./a)*dF[N]-F[N]/(r[m]-a*x[i-1]);
                  basis(m,i-1)= pow(-1.0,N+i)*(r[m]/(a*x[i-1]))*sqrt(a*x[i-1]*(1.-x[i-1]))*F[N]/(r[m]-a*x[i-1]);
                  basis_i(m,i-1)= pow(-1.0,N+i)*sqrt(a*x[i-1]*(1.-x[i-1]))*F[N]/(r[m]-a*x[i-1]);
                  del_basis(m,i-1)= pow(-1.0,N+i)*sqrt(a*x[i-1]*(1.-x[i-1]))*val/(a*x[i-1]*(r[m]-a*x[i-1]));
                  del_basis_i(m,i-1)= pow(-1.0,N+i)*sqrt(a*x[i-1]*(1.-x[i-1]))*val_i/(r[m]-a*x[i-1]);
                }
            }         
        }
    }
}


// Lagrange basis initializer from initial to initial+a
void lagrange::LagrangeBasis(double initial)
{
  
  const double EPS=1.0e-10;
  int m,j,i;
  double z1,z,pp,p3,p2,p1,value,value_i,val,val_i;
  double F[N],dF[N];
  double F_i[N],dF_i[N];
  a1=initial;
  a2=initial+a;
  m=(N+1)/2;
  for (i=0;i<m;i++) {
    z=cos(PI*(i+0.75)/(N+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=0;j<N;j++) {
        p3=p2;
        p2=p1;
        p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
      }
      pp=N*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    x[i]=(-z+1.0)/2.0;
    x[N-1-i]=(z+1.0)/2.0;
    w[i]=2./((1.0-z*z)*pp*pp);
    w[N-1-i]=w[i];
  }

  for(m=0;m<r.size();m++)
    {
      value=2.*r[m]/a-1.;
      value_i=(2.*(r[m]-a1)-a)/a;
      if (value>1.) value=1.;
      if (value<-1.) value=-1.;
      if (value_i>1.) value_i=1.;
      if (value_i<-1.) value_i=-1.;
      gsl_sf_legendre_Pl_deriv_array(N,value,F,dF);
      gsl_sf_legendre_Pl_deriv_array(N,value_i,F_i,dF_i);
      for(i=1;i<=N;i++)
        {
          val=(2.*r[m]/a)*dF[N]-a*x[i-1]*F[N]/(r[m]-a*x[i-1]);
          val_i=(2./a)*dF_i[N]-F_i[N]/(r[m]-a*x[i-1]-a1);
          basis(m,i-1)= pow(-1.0,N+i)*(r[m]/(a*x[i-1]))*sqrt(a*x[i-1]*(1.-x[i-1]))*F[N]/(r[m]-a*x[i-1]);
          basis_i(m,i-1)= pow(-1.0,N+i)*sqrt(a*x[i-1]*(1.-x[i-1]))*F_i[N]/(r[m]-a*x[i-1]-a1);
          del_basis(m,i-1)= pow(-1.0,N+i)*sqrt(a*x[i-1]*(1.-x[i-1]))*val/(a*x[i-1]*(r[m]-a*x[i-1]));
          del_basis_i(m,i-1)= pow(-1.0,N+i)*sqrt(a*x[i-1]*(1.-x[i-1]))*val_i/(r[m]-a*x[i-1]-a1);
        }
    }
}
void PolarPotential(double energy,lagrange *lag,nlpotential *v)
{
  int numpoles,n,i,j;
  numpoles=v->poles.n_elem;
  for(i=0;i<lag->N;i++)
    {
      v->nlpot_reduced(i,i)=v->pot_reduced(i);
      for(j=0;j<lag->N;j++)
        {
          if(i!=j) v->nlpot_reduced(i,j)=0.;
          for(n=0;n<numpoles;n++)
            {
              v->nlpot_reduced(i,j)-=1./(energy-v->poles(n))*
                conj(v->overlaps[n].wf_reduced(i))*v->overlaps[n].wf_reduced(j);
            }
        }
    }
}
void Expand(state *st,lagrange *lag,cx_vec c)
{
  int i,j;
  double r;
  st->radius=lag->a;
  for(i=0;i<lag->r.size();i++)
    {
      r=lag->r[i];
      st->wf(i)=0.;
      st->r(i)=r;
      for(j=0;j<lag->N;j++)
        {
          st->wf(i)+=c(j)*lag->basis(i,j);
        }
    }
}



void Write(vec r,vec f,ofstream *file)
{
  int n;
  for(n=0;n<r.n_elem;n++)
    {
      *file<<r(n)<<"  "<<f(n)<<"\n";
    }
}
void Write(vec r,cx_mat f,ofstream *file)
{
  int n;
  for(n=0;n<r.n_elem;n++)
    {
      *file<<r(n)<<"  "<<real(f(n,n))<<"  "<<imag(f(n,n))<<"  "<<abs(f(n,n))<<"\n";
    }
}

void Write(vec r,cx_vec f,ofstream *file)
{
  int n;
  for(n=0;n<r.n_elem;n++)
    {
      *file<<r(n)<<"  "<<real(f(n))<<"  "<<imag(f(n))
          <<"  "<<abs(f(n))<<"\n";
    }
}
double Absorption(state *st,nlpotential *v,double mass,lagrange *lag)
{
  int i,j,dim;
  static int e;
  e++;
  double k,a;
  complejo sum;
  k=sqrt(2.*mass*AMU*st->energy)/HC;
  a=2.*mass*AMU/(HC*HC*k);
  dim=st->r_reduced.n_elem;
  sum=0.;
  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
        {
          sum+=a*conj(st->wf_reduced(i))*imag(v->nlpot_reduced(i,j))*
            st->wf_reduced(j)*lag->w[i]*lag->w[j];
          // sum+=a*conj(st->wf_reduced(i))*imag(v->nlpot_reduced(i,j))*
          //   st->wf_reduced(j)*lag->w[i]*lag->w[j];
        }
    }
  //  misc5<<e<<"  "<<real(sum)<<"  "<<imag(v->nlpot_reduced(1,1))<<"\n";
  return abs(sum);
}
void Converge(nlpotential *v,int iter,lagrange *lag)
{
  int numpoles,n,m;
  optical Uopt(v->r.n_elem,v->radius);
  numpoles=v->poles.n_elem;
  complejo phase_shift;
  cx_vec c;
  Uopt.r_reduced=v->r_reduced;
  Uopt.pot_reduced=v->pot_reduced;
  c.zeros(lag->N);
  for(n=0;n<numpoles;n++)
    {
      v->overlaps[n].energy=real(v->poles[n]);
      v->overlaps[n].l=0;
      v->overlaps[n].j=0.5;
      v->overlaps[n].GenScatt(1.,0.,&Uopt,lag,c);
      v->overlaps[n].wf_reduced=v->overlaps[n].wf_reduced%v->pot_reduced;
    }
  for(n=0;n<iter;n++)
    {
      misc1<<n;
      for(m=0;m<numpoles;m++)
        {
          phase_shift=v->overlaps[m].GenScatt(1.,0.,v,lag,c);
          v->overlaps[m].wf_reduced=v->overlaps[m].wf_reduced%v->pot_reduced;
          misc1<<"  "<<real(phase_shift);
        }
      misc1<<"\n";
    }
}
void ReadParameters(const char *fname,parameters *x)
{
  char aux[500];
  int n,pots;
  ifstream fp;
  fp.open(fname);
  if (!fp) Error("Couldn't open parameter file");
  x->points=0;   // number of points for the representation of functions
  x->lagpoints=0; // number of points of the Lagrange mesh
  x->radius=0;    // radius for the representation of functions (fm)
  x->lag_box=0;   // size of the Lagrange mesh (fm)

  // Parameters for the definition of resonances 
  x->numpoles=0;  // number of resonances g
  x->gamma=0;   // width of the resonances (imaginary part of the poles) (fm)
  x->delta_e=0;   // distance between resonances (MeV)

  // Index to mean field potential and neutron-nucleus coupling
  x->mean_field=0;
  x->coupling=0;
  pots=0;
  while(fp.getline(aux,500))
	{
      ReadParD(aux,"points",&(x->points));
      ReadParD(aux,"lagpoints",&(x->lagpoints));
      ReadParD(aux,"numpoles",&(x->numpoles));
      ReadParD(aux,"mean_field",&(x->mean_field));
      ReadParD(aux,"coupling",&(x->coupling));
      ReadParF(aux,"radius",&(x->radius));
      ReadParF(aux,"lag_box",&(x->lag_box));
      ReadParF(aux,"gamma",&(x->gamma));
      ReadParF(aux,"delta_e",&(x->delta_e));
      ReadParF(aux,"radius",&(x->radius));
	}
  fp.clear();
  fp.seekg(0, std::ios_base::beg );
  while(fp.getline(aux,500))
	{
      if (!pots) pots=ReadPotentials(aux,"Potentials-Start",x->potentials,&fp,
                                     x->points,x->lagpoints,x->radius);
	}
  for(n=0;n<2;n++)
    {
      x->potentials[n].Fill(x->points,x->radius);
    }
  fp.close();

}
int ReadPotentials(char *s,const char key[100],vector <optical> &pot,
                   ifstream *fp,int points,int lagpoints,double radius)
{
  int l,r,i,l2,l3,id;
  const char fin[]="Potentials-End";
  const char flag[]="***************";
  char aux[500]="\0";
  l = strlen(key);
  l2 = strlen(fin);
  l3 = strlen(flag);
  double V,W,Vso,Wso,Wr,Wd,radV,radW,radcoul,
    aV,aW,radso,aso,radWd,aWd;
  i=0;
  if (!strncmp(s,key,l))
    {
      while(strncmp(aux,fin,l2))
        {
          fp->getline(aux,500);
          ReadParD(aux,"id",&(id));
          ReadParF(aux,"RealVolume",&(V));
          ReadParF(aux,"ImaginaryVolume",&(W));
          ReadParF(aux,"RealSpinOrbit",&(Vso));
          ReadParF(aux,"ImaginarySpinOrbit",&(Wso));
          ReadParF(aux,"RealSurface",&(Wr));
          ReadParF(aux,"ImaginarySurface",&(Wd));
          ReadParF(aux,"RadiusRealVolume",&(radV));
          ReadParF(aux,"RadiusImaginaryVolume",&(radW));
          ReadParF(aux,"RadiusCoulomb",&(radcoul));
          ReadParF(aux,"DifusivityRealVolume",&(aV));
          ReadParF(aux,"DifusivityImaginaryVolume",&(aW));
          ReadParF(aux,"RadiusSpinOrbit",&(radso));
          ReadParF(aux,"DifusivitySpinOrbit",&(aso));
          ReadParF(aux,"RadiusImaginarySurface",&(radWd));
          ReadParF(aux,"DifusivityImaginarySurface",&(aWd));
          if(!strncmp(aux,flag,3))
            {
              pot.push_back(optical(points,radius,lagpoints));
              pot[i].id=id;
              pot[i].V=V;
              pot[i].W=W;
              pot[i].Vso=Vso;
              pot[i].Wso=Wso;
              pot[i].Wr=Wr;
              pot[i].Wd=Wd;
              pot[i].radV=radV;
              pot[i].radW=radW;
              pot[i].radcoul=radcoul;
              pot[i].aV=aV;
              pot[i].aW=aW;
              pot[i].radso=radso;
              pot[i].aso=aso;
              pot[i].radWd=radWd;
              pot[i].aWd=aWd;
              i++;
            }
        }
		cout<<i<<" Potentials have been read"<<endl;
		fflush(stdout);
		return i;
	}
	return 0;
}
/*****************************************************************************
Reads decimal number
 *****************************************************************************/
void ReadParD(char *s,const char key[20], int *par)
{
	int l,r;
	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %d",par);
		if (r<1) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros (r=%d)\n",key,r); exit(1); }
	}
	fflush(stderr);
}
/*****************************************************************************
Reads float
 *****************************************************************************/
void ReadParF(char *s,const char key[20],double *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %lf",par);
		if (r<1) { fprintf(stderr,"ERROR: No puede leerse %s en el fichero de parametros\n",key); exit(1); }
	}
	fflush(stderr);
}
/*****************************************************************************
Reads char 
 *****************************************************************************/
void ReadParS(char *s,const char key[20], char *par)
{
	int l,r;

	l = strlen(key);
	if (!strncmp(s,key,l))
	{
		r = sscanf(s,"%*s %s",par);
		if (r<1) { fprintf(stderr,"%s ",key); *par='\0'; }
	}
	fflush(stderr);
}
/*****************************************************************************
Reads  string
 *****************************************************************************/
void ReadParStr(string s,const char key[20], string par)
{
	int l,r;
    string str=key;
	l = strlen(key);
	if (s==key)
	{
		par=s;
        cout<<str<<" = "<<par<<endl;
	}
}
void optical::Fill(int points,double radius)
{
  int n;
  double delta_r;
  delta_r=radius/double(points);
  for(n=0;n<points;n++)
	{
      r(n)=delta_r*(n+1.);
      Nuclear(n)=-V/(1.+exp((r(n)-radV)/aV))-I*W/
        (1.+exp((r(n)-radW)/aW))-4.*I*Wd*
        exp((r(n)-radWd)/aWd)/((1.+exp((r(n)-radWd)/aWd))
                               *(1.+exp((r(n)-radWd)/aWd)))
        -Wr*exp((r(n)-radV)/aV)/((1.+exp((r(n)-radV)/aV))
                               *(1.+exp((r(n)-radV)/aV)));
	}
  pot=Nuclear;
}

void optical::Fill(lagrange *lag)
{
  int n;
  for(n=0;n<lag->N;n++)
	{
      r_reduced(n)=lag->x[n]*lag->a;
      pot_reduced(n)=-V/(1.+exp((r_reduced(n)-radV)/aV))-I*W/
        (1.+exp((r_reduced(n)-radW)/aW))-4.*I*Wd*
        exp((r_reduced(n)-radWd)/aWd)/((1.+exp((r_reduced(n)-radWd)/aWd))
                                       *(1.+exp((r_reduced(n)-radWd)/aWd)))
        -Wr*exp((r_reduced(n)-radV)/aV)/((1.+exp((r_reduced(n)-radV)/aV))
                               *(1.+exp((r_reduced(n)-radV)/aV)));
	}
}
