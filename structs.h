class estado;
class lagrange;
class state;
class nlpotential {
 public:
  int id;
  double radius;
  string type; // type either 'loc', 'nloc', 'locnloc'
  vec r;
  vec r_reduced;
  // reduced variables are defined at Gaussian integration points.
  // they are defined along with a lagrange object
  cx_vec pot;
  cx_vec pot_reduced;
  cx_mat nlpot;
  cx_mat nlpot_reduced; 
  string file;
  cx_vec poles;     // complex poles of S-matrix (resonances)
  vector <state> overlaps;   // overlaps for the pole decomposition
  nlpotential () {}
  nlpotential (int points,double x1):radius{x1}
  {
    double delta_r;
    int n;
    r.zeros(points);
    pot.zeros(points);
    nlpot.zeros(points,points);
    delta_r=radius/double(points);
    for(n=0;n<points;n++)
      {
        r(n)=delta_r*(n+1.);
      }
    type="locnloc";
  }
nlpotential (int points,double x1,int lagpoints):radius{x1}
  {
    double delta_r;
    int n;
    r.zeros(points);
    pot.zeros(points);
    nlpot.zeros(points,points);
    pot_reduced.zeros(lagpoints);
    nlpot_reduced.zeros(lagpoints,lagpoints);
    r_reduced.zeros(lagpoints);
    delta_r=radius/double(points);
    for(n=0;n<points;n++)
      {
        r(n)=delta_r*(n+1.);
      }
    type="locnloc";
  }
};

class estado {
 public:
  int id;
  nlpotential* pot;
  nlpotential* pot_reduced;
  int puntos;
  double radio;
  int l;
  double j;
  int nodos;
  double r[MAX_PTS];
  complejo wf[MAX_PTS];
  double energia;
  double spec;
  char file[100];
  estado() {};
  
};

class lagrange
{
 public:
  vector_dbl x;   // Lagrange points
  vector_dbl w;   // Lagrange weights
  vector_dbl r;   // radial grid (length pts), from 0 to a
  int N;          // size of Lagrange basis
  mat basis;     // Lagrange basis functions (pts x N)
  mat basis_i;     // Lagrange irregular basis functions (pts x N)
  mat del_basis;     // derivative of Lagrange basis functions (pts x N)
  mat del_basis_i;     // derivative of Lagrange irregular basis functions (pts x N)

  double a;   // size of box;
  double a1;  // initial point for propagation method
  double a2;   //final point for propagation method
  lagrange() {};
  lagrange (int pts,int Nl, double box):N{Nl},a{box}
  {
    int i;
    double rn;
    double step=box/double(pts);
    for(i=0;i<N;i++)
      {
        x.push_back(0.);
        w.push_back(0.);
      }
    rn=step;
    for(i=0;i<pts;i++)
      {
        r.push_back(rn);
        rn+=step;
      }
    basis.zeros(pts,N);
    basis_i.zeros(pts,N);
    del_basis.zeros(pts,N);
    del_basis_i.zeros(pts,N);
    a1=0.;
  }
  lagrange (int pts,int Nl,double box,double initial):N{Nl},a{box},a1{initial}
  {
    a2=a1+a;
    int i;
    double rn;
    double step=box/double(pts);
    for(i=0;i<N;i++)
      {
        x.push_back(0.);
        w.push_back(0.);
      }
    rn=step;
    for(i=0;i<pts;i++)
      {
        r.push_back(rn);
        rn+=step;
      }
    basis.zeros(pts,N);
    basis_i.zeros(pts,N);
    del_basis.zeros(pts,N);
    del_basis_i.zeros(pts,N);
  }
  lagrange (int Nl, double box,double step):N{Nl},a{box}
  {
    double rn;
    int i;
    for(i=0;i<N;i++)
      {
        x.push_back(0.);
        w.push_back(0.);
      }
    rn=step;
    while(rn<=a)
      {
        r.push_back(rn);
        rn+=step;
      }
    basis.zeros(r.size(),N);
    basis_i.zeros(r.size(),N);
    del_basis.zeros(r.size(),N);
    del_basis_i.zeros(r.size(),N);
  }
  void LagrangeBasis();
  void LagrangeBasis(double initial);
};
class MeanField : public nlpotential{
 public:
  string ws[5];
  double V;
  double VSO;
  double aV;
  double aSO;
  double RV;
  double RSO;
  double k;
  double rhc;
  vec Coulomb;
  cx_vec SpinOrbit; 
  cx_vec Nuclear;
  MeanField(int points,double radiusi):nlpotential(points,radiusi),V{0.},VSO{0.},aV{-1.},
                                       aSO{-1.},RV{-1.},RSO{-1.},k{-1.},rhc{-1.}
  {
    type="loc";
    Coulomb.zeros(points);
    SpinOrbit.zeros(points);
    Nuclear.zeros(points);
   }
 MeanField(int points,double radiusi,double Vi,
           double VSOi,double RVi,double RSOi,double aVi,
           double aSOi):nlpotential(points,radiusi),V{Vi},VSO{VSOi},aV{aVi},
                        aSO{aSOi},RV{RVi},RSO{RSOi},k{-1.},rhc{-1.}
  {
    int n;
	double delta_r;
    Coulomb.zeros(points);
    SpinOrbit.zeros(points);
    Nuclear.zeros(points);   
    delta_r=radius/double(points);
	for(n=0;n<points;n++)
      {
        r(n)=delta_r*(n+1.);
        Nuclear(n)=-V/(1.+exp((r(n)-RV)/aV));
      }
    pot=Nuclear;
  }
  MeanField(int points,double radiusi,double Vi,
            double VSOi,double RVi,double RSOi,double aVi,
            double aSOi,lagrange *lag):nlpotential(points,radiusi),V{Vi},VSO{VSOi},aV{aVi},
                                       aSO{aSOi},RV{RVi},RSO{RSOi},k{-1.},rhc{-1.}
  {
    int n;
	double delta_r;
	type="loc";
    delta_r=radius/double(points);
    r_reduced.zeros(lag->N);
    pot_reduced.zeros(lag->N);
	for(n=0;n<lag->N;n++)
      {
        r_reduced(n)=lag->x[n]*lag->a;
        pot_reduced(n)=-V/(1.+exp((r_reduced(n)-RV)/aV));
      }
  }
  void AddCoulomb(double q1q2);
  void AddSpinOrbit(int l,double j,double s);
  void Set();
};

class optical : public nlpotential {
 public:
  double V;
  double W;
  double Vso;
  double Wso;
  double Wd;
  double Wr;
  double radV;
  double radW;
  double radso;
  double radWd;
  double radcoul;
  double aV;
  double aW;
  double aso;
  double aWd;
  vec Coulomb;
  cx_vec SpinOrbit; 
  cx_vec Nuclear;
  optical():nlpotential()
  {
    type="loc";
  }
  optical(int points,double radiusi):nlpotential(points,radiusi)
  {
    type="loc";
    Coulomb.zeros(points);
    SpinOrbit.zeros(points);
    Nuclear.zeros(points);
  }
  optical(int points,double radiusi,int lagpoints):nlpotential(points,radiusi,lagpoints)
  {
    type="loc";
    Coulomb.zeros(points);
    SpinOrbit.zeros(points);
    Nuclear.zeros(points);
  }
  optical(int points,double radiusi,double Vi,double Wi,
          double Vsoi,double Wsoi,double Wdi,double radVi,double radWi,
          double radsoi,double radWdi,double radcouli,double aVi,
          double aWi,double asoi,double aWdi):nlpotential(points,radiusi),V{Vi},W{Wi},Vso{Vsoi},Wso{Wsoi},
                                     Wd{Wdi},radV{radVi},radW{radWi},radso{radsoi},radWd{radWdi},
                                     radcoul{radcouli},aV{aVi},aW{aWi},aso{asoi},aWd{aWdi}
  {
    int n;
	double delta_r;
	type="loc";
    delta_r=radius/double(points);
    //    cout<<"V: "<<V<<"  R: "<<radV<<"  a: "<<aV<<"\n";
	for(n=0;n<points;n++)
	{
      r(n)=delta_r*(n+1.);
      Nuclear(n)=-V/(1.+exp((r(n)-radV)/aV))-I*W/
        (1.+exp((r(n)-radW)/aW))-4.*I*Wd*
        exp((r(n)-radWd)/aWd)/((1.+exp((r(n)-radWd)/aWd))
                               *(1.+exp((r(n)-radWd)/aWd)));
	}
    pot=Nuclear;
  }
  optical(int points,double radiusi,double Vi,double Wi,
          double Vsoi,double Wsoi,double Wdi,double radVi,double radWi,
          double radsoi,double radWdi,double radcouli,double aVi,
          double aWi,double asoi,double aWdi,lagrange* lag):nlpotential(points,radiusi),V{Vi},W{Wi},Vso{Vsoi},Wso{Wsoi},
                                     Wd{Wdi},radV{radVi},radW{radWi},radso{radsoi},radWd{radWdi},
                                     radcoul{radcouli},aV{aVi},aW{aWi},aso{asoi},aWd{aWdi}
  {
    int n;
	double delta_r;
	type="loc";
    delta_r=radius/double(points);
    r_reduced.zeros(lag->N);
    pot_reduced.zeros(lag->N);
    for(n=0;n<points;n++)
	{
      r(n)=delta_r*(n+1.);
      pot(n)=-V/(1.+exp((r(n)-radV)/aV))-I*W/
        (1.+exp((r(n)-radW)/aW))-4.*I*Wd*
        exp((r(n)-radWd)/aWd)/((1.+exp((r(n)-radWd)/aWd))
                               *(1.+exp((r(n)-radWd)/aWd)));
	}
	for(n=0;n<lag->N;n++)
	{
      r_reduced(n)=lag->x[n]*lag->a;
      pot_reduced(n)=-V/(1.+exp((r_reduced(n)-radV)/aV))-I*W/
        (1.+exp((r_reduced(n)-radW)/aW))-4.*I*Wd*
        exp((r_reduced(n)-radWd)/aWd)/((1.+exp((r_reduced(n)-radWd)/aWd))
                               *(1.+exp((r_reduced(n)-radWd)/aWd)));
      // misc1<<lag->x[n]<<"  "<<r_reduced(n)<<"  "<<real(pot_reduced(n))<<"\n";
	}
    //exit(0);
  }
  void AddCoulomb(double q1q2);
  void AddSpinOrbit(int l,double j,double s);
  void Set();
  void Fill(int points,double radius);
  void Fill(lagrange *lag);
};





class state {
 public:
  int id;
  int nodes;
  double spec;
  double s;
  double radius;
  nlpotential* potential;
  nlpotential* potential_reduced;
  // reduced variables are defined at Gaussian integration points.
  // they are defined along with a lagrange object
  int l;
  double j;
  vec r;
  vec r_reduced;
  cx_vec wf;
  cx_vec wf_reduced;
  cx_vec vertex;
  cx_vec vertex_reduced;
  double energy;
  double D0;
  string file;
  state(){}
  state(int points):r(points),wf(points),vertex(points){}
  state(int points,lagrange *lag):r(points),wf(points),vertex(points),
                                  r_reduced(lag->N),wf_reduced(lag->N),vertex_reduced(lag->N){}
  state(int points,double radiusi):radius{radiusi},r(points),wf(points),vertex(points)
  {
    double step;
    int n;
    step=radiusi/double(points);
    for(n=0;n<points;n++)
      {
        r(n)=(n+1)*step;
      }
  }
  void Set(estado* st, double si){
    id=st->id;
    l=st->l;
    j=st->j;
    nodes=st->nodos;
    spec=st->spec;
    energy=st->energia;
    file=st->file;
    s=si;
  }

  complejo GenScatt(double mass,double q1q2, nlpotential* potentiali, lagrange* lag,cx_vec &c);
  complejo GenScatt(double q1q2,double mass,optical* pot);
  void GenBound(double radiusi, double mass, double q1q2, double energyi, MeanField* potentiali);
  void GenBound(double radiusi, double mass, double q1q2, MeanField* potentiali);
  void Normalize(int rule);
};
class parameters {
public:
  // Numerical parameters
  int points,lagpoints;
  double radius,box,lag_box,gamma,delta_e;
  int numpoles,mean_field,coupling;
  vector <optical> potentials;
};
