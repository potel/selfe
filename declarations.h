void TestBuild(parameters *parm);
void Error(const string message);
void PolarPotential(double energy,lagrange *lag,nlpotential *v);
void Expand(state *st,lagrange *lag,cx_vec c);
void Write(vec r,cx_vec f,ofstream *file);
void Write(vec r,vec f,ofstream *file);
void Write(vec r,cx_mat f,ofstream *file);
double Absorption(state *st,nlpotential *v,double mass,lagrange *lag);
void Converge(nlpotential *v,int iter,lagrange *lag);
void ReadParameters(const char *fname,parameters *x);
int ReadPotentials(char *s,const char key[100],vector <optical> &pot,
                   ifstream *fp,int poins,int lagpoints,double radius);
void ReadParD(char *s,const char key[20], int *par);
void ReadParF(char *s,const char key[20],double *par);
void ReadParS(char *s,const char key[20], char *par);
void ReadParStr(string s,const char key[20], string par);
