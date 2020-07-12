#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

std::ofstream fout ("Ergebnisse.txt");

// Testfunktionen als Funktor
class Pol1 {
public:
  double operator()(double x) { return 3 * x + 2; }
};

class Pol2 {
public:
  double operator()(double x) { return -2*x*x+3*x+1; }
};

class Gauss {
public:
  double operator()(double x) { return 1 / (sqrt(M_PI * 2)) * exp(-x * x / 2); }
};

// berechnet Werte nach Trapezformel von I_0 bis I_N
template < class Functor > std :: vector < double >
trapez ( Functor f, double a, double b, int N) {
  //Pol1 f;
  std::vector<double> I(N + 1); // Feld mit N+1 Eintraegen
  const double h = b - a;
  I[0] = h / 2 * (f(a) + f(b));
  for (int k = 1; k <= N; ++k) {
    double sum = 0;
    int n = pow(2, k);
    for (int i=1; i<n; ++i){
      sum += f(a+i*h/n);
    }
    I[k] = h/n*(0.5*f(a)+0.5*f(b)+sum); // setze k-ten Wert im Feld
  }
  return I;
}

// berechnet die Richardsonextrapolation aus I(k-1)  und I(k)
double richardson(double Iprev, double I) { return (4*I-Iprev)/3; }

// berechet Naeherungen ueber das Romberg-Verfahren
// I: Ergebnis von trapez()
std::vector<std::vector<double>> romberg(std::vector<double> I) {
  const int N = I.size() - 1;
  std::vector<std::vector<double>> R(N + 1);

  for (int k = 0; k <= N; ++k) {
    R[k].push_back(I[k]);
  }
 
  for (int j = 1; j <= N; ++j) {
    for(int k = 0; k<=N-j; ++k) {
      double extr = R[k+1][j-1]+((R[k+1][j-1]-R[k][j-1])/(pow(2,2*j)-1));
      R[k].push_back(extr);
     }   
  }
  
  return R;
}


void testeAufgabe1() {
  Pol1 f;
  std::vector<double> I_f = trapez(f, 0, 3, 3);
  for (double tf : I_f) {
    std::cout << "A1: f:" << tf << " == " << 19.5 << ":" << (tf == 19.5 ? "ja" : "nein") << std::endl;
  }
  Pol2 g;
  std::cout << "A1: g(1) = 2 ?" << (g(1) == 2 ? "ja" : "nein") << std::endl;
  std::vector<double> I_g = trapez(g, 0, 3, 3);
  std::cout << "A1: g0:" << I_g[0] << " == " << -10.5 << ":" << (I_g[0] == -10.5 ? "ja" : "nein") << std::endl;
  std::cout << "A1: g1:" << I_g[1] << " == " << -3.75 << ":" << (I_g[1] == -3.75 ? "ja" : "nein") << std::endl;
  std::cout << "A1: g2:" << I_g[2] << " == " << -2.0625 << ":" << (I_g[2] == -2.0625 ? "ja" : "nein") << std::endl;
  double rich = richardson(I_g[0], I_g[1]);
  std::cout << "A1: Richardson : " << rich << " : " << (rich == -1.5 ? "ja " : "nein") << std::endl;
}

void testeAufgabe2() {
  Pol1 f;
  std::vector<std::vector<double>> Rf = romberg(trapez(f, 0, 3, 3));
  bool alle_richtig = true;
  int entries = 0;
  for (auto row : Rf) {
    for (double val : row) {
      alle_richtig &= val == 19.5;
      ++entries;
    }
  }
  std::cout << "A2: alle Eintraege für f sind 1.5:" << (alle_richtig ? "ja" : "nein") << std::endl;
  std::cout << "A2: korrekte Zahl an Einträgen:" << (entries == 10 ? "ja" : "nein") << std::endl;
  Pol2 g;
  std::vector<std::vector<double>> Rg = romberg(trapez(g, 0, 3, 3));
  std::cout << "A2: R[1][1] und R[2][1] für g gleich -1.5: " << ((Rg[1][1] == -1.5) && (Rg[2][1] == -1.5) ? " ja " : " nein") << std::endl;
}


int main() {

  fout<<"Grenzen der Integrale a=0, b=3"<<std::endl;

  fout<<"f(x)          || "<<" I_0(f)  | "<<" I_1(f)  | "<<" I_2(f)  | "<<" I_3(f)  | "<<"Integral(f(x))"<<std::endl;

  Pol1 f1;
  std::vector<double> tf1 = trapez(f1, 0., 3., 3);
  fout<<"3x+2          ||   "<<tf1[0]<<"   |   "<<tf1[1]<<"   |   "<<tf1[2]<<"   |   "<<tf1[3]<<"   |   "<<"19,5"<<std::endl;

  Pol2 f2;
  std::vector<double> tf2 = trapez(f2, 0., 3., 3);
  fout<<"-2*x^2+3x+1   ||  "<<tf2[0]<<"   |  "<<tf2[1]<<"   | "<<tf2[2]<<"  | "<<tf2[3]<<" |  "<<" -1,5"<<std::endl;

  Gauss f3;
  std::vector<double> tf3 = trapez(f3, 0., 3., 3);
  fout<<"Gauss Fnkt.   || "<<tf3[0]<<" | "<<tf3[1]<<" | "<<tf3[2]<<" | "<<tf3[3]<<" |   "<<"0,49865"<<std::endl;

  fout<<std::endl;

  fout<<"f(x)          || "<<"Rich(I_1,I_0) | "<<"Rich(I_2,I_1) | "<<"Rich(I_3,I_2) | "<<"Integral(f(x))"<<std::endl;
  fout<<"3x+2          ||     "<<richardson(tf1[0],tf1[1])<<"      |     "<<richardson(tf1[1],tf1[2])<<"      |     "<<richardson(tf1[2],tf1[3])<<"      |   "<<"19,5"<<std::endl;
  fout<<"-2*x^2+3x+1   ||     "<<richardson(tf2[0],tf2[1])<<"      |     "<<richardson(tf2[1],tf2[2])<<"      |     "<<richardson(tf2[2],tf2[3])<<"      |   "<<" -1,5"<<std::endl;
  fout<<"Gauss Fnkt.   ||   "<<richardson(tf3[0],tf3[1])<<"    |   "<<richardson(tf3[1],tf3[2])<<"    |   "<<richardson(tf3[2],tf3[3])<<"    |   "<<"0,49865"<<std::endl;



  testeAufgabe1();
  testeAufgabe2();

}
