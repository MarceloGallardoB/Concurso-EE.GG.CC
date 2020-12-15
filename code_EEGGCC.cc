#include <iostream>      // inpuyts básicos
#include <fstream>        // input output file stream class
#include <cmath>          // librería matemática de base
#include <iomanip>        // manipuladores
#include <valarray>       // funcion valarray para simplificar operaciones
#include <vector>       // funciones vectoriales
#include "ConfigFile.h" // Archivo anexo que contiene el código que permite importar las simulaciones
                          

using namespace std; // para usar la libería de C++ estándar

//La clase engine es la que contiene el motor principal de las simulaciones, en este caso
//la clase que nos permite simlular el péndulo
class Engine
{

private:
  //definición de variables internas: propias al sistema dinámico
  
  double const pi=3.14159265358979323846264338327950288419716939937510582097494459230e0;
  double tfin=0.e0;		  // Tiempo final de la simulacion
  double g=0.e0;		  // aceleracion gravitacional 
  double mass=1.e0,q=0.e0,L=0.e0; //masa, carga eléctrica, longitud del pendulo
  double k=0.e0;		  // constante de fricción
  double omega=0.e0,phi=0.e0;	  // velocidad angular y desfase
  double E0=0.e0; 	  // intensidad del campo eléctrico
  valarray<double> winit=valarray<double>(0.e0,2); // posiciones y velocidades iniciales
  unsigned int sampling=1; // # de pasos cada cuanto se escribe en el output (t + N * samplig dt)
  unsigned int last;       
  ofstream *outputFile;    // puntero hacia el archivo de salida

  //funciones internas a la clase Engine 
   
  void printOut(bool write)
  
  {
    if((!write && last>=sampling) || (write && last!=1))
    { 
      *outputFile << t << " " << w[0] << " " << w[1] << " " \
      << energy() << " " << nonConservativePower() << " " << \
      (w[0]+pi)-(2.*pi)*floor((w[0]+pi)/(2.*pi))-pi << endl; // escritura en el output file
      last = 1;
    }
    else
    {
      last++;
    }
  }

 //Cálculo de la fuerzas no conservativas
  double nonConservativePower()
  {
	  double P = L*w[1]*(q*E0*sin(omega*t)*sin(w[0]) - L*k*w[1]);
    
    return P;

  }

 //energia mecánica del sistema (ver calculos analiticos)
  double energy() const
 {
	double E = 0.5*mass*pow(L,2)*pow(w[1],2) + mass*g*L*(1-cos(w[0])); 
   
   return E;
 }

  // iteración temporal: función redifinida en la subclase integrador
  virtual void step()=0;

protected:

  //Variables libres, se puede acceder desde otra sección del código que no sea la clase Engine

  double t=0.e0;  // tiempo discretizado
  double dt=1.e5; //paso de tiempo
  valarray<double> w=valarray<double>(0.0e0,2); //velocidades, posiciones
  valarray<double> constants=valarray<double>(0.0e0,3);  

  //Función accessible desde el integrador numérico (para este caso, Stormer Verlet)
  
  double electricAcceleration(double const& t_,double const& theta_)
  {
	  
	  double a_E= (q/(mass*L))*E0*sin(omega*t_ + phi)*sin(theta_);
   //se calcula la fuerza del campo eléctrico
    return a_E;
  }

 //acceleracion gravitacional
  double conservativeAcceleration(double const& theta_)
  {
	double a_P = -g/L*sin(theta_);
    // calculer l acceleration de gravite
    return a_P;
  }

//acceleracion por las fuerzas de fricción
  double dragAcceleration(double const& thetaDot_)
  {
	  double a_F = (-k/mass)*thetaDot_;
	 
    // calculer l acceleration de la trainee lineaire
    return a_F;
  }
  //acceleracion total
	double acceleration_totale (double const& t, double const&  theta, double const& thetaDot){
		
		return electricAcceleration(t, theta) + conservativeAcceleration(theta) + dragAcceleration(thetaDot);
		
	}
public:
  //Constructor de la clase: 
  Engine(ConfigFile configFile)
  {
    // variable locale
    tfin     = configFile.get<double>("tfin");		
    dt	     = tfin/configFile.get<double>("nsteps");		
    mass     = configFile.get<double>("mass");		
    q        = configFile.get<double>("q");		
    L        = configFile.get<double>("L");		
    g        = configFile.get<double>("g");	
    k        = configFile.get<double>("k");		
    E0       = configFile.get<double>("E0");	
    omega    = configFile.get<double>("omega");		
    phi      = configFile.get<double>("phi");		
    winit[0] = configFile.get<double>("theta0");
    winit[1] = configFile.get<double>("thetaDot0");	
    sampling = configFile.get<unsigned int>("sampling");

    // Apertura del archivo de salida
    outputFile = new ofstream(configFile.get<string>("output").c_str()); 
    outputFile->precision(15); //Los números se escriben con 15 decimales de precisión
  };

  // Destructor virtual
  virtual ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };

  // Simulación completa
  void run()
  {
    t = 0.e0; // t0
    w = winit;   // vector posicion, velocidad angular en t0
    last = 0; 
    printOut(true); 
    while(t<=tfin)
    {
      step();  // metodo step: integracion 
      printOut(false); 
    }
    printOut(true); 
  };

};

// Integrador Stormer Verlet:
class EngineStormerVerlet: public Engine
{
protected:
public:

  // Constructor
  EngineStormerVerlet(ConfigFile configFile): Engine(configFile){}
  void step()
  {
    double xmem = w[0] ;
    w[0]=w[0]+w[1]*dt + 0.5*acceleration_totale(t,w[0],w[1])*dt*dt;
    double v = w[1] + 0.5*acceleration_totale(t,w[0],w[1])*dt;
    w[1] = w[1]+0.5*(acceleration_totale(t+dt,w[0],v) +  acceleration_totale(t,xmem,v))*dt;
    t+=dt;  
  }
};

int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // 
  if(argc>1)
    inputPath = argv[1];

  ConfigFile configFile(inputPath); 

  for(int i(2); i<argc; ++i) 
    configFile.process(argv[i]);
  string schema(configFile.get<string>("schema"));

  Engine* engine; 

  if(schema == "StormerVerlet" || schema == "SV")
  {

    engine = new EngineStormerVerlet(configFile);
  }
  else
  {
    cerr << "Schema inconnu" << endl;
    return -1;
  }

  engine->run(); 

  delete engine; 
  cout << "Fin de la simulation." << endl;
  return 0;
}
