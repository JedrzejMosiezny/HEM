
 #include "udf.h"
 #include "stdio.h"
 #include "ctype.h"
 #include "stdarg.h"
 
 #define MW 44   

 
#define Tc_co2 304.13 
#define pc_co2 467.44 
#define rhoc_co2 7377300 

#define cp_co2_b1 20.01502241
#define cp_co2_b2 -0.064337181
#define cp_co2_b3 59.00656609
#define cp_co2_b4 124.6784428
#define cp_co2_b5 -550.5406959

#define R_co2_b1 0.287858279
#define R_co2_b2 -1.054857595
#define R_co2_b3 944.5124684
#define R_co2_b4 -1575.56863
#define R_co2_b5 -1120.486766

/* zakres temperatur dla cp */
#define Tmincp 150
#define Tmaxcp 288
 
 /* zakres temperatur dla R */
#define TminR 160
#define TmaxR 184.75
 
 
 static int (*usersMessage)(char *,...);
 static void (*usersError)(char *,...);
 
 DEFINE_ON_DEMAND(I_do_nothing)
 {
    /* This is a dummy function to allow us to use */
    /* the Compiled UDFs utility      */
 }
 
 
 void IDEAL_error(int err, char *f, char *msg)
 {
  if (err)
    usersError("IDEAL_error (%d) from function: %s\n%s\n",err,f,msg);
 }
 
 void IDEAL_Setup(Domain *domain, cxboolean vapor_phase, char *filename,
     int (*messagefunc)(char *format, ...),
     void (*errorfunc)(char *format, ...))
 {  
      usersMessage = messagefunc;
      usersError  = errorfunc;
      usersMessage("\nLoading Real-Ideal Library: %s\n", filename);
 }
 
 double IDEAL_density(cell_t cell, Thread *thread, cxboolean vapor_phase, double Temp, double press, double yi[])
 {
      double x=1-Temp/Tc_co2;
	  
	  double xminR=1-TminR/Tc_co2;

      double RGAS;
      
      if (Temp<TminR) RGAS = R_co2_b1*(1+R_co2_b2/xminR+R_co2_b3*xminR+R_co2_b4*xminR*xminR+R_co2_b5*xminR*xminR*xminR);
      else RGAS = R_co2_b1*(1+R_co2_b2/x+R_co2_b3*x+R_co2_b4*x*x+R_co2_b5*x*x*x);
      
      double r = press/(RGAS*Temp); /* Density at Temp & press */
      
      return r;      /* (Kg/m^3) */
 }
 
 double IDEAL_specific_heat(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      double x=1-Temp/Tc_co2;
	  double xmincp=Tmincp/Tc_co2;
	  double xmaxcp=Tmaxcp/Tc_co2;
      double cp;
      
      if (Temp<Tmincp) {cp = cp_co2_b1*(1+cp_co2_b2/xmincp+cp_co2_b3*xmincp+cp_co2_b4*xmincp*xmincp+cp_co2_b5*xmincp*xmincp*xmincp);}
      if (Temp>Tmaxcp) {cp = cp_co2_b1*(1+cp_co2_b2/xmaxcp+cp_co2_b3*xmaxcp+cp_co2_b4*xmaxcp*xmaxcp+cp_co2_b5*xmaxcp*xmaxcp*xmaxcp);}
      else {cp = cp_co2_b1*(1+cp_co2_b2/x+cp_co2_b3*x+cp_co2_b4*x*x+cp_co2_b5*x*x*x);}
      
      return cp;     /* (J/Kg/K) */
 }
 
 double IDEAL_enthalpy(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      double x=1-Temp/Tc_co2;
    
      double Tref = 293.15;
      
      double xr=1-Tref/Tc_co2;
      
	  double hs = - 7034295.45; //standard heat of formation
      
      double h = -Tc_co2*cp_co2_b1*((x+cp_co2_b2*log(x)+cp_co2_b3*0.5*x*x+1./3*cp_co2_b4*x*x*x+cp_co2_b5*0.25*pow(x,4))
      -(xr+cp_co2_b2*log(xr)+cp_co2_b3*0.5*xr*xr+1./3*cp_co2_b4*xr*xr*xr+cp_co2_b5*0.25*pow(xr,4))) + hs;

      return h; 
 }
 
 #define TDatum 288.15
 #define PDatum 1.01325e5 
 

 double IDEAL_entropy(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
     
      double x=1-Temp/Tc_co2;
	  
	  double xminR=1-TminR/Tc_co2;

      double RGAS;
      
      if (Temp<TminR) RGAS = R_co2_b1*(1+R_co2_b2/xminR+R_co2_b3*xminR+R_co2_b4*xminR*xminR+R_co2_b5*xminR*xminR*xminR);
      else RGAS = R_co2_b1*(1+R_co2_b2/x+R_co2_b3*x+R_co2_b4*x*x+R_co2_b5*x*x*x);
      
	  double xmincp=Tmincp/Tc_co2;
	  
	  double xmaxcp=Tmaxcp/Tc_co2;
	  
      double cp;
      
      if (Temp<Tmincp) {cp = cp_co2_b1*(1+cp_co2_b2/xmincp+cp_co2_b3*xmincp+cp_co2_b4*xmincp*xmincp+cp_co2_b5*xmincp*xmincp*xmincp);}
      if (Temp>Tmaxcp) {cp = cp_co2_b1*(1+cp_co2_b2/xmaxcp+cp_co2_b3*xmaxcp+cp_co2_b4*xmaxcp*xmaxcp+cp_co2_b5*xmaxcp*xmaxcp*xmaxcp);}
      else {cp = cp_co2_b1*(1+cp_co2_b2/x+cp_co2_b3*x+cp_co2_b4*x*x+cp_co2_b5*x*x*x);}
       
     
      double s=cp*log(fabs(Temp/TDatum))+RGAS*log(fabs(PDatum/P));
      
      return s;      /* (J/Kg/K) */
 }
 
 
 double IDEAL_mw(double yi[])
 {
      return MW;     /* (Kg/Kmol) */
 }
 
 double IDEAL_speed_of_sound(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      double x=1-Temp/Tc_co2;
	  
	  double xminR=1-TminR/Tc_co2;

      double RGAS;
      
      if (Temp<TminR) RGAS = R_co2_b1*(1+R_co2_b2/xminR+R_co2_b3*xminR+R_co2_b4*xminR*xminR+R_co2_b5*xminR*xminR*xminR);
      else RGAS = R_co2_b1*(1+R_co2_b2/x+R_co2_b3*x+R_co2_b4*x*x+R_co2_b5*x*x*x);
      
	  double xmincp=Tmincp/Tc_co2;
	  
	  double xmaxcp=Tmaxcp/Tc_co2;
	  
      double cp;
      
      if (Temp<Tmincp) {cp = cp_co2_b1*(1+cp_co2_b2/xmincp+cp_co2_b3*xmincp+cp_co2_b4*xmincp*xmincp+cp_co2_b5*xmincp*xmincp*xmincp);}
      if (Temp>Tmaxcp) {cp = cp_co2_b1*(1+cp_co2_b2/xmaxcp+cp_co2_b3*xmaxcp+cp_co2_b4*xmaxcp*xmaxcp+cp_co2_b5*xmaxcp*xmaxcp*xmaxcp);}
      else {cp = cp_co2_b1*(1+cp_co2_b2/x+cp_co2_b3*x+cp_co2_b4*x*x+cp_co2_b5*x*x*x);}
       
      
      return sqrt(Temp*cp*RGAS/(cp-RGAS)); /* m/s */
 }
 
 double IDEAL_viscosity(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      double mu=1.7894e-05;
      return mu;     /* (Kg/m/s) */
 }
 
 double IDEAL_thermal_conductivity(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      double ktc=0.0242;
      return ktc;      /* W/m/K */
 } 
 double IDEAL_rho_t(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      /* derivative of rho wrt. Temp at constant p */
      double rho_t=-density/Temp;
      return rho_t;     /* (Kg/m^3/K) */
 } 
 double IDEAL_rho_p(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      /* derivative of rho wrt. pressure at constant T */
      
      double x=1-Temp/Tc_co2;
	  
	  double xminR=1-TminR/Tc_co2;

      double RGAS;
      
      if (Temp<TminR) RGAS = R_co2_b1*(1+R_co2_b2/xminR+R_co2_b3*xminR+R_co2_b4*xminR*xminR+R_co2_b5*xminR*xminR*xminR);
      else RGAS = R_co2_b1*(1+R_co2_b2/x+R_co2_b3*x+R_co2_b4*x*x+R_co2_b5*x*x*x);
 
      
      double rho_p=1.0/(RGAS*Temp);
      return rho_p;    /* (Kg/m^3/Pa) */
 }
 
 double IDEAL_enthalpy_t(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      /* derivative of enthalpy wrt. Temp at constant p */
      
      double x=1-Temp/Tc_co2;
      
	  double xmincp=Tmincp/Tc_co2;
	  
	  double xmaxcp=Tmaxcp/Tc_co2;
	  
      double cp;
      
      if (Temp<Tmincp) {cp = cp_co2_b1*(1+cp_co2_b2/xmincp+cp_co2_b3*xmincp+cp_co2_b4*xmincp*xmincp+cp_co2_b5*xmincp*xmincp*xmincp);}
      if (Temp>Tmaxcp) {cp = cp_co2_b1*(1+cp_co2_b2/xmaxcp+cp_co2_b3*xmaxcp+cp_co2_b4*xmaxcp*xmaxcp+cp_co2_b5*xmaxcp*xmaxcp*xmaxcp);}
      else {cp = cp_co2_b1*(1+cp_co2_b2/x+cp_co2_b3*x+cp_co2_b4*x*x+cp_co2_b5*x*x*x);}
       
      
      return cp;
 }
 
 double IDEAL_enthalpy_p(cell_t cell, Thread *thread,double Temp, double density, double P, double yi[])
 {
      /* derivative of enthalpy wrt. pressure at constant T  */
      /* general form dh/dp|T = (1/rho)*[ 1 + (T/rho)*drho/dT|p] */
      /* but for ideal gas dh/dp = 0        */
      return 0.0 ;
 }
 
 UDF_EXPORT RGAS_Functions RealGasFunctionList =
 {
      IDEAL_Setup,                   /* initialize */
      IDEAL_density,                 /* density */
      IDEAL_enthalpy,                /* enthalpy */
      IDEAL_entropy,                 /* entropy */
      IDEAL_specific_heat,           /* specific_heat */
      IDEAL_mw,                      /* molecular_weight */
      IDEAL_speed_of_sound,          /* speed_of_sound */
      IDEAL_viscosity,               /* viscosity */
      IDEAL_thermal_conductivity,    /* thermal_conductivity */
      IDEAL_rho_t,                   /* drho/dT |const p */  
      IDEAL_rho_p,                   /* drho/dp |const T */
      IDEAL_enthalpy_t,              /* dh/dT |const p  */
      IDEAL_enthalpy_p        /* dh/dp |const T  */
 };
 /**************************************************************/