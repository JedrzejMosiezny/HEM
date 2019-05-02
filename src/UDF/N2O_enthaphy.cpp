#include <cstdio>
#include <cmath>
#include <iostream>
using namespace std;

double pc, Tc;

double h_s(double T)
{
	double b1, b2, b3, b4, b5, x;
	b1 = -200.0;
	b2 = 440.055;
	b3 = -459.701;
	b4 = 434.081;
	b5 = -485.338;
	
	x = 1-T/Tc;
	
	
	return (b1+b2*pow(x,1./3)+b3*pow(x,2./3)+b4*x+b5*pow(x,4./3))*1000;  //[J/kg]
}

double p_s(double T)
{
	double b1, b2, b3, b4, x;
	b1 = -6.71893;
	b2 = 1.35966;
	b3 = -1.3779;
	b4 = -4.051;
	
	x = 1-T/Tc;
	
	return exp(Tc/T*(b1*x+b2*pow(x,3./2)+b3*pow(x,5./2)+b4*pow(x,5)))*pc;
}

double cp_s(double T)
{
	double b1, b2, b3, b4, b5, x;
	b1 = 132.632;
	b2 = 0.052187;
	b3 = -0.364923;
	b4 = -1.20233;
	b5 = 0.536141;
	
	x = 1-T/Tc;
	
	return b1*(1+b2*pow(x,-2./3)+b3*pow(x,-1./3)+b4*pow(x,1./3)+b5*pow(x,2./3))*1000;  //[J/kgK]
}

double T_s(double p) ///dobrze by bylo by jakos odroznic nazwe funkcji od parametru, i zeby bylo wiadomo co bierze za argument bo zazwyczaj jest temp a tu na odwrot
{
	double a, b, c, d, x;
	a = 0.359857348;
	b = 5.960225154;
	c = 48.17054333;
	d = 309.6838723;
	
	x = log(p/pc);
	
	return a*pow(x,3)+b*pow(x,2)+c*x+d;
}

double dpdT(double T)
{
	double deltaT=1e-3;
	
	return (p_s(T+deltaT)-p_s(T-deltaT))/(2*deltaT);
	
}

double dhs_dTs(double T) //analitycznie wyznaczona pochodna hs wzd³u¿ linii saturacji, ka¿dy wspó³czynnik uzyskaæ mo¿na ze wspo³czynników oryginalnego wielomianu hs(Tr) jako  Bn = -b1*bn*potegaPrzy(1-Tr)
{
	double B2, B3, B4, B5, x;
	B2 = -146.685;
	B3 = 306.4673333;
	B4 = -434.081;
	B5 = 647.1173333;
	
	x = 1-T/Tc;
	
	
	return (B2*pow(x,-2./3)+B3*pow(x,-1./3)+B4+B5*pow(x,1./3))*1000/Tc;  //[J/kg/K]
	
}

double h(double p, double T)
{
	double Tsat_p, hT, hp, delta_hp, delta_hT,hs_p,hs_T;
	
	Tsat_p = T_s(p);// znajdz Ts dla zadanego cisnienia
	delta_hp = cp_s(Tsat_p)*(T-Tsat_p); //izobaryczna korekta entalpii
	delta_hT = (dhs_dTs(T)-cp_s(T))*(p-p_s(T))/dpdT(T);// wczesniej bylo delta_hT = cp_s(T)*pow(dpdT(T),-1)*(p-pc);
	hp = h_s(Tsat_p)+delta_hp;
	hT = h_s(T)+delta_hT;
	cout << "delta hp = " << delta_hp << endl;
	cout << "delta hT = " << delta_hT << endl;
	return (hp*delta_hT*delta_hT+hT*delta_hp*delta_hp)/(delta_hp*delta_hp+delta_hT*delta_hT);//zamienilem kwadraty na mnozenia bo to pewnie mniej wymagajkace obliczen niz pow
}

int main()
{
	pc = 73e5; //[Pa]
	Tc = 309; // [K]
	
	double p_input, T_input;
	
	p_input=30e5;
	T_input=273.15;
		
	cout << "Tcrit = " << pc << endl;
	cout << "pcrit = " << Tc << endl<<endl;
	
	cout << "p = " << p_input << endl;
	cout << "T = " << T_input << endl<<endl;
	double hsat_p = h_s(T_s(p_input));
	double p_sat_1 = p_s(T_input);
	cout << "h_s = " << h_s(T_input) << endl;
	cout << "p_s = " << p_sat_1 << endl;
	cout << "cp_s = " << cp_s(T_input) << endl;
	cout << "T_s(p) = " << T_s(p_input) << endl;
	cout << "dpdT = " << dpdT(T_input) << endl;
	cout << "sqrt4 = " << pow(4, 1./2) << endl;
	cout << "dhs/dT = " << dhs_dTs(T_input) << endl;
	cout << "1-Tr = " << 1-T_input/Tc << endl;
	double hgiven = h(p_input, T_input);
	cout << "h (p,T) = " << hgiven << endl;
	cout << "average isobaric Cp from saturation pressure = " << (hgiven - hsat_p)/(T_input-T_s(p_input)) << endl;
    return 0;
    
}
