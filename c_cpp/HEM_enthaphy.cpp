#include <cstdio>
#include <cmath>
#include <iostream>
using namespace std;

double pc, Tc;

double h_s(double T)
{
	double b1, b2, b3, b4, b5, x;
	b1 = -200;
	b2 = 440.055;
	b3 = -459.701;
	b4 = 434.081;
	b5 = -485.338;
	
	x = 1-T/Tc;
	
	
	return b1*(1+b2*pow(x,1/3)+b3*pow(x,2/3)+b4*x+b5*pow(x,4/3))*1000;  //[J/kg]
}

double p_s(double T)
{
	double b1, b2, b3, b4, x;
	b1 = -6.71893;
	b2 = 1.35966;
	b3 = -1.3779;
	b4 = -4.051;
	
	x = 1-T/Tc;
	
	return exp(Tc/T*(b1*x+b2*pow(x,3/2)+b3*pow(x,5/2)+b4*pow(x,5)))*pc;
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
	
	return b1*(1+b2*pow(x,-2/3)+b3*pow(x,-1/3)+b4*pow(x,1/3)+b5*pow(x,2/3))*1000;  //[J/kgK]
}

double T_s(double p)
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

double h(double p, double T)
{
	double hT, hp, delta_hp, delta_hT;
	
	hp = h_s(T_s(p))+cp_s(T)*(T-Tc);
	hT = h_s(T);
	delta_hp = cp_s(T)*(T-Tc);
	delta_hT = cp_s(T)*pow(dpdT(T),-1)*(p-pc);
	
	return hp*pow(delta_hT,2)/(pow(delta_hT,2)+pow(delta_hp,2))+hT*pow(delta_hp,2)/(pow(delta_hT,2)+pow(delta_hp,2));
}

int main()
{
	pc = 73e5; //[Pa]
	Tc = 309; // [K]
	
	double p_input, T_input;
	
	p_input=70e5;
	T_input=300;
		
	cout << "Tsat = " << pc << endl;
	cout << "psat = " << Tc << endl<<endl;
	
	cout << "p = " << p_input << endl;
	cout << "T = " << T_input << endl<<endl;
	
	cout << "h_s = " << h_s(T_input) << endl;
	cout << "p_s = " << p_s(T_input) << endl;
	cout << "cp_s = " << cp_s(T_input) << endl;
	cout << "T_s(p) = " << T_s(p_input) << endl;
	cout << "dpdT = " << dpdT(T_input) << endl;
	cout << "h = " << h(p_input, T_input) << endl;

    return 0;
    
}
