#include<iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#define ROZMIAR 30
#define ROZMIAR_PLANSZY ROZMIAR+6
#define ILOSC_ROWNAN ROZMIAR*ROZMIAR*2
#define TEST 0 
//test=0 brak zbêdnych wiadomoœci
//test=1 du¿ych plansze
//test=2 ma³e plansze, 
const double eps = 1e-8;

#define MAXN ILOSC_ROWNAN

double err = eps;



bool converge(double *xk, double *xkp,int ile_mam_rownan)
{
	double norm = 0;
	for (int i = 0; i<ile_mam_rownan; i++)
	{
		norm += (xk[i] - xkp[i])*(xk[i] - xkp[i]);
	}
	if (sqrt(norm) >= eps)
		return false;
	err = eps;
	return true;
}


bool gauss(int n, double ** AB, double * X)
{

	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}


struct rownanie2 {
	int pozycja_gracza_1;
	int pozycja_gracza_2;
	int kogo_ruch;

	int tab_row[6];
	int ca³osc;//ile czesci wygrywajacych/przegrywajacych
};



double monte_carlo (int ilosc_pol, int tablica[ROZMIAR],int ilosc_gier) {
	double wynik = 0;
	int pozycja_1 = 0;
	int pozycja_2 = 0;
	int licznik_wygranych = 0;
	//int ilosc_gier = 10000000;
	int wygrany = 0;//0 gra trwa, 1 gracz 1 wygra³, -1 gracz 1 przegra³
	for (int i = 0; i < ilosc_gier; i++)
	{
		pozycja_1 = 0;
		pozycja_2 = 0;
		wygrany = 0;
		while (wygrany == 0) {
			pozycja_1 = pozycja_1 + (rand() % 6 + 1);
			if (pozycja_1 >= ilosc_pol)
			{
				wygrany = 1;
				break;
			}
			pozycja_1 += tablica[pozycja_1];
			pozycja_2 = pozycja_2 + (rand() % 6 + 1);
			if (pozycja_2 >= ilosc_pol)
			{
				wygrany = -1;
				break;
			}
			pozycja_2 += tablica[pozycja_2];
		}
		if (wygrany == 1)licznik_wygranych++;
	}
	return (double)(licznik_wygranych / (double)ilosc_gier);

}
int zainicjuj_opis(rownanie2 tabelka_rownan[ILOSC_ROWNAN], int plansza[ROZMIAR_PLANSZY],int ile_opisow, int A,int B, int gracz) {
	for (int i = 0; i < ile_opisow;i++)
	{
		if (tabelka_rownan[i].kogo_ruch == gracz &&
			tabelka_rownan[i].pozycja_gracza_1 == A &&
			tabelka_rownan[i].pozycja_gracza_2 == B) {
			//cout << "nie inicjuje bo juz jest" << endl;
			return i;
		}
	}
	
	tabelka_rownan[ile_opisow].kogo_ruch = gracz;
	tabelka_rownan[ile_opisow].pozycja_gracza_1 = A;
	tabelka_rownan[ile_opisow].pozycja_gracza_2 = B;
	//cout << "zainicjowano (" << A <<","<< B <<","<< gracz << ") jako numer " << ile_opisow << endl;
	
	return ile_opisow;


}
	


int	uzupelniaj_tabelke(rownanie2 tabelka_rownan[ILOSC_ROWNAN], int plansza[ROZMIAR_PLANSZY],int ile_mam_pol) {
	int ile_opisow = 1;
	int ile_uzupelnionych_rownan = 0;

	int pole_A;
	int pole_B;
	int pole_A_prim;
	int pole_B_prim;
	int gracz;
	int numer_opisu;
	while (ile_opisow != ile_uzupelnionych_rownan)
	{
		pole_A = tabelka_rownan[ile_uzupelnionych_rownan].pozycja_gracza_1;
		pole_B = tabelka_rownan[ile_uzupelnionych_rownan ].pozycja_gracza_2;
		gracz = tabelka_rownan[ile_uzupelnionych_rownan ].kogo_ruch;
		tabelka_rownan[ile_uzupelnionych_rownan ].ca³osc = 0;
		for (int i = 1; i <= 6; i++)
		{
			if (gracz == 1) {
				pole_A_prim = pole_A + i + plansza[pole_A + i];
				pole_B_prim = pole_B;
				if (pole_A_prim >= ile_mam_pol)
				{
					tabelka_rownan[ile_uzupelnionych_rownan ].tab_row[i - 1] = -1;
					tabelka_rownan[ile_uzupelnionych_rownan ].ca³osc++;
				}
				else
				{
					numer_opisu = zainicjuj_opis(tabelka_rownan, plansza, ile_opisow, pole_A_prim, pole_B_prim, 3 - gracz);
					if (numer_opisu == ile_opisow)ile_opisow++;
					tabelka_rownan[ile_uzupelnionych_rownan ].tab_row[i - 1] = numer_opisu;
				}

			
			}
			else //gracz==2
			{
				pole_B_prim = pole_B + i + plansza[pole_B + i];
				pole_A_prim = pole_A;
				if (pole_B_prim >= ile_mam_pol)
				{
					tabelka_rownan[ile_uzupelnionych_rownan ].tab_row[i - 1] = -1;
					tabelka_rownan[ile_uzupelnionych_rownan].ca³osc = 0;
				}
				else
				{
					numer_opisu = zainicjuj_opis(tabelka_rownan, plansza, ile_opisow, pole_A_prim, pole_B_prim, 3 - gracz);
					if (numer_opisu == ile_opisow)ile_opisow++;
					tabelka_rownan[ile_uzupelnionych_rownan ].tab_row[i - 1] = numer_opisu;
				}



			}
			
			//cout << "uzupelniono fragment rownania "<<ile_uzupelnionych_rownan  << endl;
			
			
		}
		if (TEST == 1 || TEST==2)
		{
			cout << "rownanie numer:  " << ile_uzupelnionych_rownan << ". (" << pole_A << "," << pole_B << "," << gracz << ")  ";
			for (int i = 1; i <= 6; i++)
				cout << "\t  " << tabelka_rownan[ile_uzupelnionych_rownan].tab_row[i - 1];
			cout << "\t cALOSC " << tabelka_rownan[ile_uzupelnionych_rownan].ca³osc;
			cout << endl;
		}
	ile_uzupelnionych_rownan++;
	}


	return ile_uzupelnionych_rownan ;
}

void tworz_maciez_AB(rownanie2 tabelka_rownan[ILOSC_ROWNAN],int  ile_mam_rownan, double **tab2) {
	for (int i = 0; i < ile_mam_rownan; i++)
	{
		tab2[i][i] = -6;
		
	}
	for (int i = 0; i < ile_mam_rownan; i++)
	{
		tab2[i][ile_mam_rownan] = -tabelka_rownan[i].ca³osc;
		//b[i][0]= -tabelka_rownan[i].ca³osc;
	}
	for (int i = 0; i < ile_mam_rownan; ++i)
	{
		for (int j = 0; j < 6; ++j)
			if (tabelka_rownan[i].tab_row[j] != -1)
			{
				tab2[i][tabelka_rownan[i].tab_row[j]]++;

			}
		
	}
	//for (int i = 0; i < ile_mam_rownan; i++)
	//	for (int j = 0; j < ile_mam_rownan; j++)
	//		A[i][j] = tab2[i][j];

}


int main() {
	srand(time(NULL));

	int ilosc_pol;
	char znak;
	int plansza[ROZMIAR_PLANSZY] = {0};
	int pole;
	int cofacz;
	cin >> ilosc_pol;
	while (getchar() != '.') {
		getchar();
		cin >> pole;
		getchar();
		cin >> cofacz;
		getchar();
		plansza[pole] = cofacz;
	}
	if(TEST==1 || TEST==2)
	for (int i = 0; i < ilosc_pol; i++)
	{
		cout << "pole nr " << i << " cofa o " << plansza[i] << "\t"<< (rand() % 6 + 1)<<endl;
	}

	rownanie2 tabelka_rownan[ILOSC_ROWNAN];
	tabelka_rownan[0].pozycja_gracza_1 = 0;
	tabelka_rownan[0].pozycja_gracza_2 = 0;
	tabelka_rownan[0].kogo_ruch= 1;

	int ile_mam_rownan;
	ile_mam_rownan=uzupelniaj_tabelke(tabelka_rownan, plansza,ilosc_pol);
	cout << "wygenerowalam " << ile_mam_rownan << " rownan" << endl;

	int w = ile_mam_rownan;
	int k = ile_mam_rownan + 1;
	double **tabAB = new double *[w]; //alokacja pamieci
	for (int i = 0; i < w; ++i)
	{
		tabAB[i] = new double[k]; //alokacja pamieci
		for (int j = 0; j < k; ++j)
			tabAB[i][j] = 0;
	}
	tworz_maciez_AB(tabelka_rownan, ile_mam_rownan, tabAB);

	//////////////////////////
	
	double **A = new double *[ile_mam_rownan]; //alokacja pamieci
	for (int i = 0; i < ile_mam_rownan; ++i)
	{
		A[i] = new double[ile_mam_rownan]; //alokacja pamieci
		for (int j = 0; j < ile_mam_rownan; ++j)
			A[i][j] = tabAB[i][j];
	}

	double **b = new double *[ile_mam_rownan]; //alokacja pamieci
	for (int i = 0; i < ile_mam_rownan; ++i)
	{
		b[i] = new double[ile_mam_rownan]; //alokacja pamieci
		for (int j = 0; j < ile_mam_rownan; ++j)
		{
			b[i][j] = 0;
			if (j == 0)
				b[i][0]  = tabAB[i][ile_mam_rownan];
		}
		
	}

	double **At = new double *[ile_mam_rownan]; //alokacja pamieci
	for (int i = 0; i < ile_mam_rownan; ++i)
	{
		At[i] = new double[ile_mam_rownan]; //alokacja pamieci
		for (int j = 0; j < ile_mam_rownan; ++j)
			At[i][j] = 0;
	}


	double **C = new double *[ile_mam_rownan]; //alokacja pamieci
	for (int i = 0; i < ile_mam_rownan; ++i)
	{
		C[i] = new double[ile_mam_rownan]; //alokacja pamieci
		for (int j = 0; j < ile_mam_rownan; ++j)
			C[i][j] = 0;
	}

	double **d = new double *[ile_mam_rownan]; //alokacja pamieci
	for (int i = 0; i < ile_mam_rownan; ++i)
	{
		d[i] = new double[ile_mam_rownan]; //alokacja pamieci
		for (int j = 0; j < ile_mam_rownan; ++j)
			d[i][j] = 0;
	}



/////////////////////////////////////////////////////////
	if (TEST == 2)
	for (int i = 0; i < w; ++i)
	{
		for (int j = 0; j < k; ++j)
		
			cout<<" "<<tabAB[i][j] ;
		cout << endl;
	}

	double *X;
	X = new double[ile_mam_rownan];

	
	if (gauss(ile_mam_rownan, tabAB, X))
	{
		if (TEST==2 || TEST==1)
		for (int i = 0; i < ile_mam_rownan; i++)
			cout << "x" << i + 1 << " = " << setw(9) << X[i]
			<< endl;
	}
	else
		cout << "DZIELNIK ZERO\n";
	//////////////////////Gauss-Seidl///////////////////
	// A*x=b
	/////////////////////////////////

	
	//transpose A
	for (int i = 0; i<ile_mam_rownan; i++)
	{
		for (int j = 0; j<ile_mam_rownan; j++)
		{
			
			At[i][j] = A[j][i];
		}
	}

	//multiply At*A = C
	for (int i = 0; i<ile_mam_rownan; i++)
	{
		for (int j = 0; j<ile_mam_rownan; j++)
		{
			C[i][j] = 0;
			for (int k = 0; k<ile_mam_rownan; k++)
			{
				C[i][j] = C[i][j] + At[i][k] * A[k][j];
			}
		}
	}

	//multiply At*b = d
	for (int i = 0; i<ile_mam_rownan; i++)
	{
		for (int j = 0; j<1; j++)
		{
			d[i][j] = 0;
			for (int k = 0; k<ile_mam_rownan; k++)
			{
				d[i][j] = d[i][j] + At[i][k] * b[k][j];
			}
		}
	}

	//Solve Gauss-Seidel method for C*x = d
	//double x[N]; 
	//double p[N]; 
	double* x = new double[ile_mam_rownan];//current values
	double* p = new double[ile_mam_rownan];//previous values


	for (int i = 0; i<ile_mam_rownan; i++)
	{
		x[i] = 0;
		p[i] = 0;
	}

	do
	{
		for (int i = 0; i<ile_mam_rownan; i++)
			p[i] = x[i];

		for (int i = 0; i<ile_mam_rownan; i++)
		{
			double v = 0.0;
			for (int j = 0; j<i; j++)
			{
				double cij;
				if (i == j) cij = 0;
				else cij = -1.0*C[i][j] / C[i][i];
				v += cij*x[j];
			}

			for (int j = i; j<ile_mam_rownan; j++)
			{
				double cij;
				if (i == j) cij = 0;
				else cij = -1.0*C[i][j] / C[i][i];
				v += cij*p[j];
			}
			v += 1.0*d[i][0] / C[i][i];
			x[i] = v;
		}

	} while (!converge(x, p,ile_mam_rownan));

	//print solution
	printf("Err. val.: %.10lf\n", err);
	if(TEST==2)
	for (int i = 0; i<ile_mam_rownan; i++)
	{
		printf("%lf ", x[i]);
	}
	cout << endl;


	////////////////////////////////////////////////////
	printf("Gauss \t\t%.10lf \n", X[0]);
	printf("Gauss-Seidlera \t%.10lf \n", x[0]);
	for (int i = 1; i <= 100000; i *= 10)
	printf("%d Monte Carlo \t\t %.10lf \n",i, monte_carlo(ilosc_pol, plansza,i));



	system("pause");
	return 0;
}
