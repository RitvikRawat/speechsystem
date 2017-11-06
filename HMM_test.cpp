#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <cfloat>
#include <algorithm>
#include <numeric>
#include <cmath>  

#define NO_OF_STATES 5
#define CODE_BOOK_SIZE 8

using namespace std;

struct HMM_instance{
	int N;		// no of states
	int M;		// no of distinct observations
	int T;		// length of observation sequence
	vector<vector<long double> > A;			//	(N+1)by(N+1)
	vector<vector<long double> > B;			//	(N+1)by(M+1)
	vector<vector<long double> > alpha;		//	(T+1)by(N+1)
	vector<vector<long double> > beta;		//	(T+1)by(N+1)
	vector<vector<long double> > gamma;		//	(T+1)by(N+1)
	vector<vector<long double> > delta;		//	(T+1)by(N+1)
	vector<vector<int> > psi;		//	(T+1)by(N+1)
	vector<long double> pi;
	vector<vector<vector<long double> > > zi; // (T+1)by(N+1)by(N+1)
};

// HMM_instance *HMM_1 = (HMM_instance*)malloc( sizeof(HMM_instance) );
HMM_instance *HMM_1 = new HMM_instance();

// Helper function to print a vector.
void print_vector(vector <long double> V)
{
    for(int i = 1 ; i < V.size() ; i++)
    {
        cout << V[i] << "\t";
    }
    cout << endl <<"-------------------" << endl;
    // cout << "Size of vector is: " << V.size() << endl;
}

void print_vector_INT(vector <int> V)
{
    for(vector <int> :: iterator it = V.begin(); it != V.end(); it++)
    {
        cout << *it << "\t";
    }
    cout << "-------------------" << endl;
    //cout << "Size of vector is: " << V.size() << endl;
}

void print_2D_vector(vector<vector<long double> > V){
	for(int i = 1 ; i < V.size() ; i++){
		for(int j = 1 ; j < V[i].size() ; j++){
			cout << V[i][j] << "\t" ;
		}
		cout << endl;
	}
	cout <<"_________________"<< endl;
}

void print_3D_vector(vector<vector<vector<long double> > > Z){
	for (int i = 1; i < Z.size(); ++i)
	{
		for (int j = 1; j < Z[0].size(); ++j)
		{
			for (int k = 1; k < Z[0][0].size(); ++k)
			{
				cout << Z[i][j][k] << "\t";
			}
			cout << endl;
		}
		cout << "................" << endl;
	}
}

void print_2D_vector_INT(vector<vector<int> > V){
	for(int i = 1 ; i < V.size() ; i++){
		for(int j = 1 ; j < V[i].size() ; j++){
			cout << V[i][j] << "\t\t" ;
		}
		cout << endl;
	}
}

// Normalisation to values between 5000 and -5000
void normalise(vector <long double>& V)
{
	long double max = 0;

	for(int i = 0; i<V.size(); i++)
	{
		if( abs( V[i] ) >= max )
		{
			max = abs( V[i] );
		}
	}

	long double scaling_factor = 5000 / max;
	
	for(int i = 0; i<V.size(); i++)
	{
		V[i] = V[i] * scaling_factor;
	}
}

// DC Shift correction for the input
void DC_correction(vector <long double>& V)
{
	long double mean = ( accumulate(V.begin(),V.end(),0.0) / V.size() );

	for(int i = 0; i<V.size(); i++)
	{
		V[i] = V[i] - mean;
	}
}

vector<long double> readtext(string s){
	ifstream infile( s.c_str() );

	vector <long double> input_vector;
	long double val;
	// read the input as a text file line by line
	while (infile >> val)
	{
	    input_vector.push_back(val);
	}
	//DC correction for the input
	DC_correction(input_vector);
	// Normalisation
	normalise(input_vector);
	return input_vector;
}

vector<vector<long double> > readCSV(string s){
	vector<vector<long double> > A;
	vector<long double> row;
	ifstream infile( s.c_str() );
    string line = "";
    while (getline(infile, line)){
        stringstream strstr(line);
        string word = "";
		while (getline(strstr,word, ',')) row.push_back(stof(word));
		A.push_back(row);
		row.clear();
    }
	return A;
	
}

vector<vector<long double> > getframes(vector<long double> in){
	vector<vector<long double> > frame;
	for(int i = 0 ; ((i*80) + 320) < in.size() ; i++)
	{
		vector<long double> tempv;
		for(int j = i*80 ; j < (i*80) + 320 ; j++)
		{
			tempv.push_back( in[j] );
		}
		frame.push_back( tempv );
	}
	return frame;
}

void apply_hamming(vector<vector<long double> > &frame){
	long double hamming[320];
	for(int j = 0 ; j < frame[0].size() ; j++)
	{
		long double weight = 0.54 - ( 0.46 * ( cos(44.0*j/7.0*319.0) ) );
		hamming[j] = weight;
	}
	for(int i = 0 ; i < 10 ; i++)
	{
		for(int j = 0 ; j < frame[i].size() ; j++)
		{
			frame[i][j] = frame[i][j] * hamming[j];
		}
	}
}

vector<vector<long double> > getRis(vector<vector<long double> > frame){
	//long double r[10][13];
	vector<vector<long double> > r;
	vector<long double> row(13, 0.0);
	for(int i = 0 ; i < frame.size() ; i++){
		r.push_back( row );
	}

	for( int i = 0; i < frame.size() ; i++)
	{
		for( int j = 0 ; j <= 12 ; j++)
		{
			long double tempr = 0;
			for(int k = 0 ; k <= 319 - j ; k++)
			{
				tempr = tempr + frame[i][k]*frame[i][k + j];
			}
			r[i][j] = tempr;
		}
	}
	return r;
}

// Levinson Durbin algorithm for Toeplitz matrices
vector<long double> LD( vector<long double> R)
{
    int p = R.size() - 1;
    
    vector<long double> E( R.size(), 0);

    E[0]=R[0];
	int k=R.size();
    vector< vector<long double> > a;
	vector<long double> rowt;
	for(int i=0;i<k;i++){
		rowt.push_back(0.0);
	}
	for(int i=0;i<k;i++){
		a.push_back(rowt);
	}
    E[1]=(1-a[1][1]*a[1][1])*E[0];
    
    for(int i=2;i<=p;i++)
    {
        a[i][i]=R[i];
        for(int j=1;j<=i-1;j++)
        {
            a[i][i]-=a[i-1][j]*R[i-j];
        }
        
        if(E[i-1]!=0)
            a[i][i]/=E[i-1];
        
        for(int j=1;j<=i-1;j++)
        {
            a[i][j]=a[i-1][j]-a[i][i]*a[i-1][i-j];
        }
        
        E[i]=(1-a[i][i]*a[i][i])*E[i-1];
    
    }
    
    vector<long double> alpha;
    for(int i=1;i<=p;i++)
    {
    	alpha.push_back( a[p][i] );
    }

    return alpha;

}

vector<long double> getCepstral(vector<long double> a)
{	
	vector<long double> c(13,0.0);
	vector<long double> A(13,0.0);
	for(int i=1;i<=12;i++){
		A[i]=a[i-1];
	}
	c[1] = A[1];
	for(int m = 2; m <= 12 ; m++)
	{
		long double term = 0.0;
		for(int k = 1 ; k <= m-1 ; k++)
		{
			term = term +  (long double)((k*c[k]*A[m-k])/m) ;	
		}
			c[m] = A[m] + term;
	}
	c.erase(c.begin());
	return c;

}

vector<long double> raised_sine(vector<long double> C)
{
	vector<long double> wt;
	for(int i = 1 ; i <= 12 ; i++)
	{
		long double arg = (3.14/12.0) * i; 
		long double value = 1.0 + (6.0 * sin(arg) );
		wt.push_back( value );
	}

	//print_vector(wt);
	vector<long double> C1( C.size() , 0);
	for(int i = 0 ; i < 12 ; i++)
	{
		C1[i] = C[i] * wt[i];
	}
	//print_vector(C1);
	return C1;
}

float tokura_dist(vector<long double> v1, vector<long double> v2)
{
	float w[13];
	w[0]=0;
	w[1]=1;
	w[2]=3;
	w[3]=7;
	w[4]=13;
	w[5]=19;
	w[6]=22;
	w[7]=25;
	w[8]=33;
	w[9]=42;
	w[10]=50;
	w[11]=56;
	w[12]=61;

	float dist = 0;
	for(int i = 0 ; i < 12 ; i++)
	{
		dist += ((float) w[i+1] * (v1[i]-v2[i]) * (v1[i]-v2[i]) );
	}
	return dist;
}

int getMinDist(vector<long double> v, vector<vector<long double> > c_book){
	int min_index = 0;
	long double dist = FLT_MAX, min_dist = FLT_MAX;
	for(int i = 0 ; i < c_book.size() ; i++){
		dist = tokura_dist( v, c_book[i] );
		if(dist < min_dist){
			min_index = i;
			min_dist = dist;
		}
	}
	return min_index;
}

vector<int> vector_quantize(vector<long double> in){
	vector<vector<long double> > codebook = readCSV("LBG_codebook_8.csv");
	vector<vector<long double> > frame = getframes(in);
	apply_hamming(frame);
	vector<vector<long double> > r = getRis(frame);
	
	//Levinson's Durbin Algorithm for ai's
	vector<vector<long double> > Ais; vector<long double> Atemp;
	for( int k = 0; k < frame.size() ; k++ )
	{
		vector<long double> R;
		for(int j = 0 ; j <= 12 ; j++)
		{
			R.push_back( r[k][j] );
		}
		//function returns the list of ai's
		Atemp = LD(R);
		//algo ends
		R.clear();

		Ais.push_back( Atemp );
		Atemp.clear();
	}

	vector<vector<long double> > Cis;
	vector<long double> Ctemp, Ctemp2;
	for( int k = 0; k < frame.size() ; k++ )
	{
		Ctemp = getCepstral( Ais[k] );
		Ctemp2 = raised_sine( Ctemp );
		Cis.push_back( Ctemp2 );
		Ctemp.clear();
		Ctemp2.clear();
	}

	vector<int> O;
	// NOTE that codebook index is from 0 but I use obs from 1 to 8
	for(int i = 0 ; i < Cis.size() ; i++){
		O.push_back( ( getMinDist( Cis[i] , codebook ) )+1 );
	}

	// for(int i = 0 ; i < O.size() ; i++){
	// 	cout << O[i] << endl;
	// }
	//print_vector_INT(O);
	//cout << O.size() <<"_"<< Cis.size() << endl;
	// print_2D_vector(Cis);
	// cout << Cis.size() << endl;

	return O;
}


vector<int> read_obs_seq(){
	vector<long double> in = readtext("input_signal.txt");
	vector<int> O = vector_quantize(in);
	vector<int> V(3, 3);
	// V[0] = 1; V[1] = 2; V[2] = 1;
	return O;
}

vector<vector<long double> > zero_2D(int r, int c){
	vector<long double> row(c,0.0);
	vector<vector<long double> > V;
	for(int i = 0 ; i < r ; i++){
		V.push_back( row );
	}
	return V;
}

vector<vector<vector<long double> > > zero_3D(int x, int y, int z){
	vector<long double> row(z,0.0);
	vector<vector<long double> > grid;
	for(int i = 0 ; i < y ; i++){
		grid.push_back( row );
	}
	vector<vector<vector<long double> > > cube;
	for (int i = 0; i < x; i++){
		cube.push_back( grid );
	}
	return cube;
}

vector<vector<int> > zero_2D_INT(int r, int c){
	vector<int> row(c,0);
	vector<vector<int> > V;
	for(int i = 0 ; i < r ; i++){
		V.push_back( row );
	}
	return V;
}

vector<long double> zero_1D(int size){
	vector<long double> V(size, 0);
	return V;
}

void populate_pi( vector<long double> &pi){
	pi[1] = 1; pi[2] = 0; pi[3] = 0; pi[4] = 0; pi[5] = 0;
}

void populate_A( vector<vector<long double> > &A ){
	int N = A.size()-1;
	for(int i = 1 ; i <= N ; i++){
		for(int j = 1 ; j <= N ; j++){
			if(i == 5 && j == 5){
				A[i][j] = 1;
			}
			else if(i == j){
				A[i][j] = 0.9;
				A[i][j+1] = 0.1;
			}
		}
	}
}

void populate_B( vector<vector<long double> > &B ){
	int N = B.size()-1, M = B[0].size()-1;
	 for(int i = 1 ; i <= N ; i++){
		for(int j = 1 ; j <= M ; j++){
			B[i][j] = 0.125;
		}
	}
}

long double forward(vector<long double> pi, vector<vector<long double> > A, vector<vector<long double> > B, vector<vector<long double> > &alpha, vector<int> o){
	int N = alpha[0].size()-1, T = alpha.size()-1;
	for(int i = 1 ; i <= N ; i++){
		alpha[1][i] = pi[i] * B[i][ o[1] ];
	}
	for( int t = 1 ; t <= T-1 ; t++){
		for(int j = 1 ; j <= N ; j++){
			long double temp = 0;
			for(int i = 1 ; i <= N ; i++){
				temp += ( alpha[t][i] * A[i][j] );
			}
			alpha[t+1][j] = temp * B[j][ o[t+1] ];
		}
	}
	long double prob = 0;
	for(int i = 1 ; i <= N ; i++){
		prob += ( alpha[T][i] );
	}
	return prob;
}

void backward(vector<vector<long double> > A, vector<vector<long double> > B, vector<vector<long double> > &beta, vector<int> o){
	int N = beta[0].size()-1, T = beta.size()-1;
	for(int i = 1; i <= N ; i++){
		beta[T][i] = 1;
	}
	for(int t = T-1 ; t >= 1 ; t--){
		for(int i = 1; i <= N ; i++){
			long double temp = 0;
			for(int j = 1 ; j <= N ; j++){
				temp += ( A[i][j] * B[j][ o[t+1] ] * beta[t+1][j] );
			}
			beta[t][i] = temp;
		}
	}
}

void gamma_calculator(vector<vector<long double> > alpha, vector<vector<long double> > beta, vector<vector<long double> > &gamma){
	int N = alpha[0].size()-1, T = alpha.size()-1;
	for(int t = 1; t <= T ; t++){
		long double norm_factor = 0;
		for(int i = 1; i <= N ; i++){
			norm_factor += ( alpha[t][i] * beta[t][i] );
		}
		for(int i = 1; i <= N ; i++){
			gamma[t][i] = ( alpha[t][i] * beta[t][i] ) / norm_factor;
		}
	}
}

vector<int> viterbi(vector<vector<long double> > A, vector<vector<long double> > B, vector<long double> pi, vector<vector<long double> > &delta, vector<vector<int> > &psi, vector<int> o){
	int N = delta[0].size()-1, T = delta.size()-1;
	for(int i = 1 ; i <= N ; i++){
		delta[1][i] = pi[i] * B[i][ o[1] ];
		psi[1][i] = 0;
	}
	for(int t = 2 ; t <= T ; t++){
		for(int j = 1 ; j <= N ; j++){
			long double max_prod = delta[t-1][1] * A[1][j];
			int arg_max_prod = 1;
			for(int i = 1 ; i <= N ; i++){
				if( delta[t-1][i] * A[i][j] >= max_prod ){
					max_prod = delta[t-1][i] * A[i][j];
					arg_max_prod = i;
				}
			}
			delta[t][j] = max_prod * B[j][ o[t] ];
			psi[t][j] = arg_max_prod;
		}
	}
	vector<int> seq(T+1,-1);
	long double P_max = -1;
	for(int i = 1 ; i <= N ; i++){
		if( delta[T][i] > P_max ){
			P_max = delta[T][i];
			seq[T] = i;
		}
	}
	for(int t = T-1 ; t >= 1 ; t--){
		seq[t] = psi[t+1][ seq[t+1] ];
	}
	return seq;
}

void zi_calculator(vector<vector<long double> > A, vector<vector<long double> > B, vector<vector<long double> > alpha, vector<vector<long double> > beta, vector<vector<long double> > gamma, vector<vector<vector<long double> > > &zi, vector<int> o){
	int N = alpha[0].size()-1, T = alpha.size()-1;
	for(int t = 1 ; t <= T-1 ; t++){
		long double den = 0;
		for(int i = 1 ; i <= N ; i++){
			for(int j = 1 ; j <= N ; j++){
				den += ( alpha[t][i] * A[i][j] * B[j][ o[t+1] ] * beta[t+1][j] );
								 
			}
		}
		for(int i = 1 ; i <= N ; i++){
			for(int j = 1 ; j <= N ; j++){
				long double num = alpha[t][i] * A[i][j] * B[j][ o[t+1] ] * beta[t+1][j];
				zi[t][i][j] = num / den ;
			}
		}
	}
}

void EM_optimiser(vector<long double> &pi, vector<vector<long double> > &A, vector<vector<long double> > &B, vector<vector<long double> > gamma, vector<vector<vector<long double> > > zi, vector<int> o){
	int N = gamma[0].size()-1, T = gamma.size()-1, M = B[0].size()-1;
	for(int i = 1 ; i <= N ; i++){
		pi[i] = gamma[1][i];
	}
	for(int i = 1 ; i <= N ; i++){
		for(int j = 1 ; j <= N ; j++){
			long double num = 0;
			for(int t = 1 ; t <= T-1 ; t++){
				num += zi[t][i][j];
			}
			long double den = 0;
			for(int t = 1 ; t <= T-1 ; t++){
				den += gamma[t][i];
			}
			A[i][j] = num / den;
			// cout << A[i][j] << " is "<< num << " on " << den << endl;
		}
	}
	for(int j = 1 ; j <= N ; j++){
		for(int k = 1 ; k <= M ; k++){
			long double den = 0;
			for(int t = 1; t <= T-1 ; t++){
				den += gamma[t][j];
			}
			long double num = 0;
			for(int t = 1; t <= T-1 ; t++){
				if( o[t] == k ){
					num += gamma[t][j];
				}
			}
			B[j][k] = num / den;
		}
	}
}

// Helper function to print a vector.
void log_vector(vector <long double> V , string s){
	ofstream outfile;
	outfile.open ( s.c_str() , ofstream::out | ofstream::app);
	outfile.precision(17);

    for(int i = 1 ; i < V.size() ; i++)
    {
        outfile << V[i] << "\t";
    }
    outfile << endl <<"-------------------" << endl;
    // outfile << "Size of vector is: " << V.size() << endl;
    outfile.close();

}

void log_vector_INT(vector <int> V , string s){
	ofstream outfile;
	outfile.open ( s.c_str() , ofstream::out | ofstream::app);
	outfile.precision(17);

    for(vector <int> :: iterator it = V.begin()+1; it != V.end(); it++)
    {
        outfile << *it << endl;
    }
    outfile << "-------------------" << endl;
    //outfile << "Size of vector is: " << V.size() << endl;
    outfile.close();

}

void log_2D_vector(vector<vector<long double> > V , string s){
	ofstream outfile;
	outfile.open ( s.c_str() , ofstream::out | ofstream::app);
	outfile.precision(17);

	for(int i = 1 ; i < V.size() ; i++){
		for(int j = 1 ; j < V[i].size() ; j++){
			outfile << V[i][j] << "\t" ;
		}
		outfile << endl;
	}
	outfile <<"_________________"<< endl;
	outfile.close();

}

void log_3D_vector(vector<vector<vector<long double> > > Z , string s){
	ofstream outfile;
	outfile.open ( s.c_str() , ofstream::out | ofstream::app);
	outfile.precision(17);

	for (int i = 1; i < Z.size(); ++i)
	{
		for (int j = 1; j < Z[0].size(); ++j)
		{
			for (int k = 1; k < Z[0][0].size(); ++k)
			{
				outfile << Z[i][j][k] << "\t";
			}
			outfile << endl;
		}
		outfile << "................" << endl;
	}
	outfile.close();

}

void log_2D_vector_INT(vector<vector<int> > V , string s){
	ofstream outfile;
	outfile.open ( s.c_str() , ofstream::out | ofstream::app);
	outfile.precision(17);

	for(int i = 1 ; i < V.size() ; i++){
		for(int j = 1 ; j < V[i].size() ; j++){
			outfile << V[i][j] << "\t\t" ;
		}
		outfile << endl;
	}
	outfile.close();
}

void modify_B(vector<vector<long double> > &B){
	int N = B.size()-1, M = B[0].size()-1;
	for(int i = 1 ; i <= N ; i++){
		long double row_min = 10.0e-5;
		int noZ = 0;
		long double row_max = 0;
		int max_index = 0;
		for(int j = 1 ; j <= M ; j++){
			if(B[i][j] == 0){
				noZ++;
			}
			if(B[i][j] >= row_max){
				row_max = B[i][j];
				max_index = j;
			}
			
		}
		for(int j = 1 ; j <= M ; j++){
			if(B[i][j] == 0){
				B[i][j] = row_min;
			}
			if(j == max_index){
				B[i][j] = B[i][j] - (noZ * row_min);
			}
		}
	}
}

int main(){

	// get observation sequence
	vector<int> in_obs = read_obs_seq();
	in_obs.insert(in_obs.begin(), -1);

	// model parameters memory allocation
	HMM_1->N = NO_OF_STATES;
	HMM_1->M = CODE_BOOK_SIZE;
	HMM_1->T = in_obs.size()-1; //we have a dummy -1 in start to begin indexing from 1
	HMM_1->A = zero_2D( HMM_1->N + 1 , HMM_1->N + 1 );
	HMM_1->B = zero_2D( HMM_1->N + 1 , HMM_1->M + 1 );
	HMM_1->pi = zero_1D( HMM_1->N + 1 );
	// HMM_1->alpha = zero_2D( HMM_1->T + 1 , HMM_1->N + 1 );
	// HMM_1->beta = zero_2D( HMM_1->T + 1 , HMM_1->N + 1 );
	// HMM_1->gamma = zero_2D( HMM_1->T + 1 , HMM_1->N + 1 );
	// HMM_1->delta = zero_2D( HMM_1->T + 1 , HMM_1->N + 1 );
	// HMM_1->psi = zero_2D_INT( HMM_1->T + 1 , HMM_1->N + 1 );
	// HMM_1->zi = zero_3D( HMM_1->T + 1 , HMM_1->N + 1 , HMM_1->N + 1 );

	// populate pi, A and B
	// populate_pi( HMM_1->pi );
	// populate_A( HMM_1->A );
	// populate_B( HMM_1->B );

	//cout << "Initial Prob( O | HMM_1 ) = " << forward(HMM_1->pi, HMM_1->A, HMM_1->B, HMM_1->alpha, in_obs) << endl;
	// forward(HMM_1->pi, HMM_1->A, HMM_1->B, HMM_1->alpha, in_obs);
	//backward(HMM_1->A, HMM_1->B, HMM_1->beta, in_obs);
	//gamma_calculator(HMM_1->alpha, HMM_1->beta, HMM_1->gamma);
	//vector<int> seq;
	//seq = viterbi(HMM_1->A, HMM_1->B, HMM_1->pi, HMM_1->delta, HMM_1->psi, in_obs);

	// cout << "Most likely state sequence : " << endl;
	// for(int i = 1 ; i < seq.size() ; i++){
	// 	cout << seq[i] << endl;
	// }
	// cout <<"___________________________"<< endl;

	//zi_calculator(HMM_1->A, HMM_1->B,HMM_1->alpha, HMM_1->beta, HMM_1->gamma, HMM_1->zi, in_obs);

	cout.precision(4);

	int num_of_iter = 5;
	for(int i = 0 ; i <= num_of_iter ; i++){
		if(i == 0){//hardcoded values
			populate_pi( HMM_1->pi );
			populate_A( HMM_1->A );
			populate_B( HMM_1->B );
		}
		else{// EM algorithm
			EM_optimiser(HMM_1->pi, HMM_1->A, HMM_1->B, HMM_1->gamma, HMM_1->zi, in_obs);
			modify_B(HMM_1->B);
		}
		HMM_1->alpha.clear(); 
		HMM_1->alpha = zero_2D( HMM_1->T + 1 , HMM_1->N + 1 );
		HMM_1->beta.clear(); 
		HMM_1->beta = zero_2D( HMM_1->T + 1 , HMM_1->N + 1 );
		HMM_1->gamma.clear(); 
		HMM_1->gamma = zero_2D( HMM_1->T + 1 , HMM_1->N + 1 );
		HMM_1->delta.clear(); 
		HMM_1->delta = zero_2D( HMM_1->T + 1 , HMM_1->N + 1 );
		HMM_1->psi.clear(); 
		HMM_1->psi = zero_2D_INT( HMM_1->T + 1 , HMM_1->N + 1 );
		HMM_1->zi.clear(); 
		HMM_1->zi = zero_3D( HMM_1->T + 1 , HMM_1->N + 1 , HMM_1->N + 1 );

		cout << "Iter # "<<i<<"--> Prob( O | HMM_1 ) = " << forward(HMM_1->pi, HMM_1->A, HMM_1->B, HMM_1->alpha, in_obs) << endl;
		backward(HMM_1->A, HMM_1->B, HMM_1->beta, in_obs);
		gamma_calculator(HMM_1->alpha, HMM_1->beta, HMM_1->gamma);
		vector<int> seq;
		seq = viterbi(HMM_1->A, HMM_1->B, HMM_1->pi, HMM_1->delta, HMM_1->psi, in_obs);
		zi_calculator(HMM_1->A, HMM_1->B,HMM_1->alpha, HMM_1->beta, HMM_1->gamma, HMM_1->zi, in_obs);

		ofstream fp; fp.open( "log_test.txt" , ofstream::out | ofstream::app );
		fp.precision(17);
		fp << "______________________________________________________" << endl;
		fp << "----------->>> ITERATION NUMBER # "<<i<<" LOG <<<-----------" << endl;
		fp << "------------------------------------------------------"<<endl;
		fp << endl << "pi is : " << endl;
		log_vector(HMM_1->pi , "log_test.txt");
		fp << endl << "A is : " << endl;
		log_2D_vector(HMM_1->A , "log_test.txt");
		fp << endl << "B is : " << endl;
		log_2D_vector(HMM_1->B , "log_test.txt");
		fp << endl << "alpha is : " << endl;
		log_2D_vector(HMM_1->alpha , "log_test.txt");
		fp << endl << "beta is : " << endl;
		log_2D_vector(HMM_1->beta , "log_test.txt");
		fp << endl << "gamma is : " << endl;
		log_2D_vector(HMM_1->gamma , "log_test.txt");
		fp << endl << "delta is : " << endl;
		log_2D_vector(HMM_1->delta , "log_test.txt");
		fp << endl << "psi is : " << endl;
		log_2D_vector_INT(HMM_1->psi , "log_test.txt");
		fp << endl << "most likely sequence is :" << endl;
		log_vector_INT(seq, "log_test.txt");
		fp.close();

	}

	//logger();
	
	// print_vector(HMM_1->pi);
	// print_2D_vector(HMM_1->A);
	// print_2D_vector(HMM_1->B);
	// print_2D_vector(HMM_1->alpha);
	// print_2D_vector(HMM_1->beta);
	// print_2D_vector(HMM_1->gamma);
	// print_2D_vector(HMM_1->delta);
	// print_2D_vector_INT(HMM_1->psi);
	// print_3D_vector(HMM_1->zi);

	// EM_optimiser(HMM_1->pi, HMM_1->A, HMM_1->B, HMM_1->gamma, HMM_1->zi, in_obs);

	// print_vector(HMM_1->pi);
	// print_2D_vector(HMM_1->A);
	// print_2D_vector(HMM_1->B);

	// print_3D_vector(HMM_1->zi);


	/*
	ofstream fp; fp.open( "log_test.txt" , ofstream::out | ofstream::app );
	fp.precision(17);
	fp << endl << "pi is : " << endl;
	log_vector(HMM_1->pi , "log_test.txt");
	fp << endl << "A is : " << endl;
	log_2D_vector(HMM_1->A , "log_test.txt");
	fp << endl << "B is : " << endl;
	log_2D_vector(HMM_1->B , "log_test.txt");
	fp << endl << "alpha is : " << endl;
	log_2D_vector(HMM_1->alpha , "log_test.txt");
	fp << endl << "beta is : " << endl;
	log_2D_vector(HMM_1->beta , "log_test.txt");
	fp << endl << "gamma is : " << endl;
	log_2D_vector(HMM_1->gamma , "log_test.txt");
	fp << endl << "delta is : " << endl;
	log_2D_vector(HMM_1->delta , "log_test.txt");
	fp << endl << "psi is : " << endl;
	log_2D_vector_INT(HMM_1->psi , "log_test.txt");
	fp << endl << "most likely sequence is :" << endl;
	log_vector_INT(seq, "log_test.txt");
	fp.close();
	*/

	return 0;
}