#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>  

using namespace std;

// Helper function to print a vector.
void print_vector(vector <float> V)
{
	for(vector <float> :: iterator it = V.begin(); it != V.end(); it++)
	{
		cout << *it << endl;
	}
	cout << "-------------------" << endl;
	cout << "Size of vector is: " << V.size() << endl;
}

// Normalisation to values between 5000 and -5000
void normalise(vector <float>& V)
{
	float max = 0;

	for(int i = 0; i<V.size(); i++)
	{
		if( abs( V[i] ) >= max )
		{
			max = abs( V[i] );
		}
	}

	float scaling_factor = 5000 / max;
	
	for(int i = 0; i<V.size(); i++)
	{
		V[i] = V[i] * scaling_factor;
	}
}

// DC Shift correction for the input
void DC_correction(vector <float>& V)
{
	float mean = ( accumulate(V.begin(),V.end(),0.0) / V.size() );

	for(int i = 0; i<V.size(); i++)
	{
		V[i] = V[i] - mean;
	}
}

// Levinson Durbin algorithm for Toeplitz matrices
vector<float> LD( vector<float> R)
{
    int p = R.size() - 1;
    
    vector<float> E( R.size(), 0);

    E[0]=R[0];
	int k=R.size();
    vector< vector<float> > a;
	vector<float> rowt;
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
    
    vector<float> alpha;
    for(int i=1;i<=p;i++)
    {
    	alpha.push_back( a[p][i] );
    }

    return alpha;

}

vector<float> getCepstral(vector<float> a)
{	
	vector<float> c(13,0.0);
	vector<float> A(13,0.0);
	for(int i=1;i<=12;i++){
		A[i]=a[i-1];
	}
	/*
	if(a.size()<12){
		cout<<"found error"<<endl;
		system("PAUSE");
	}
	*/
	c[1] = A[1];
	for(int m = 2; m <= 12 ; m++)
	{
		float term = 0.0;
		for(int k = 1 ; k <= m-1 ; k++)
		{
			//if((m-k)<a.size() && k<c.size())
			term = term +  (float)((k*c[k]*A[m-k])/m) ;
			/*
			else{
				cout<<"found error"<<endl;
				system("PAUSE");
			}
			*/
		}
		//if(m<c.size() && (m-1)<a.size() && m-1>=0)
			c[m] = A[m] + term;
		/*
		else{
			cout<<"val exceed"<<endl;
			system("PAUSE");
		}
		*/

	}
	c.erase(c.begin());
	/*
	vector<float> C;
	for(int i = 1 ; i <= 12 ; i++)
	{
		C.push_back( c[i] );
	}
	*/
	return c;

}

vector<float> raised_sine(vector<float> C)
{
	vector<float> wt;
	for(int i = 1 ; i <= 12 ; i++)
	{
		float arg = (3.14/12.0) * i; 
		float value = 1.0 + (6.0 * sin(arg) );
		wt.push_back( value );
	}

	//print_vector(wt);
	vector<float> C1( C.size() , 0);
	for(int i = 0 ; i < 12 ; i++)
	{
		C1[i] = C[i] * wt[i];
	}
	//print_vector(C1);

	return C1;
}

float tokura_dist(vector<float> v1, float v2[12])
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

// program starts here
int main()
{
	ifstream infile("t_rec/u.txt");
	//system("Recording_Module.exe 2 record.wav record.txt");
	
	//ifstream infile("record.txt");


	vector<float> input_vector;
	float val;

	// read the input as a text file line by line
	while (infile >> val)
	{
	    input_vector.push_back(val);
	}

	//DC correction for the input
	DC_correction(input_vector);

	// Normalisation
	normalise(input_vector);
	
	// all frames stores the data that will be used in processing
	vector<float> allframes;
	int start = (input_vector.size()/2) - 520;
	for( int i = 1 ; i <= 1040 ; i++ )
	{
		allframes.push_back( input_vector[start + i]);
	}
	input_vector.clear();

	// frame 2d vector has all the 10 frames we will process.
	vector<vector<float> > frame;
	for(int i = 0 ; i < 10 ; i++)
	{
		vector<float> tempv;
		for(int j = i*80 ; j < (i*80) + 320 ; j++)
		{
			tempv.push_back( allframes[j] );
		}
		frame.push_back( tempv );
	}
	allframes.clear();

	//apply Hamming window on each frame
	float hamming[320];
	for(int j = 0 ; j < frame[0].size() ; j++)
	{
		float weight = 0.54 - ( 0.46 * ( cos(44.0*j/7.0*319.0) ) );
		hamming[j] = weight;
	}
	for(int i = 0 ; i < 10 ; i++)
	{
		for(int j = 0 ; j < frame[i].size() ; j++)
		{
			frame[i][j] = frame[i][j] * hamming[j];
		}
	}

	//compute Ri's for i = 0 to 12 for all the 10 frames
	float r[10][13];
	for( int i = 0; i < 10 ; i++)
	{
		for( int j = 0 ; j <= 12 ; j++)
		{
			float tempr = 0;
			for(int k = 0 ; k <= 319 - j ; k++)
			{
				tempr = tempr + frame[i][k]*frame[i][k + j];
			}
			r[i][j] = tempr;
		}
	}
	
	//Levinson's Durbin Algorithm for ai's
	vector<float> Atemp;
	vector<vector<float> > Ais;
	
	for( int k = 0; k < 10 ; k++ )
	{
		vector<float> R;
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
	
	// for( int k = 0; k < 10 ; k++ )
	// {
	// 	print_vector( Ais[k] );
	// }

	//cepstral coefficients (ci's)
	vector<vector<float> > Cis;
	vector<float> Ctemp, Ctemp2;

	for( int k = 0; k < 10 ; k++ )
	{
		Ctemp = getCepstral( Ais[k] );
		Ctemp2 = raised_sine( Ctemp );
		Cis.push_back( Ctemp2 );
		Ctemp.clear();
		Ctemp2.clear();
	}

	// for( int k = 0; k < 10 ; k++ )
	// {
	// 	print_vector( Cis[k] );
	// }
	
	
	// referance values of vowel's cepstral coefficients.
	float a_cis[10][12] = {0.361571,2.38285,-0.971228,1.1227,-0.00394481,0.971213,-0.527497,-0.0972278,-0.234143,0.0602745,-0.00197246,0.0399727,
						0.373942,2.3995,-1.01009,1.0603,-0.0346926,0.996337,-0.545393,-0.0840947,-0.227043,-0.00869257,0.00793907,0.0586152,
						0.350343,2.31821,-1.02176,1.13372,0.0714701,1.19992,-0.248752,0.0425542,-0.188708,0.0619795,-0.00405603,0.0431927,
						0.352265,2.30672,-1.0456,1.04748,0.0726253,1.18845,-0.273536,-0.00352151,-0.170982,0.039469,-0.0253326,0.0380522,
						0.354253,2.29417,-1.11841,1.08752,0.0963564,1.35556,-0.0221697,0.199909,-0.102663,0.0795047,-0.0420437,0.0293309,
						0.338642,2.21061,-1.09257,1.11406,0.152678,1.28296,-0.129378,0.0971544,-0.100256,0.0912112,-0.0260164,0.0323945,
						0.403922,2.40573,-1.11927,1.07123,0.00591035,1.24526,-0.25196,0.056764,-0.210118,0.0282704,-0.0387888,0.040592,
						0.37032,2.26417,-1.14983,1.10004,0.141303,1.26581,-0.17796,0.125657,-0.118221,0.0951962,-0.0439744,0.0226634,
						0.380218,2.30041,-1.0777,1.1346,0.0522004,1.06946,-0.463979,-0.174027,-0.276743,0.0839738,-0.0119995,0.0311666,
						0.388868,2.36451,-1.03137,1.12016,-0.0667879,0.940547,-0.55757,-0.102572,-0.204476,0.104895,0.00991217,0.0443223};
	float e_cis[10][12] = {-0.024376,1.99327,-1.10629,1.34181,0.710555,2.41517,1.27375,1.05205,0.00435359,0.188666,-0.231932,-0.0631202,
						-0.0864282,1.95935,-1.13022,1.35392,0.800262,2.61946,1.38765,1.15088,0.064091,0.218413,-0.260994,-0.0621053,
						-0.212465,1.73863,-1.02289,1.414,0.972121,2.59154,1.47472,1.23142,0.152025,0.175374,-0.269164,-0.0609757,
						-0.260378,1.68636,-1.0513,1.54763,1.03084,2.54667,1.57562,1.3711,0.1869,0.195008,-0.274208,-0.0626803,
						-0.195577,1.8075,-0.959189,1.49011,0.706928,2.28489,1.52752,1.52842,0.230958,0.120237,-0.274522,-0.0540571,
						-0.328863,1.51255,-0.912038,1.53666,1.08262,2.43695,1.61932,1.36761,0.192672,0.103027,-0.248055,-0.0522524,
						-0.287601,1.71953,-0.94655,1.619,0.984151,2.53226,1.59957,1.37608,0.18632,0.105155,-0.254414,-0.04597,
						-0.258583,1.81395,-1.14793,1.48873,1.05657,2.93103,1.65328,1.22264,0.139014,0.281359,-0.285747,-0.0594675,
						-0.221241,1.78299,-0.996251,1.4548,0.893211,2.66817,1.53644,1.19819,0.164944,0.222623,-0.257861,-0.0620068,
						-0.177617,1.85398,-1.1179,1.44194,0.971152,2.78308,1.53316,1.11955,0.0315347,0.19234,-0.263838,-0.056995};
	float i_cis[10][12] = {-0.89503,1.23861,-0.382931,2.09431,1.2093,2.21262,2.38288,2.16972,0.28718,0.0951146,-0.344056,-0.022822,
						-0.881119,1.21703,-0.486012,2.12513,1.18833,2.18964,2.38182,2.1828,0.32275,0.11482,-0.353613,-0.0278368,
						-0.662552,1.4295,-0.656715,1.9027,0.974509,2.0935,2.20086,2.1062,0.289276,0.145824,-0.313004,-0.0409516,
						-0.705904,1.44643,-0.67074,1.94104,0.987709,2.18468,2.21777,2.13398,0.364712,0.220268,-0.354593,-0.0509646,
						-0.593789,1.53492,-0.716893,1.83382,0.906386,2.22407,2.14906,2.00466,0.266707,0.176645,-0.306605,-0.0437368,
						-0.738576,1.46304,-0.74849,1.93723,1.05333,2.31305,2.23515,2.13655,0.353132,0.346522,-0.381191,-0.0602151,
						-0.779554,1.44809,-0.75862,2.02893,1.13313,2.57736,2.24767,2.05756,0.35906,0.392838,-0.41397,-0.056787,
						-0.747095,1.46681,-0.774342,1.98369,1.16008,2.59818,2.2098,1.99995,0.34584,0.335742,-0.397436,-0.0499656,
						-0.649139,1.54827,-0.837171,1.87303,1.15257,2.5827,2.07402,1.95487,0.314006,0.305398,-0.362896,-0.0469099,
						-0.650556,1.55373,-0.814352,1.82818,1.04281,2.39273,2.08696,2.04581,0.389094,0.31084,-0.371383,-0.0575373};
	float o_cis[10][12] = {0.390147,2.63691,-1.26594,1.15677,0.164711,2.02194,0.510576,0.385878,-0.587139,0.0773115,-0.183306,-0.0303276,
						0.384297,2.65938,-1.29733,1.13787,0.145581,2.13979,0.691292,0.625865,-0.416612,0.126747,-0.19836,-0.0297558,
						0.325044,2.61034,-1.31867,1.21142,0.274931,2.33006,0.869054,0.783615,-0.335772,0.142464,-0.223591,-0.0344721,
						0.273969,2.59616,-1.33995,1.22623,0.33211,2.49198,1.05656,0.989741,-0.226991,0.198882,-0.244617,-0.0382301,
						0.261631,2.57824,-1.3294,1.24245,0.361961,2.49721,1.08472,1.01327,-0.228644,0.15009,-0.235838,-0.026317,
						0.231945,2.54014,-1.35906,1.22074,0.388701,2.56189,1.17553,1.05632,-0.180484,0.186157,-0.233773,-0.0302407,
						0.241225,2.56516,-1.37027,1.21471,0.37,2.55866,1.15519,1.04751,-0.177493,0.196558,-0.248217,-0.0394812,
						0.226666,2.55025,-1.37454,1.24852,0.427567,2.59127,1.16003,1.05647,-0.155335,0.164117,-0.239357,-0.0296077,
						0.251831,2.57595,-1.35282,1.25501,0.404991,2.53364,1.07979,0.997721,-0.19932,0.124528,-0.232559,-0.0245811,
						0.311565,2.62729,-1.34381,1.1867,0.279378,2.39918,0.96815,0.927499,-0.236974,0.153847,-0.21763,-0.0264214};
	float u_cis[10][12] = {-0.773246,1.97172,-1.24019,1.95748,1.3364,3.50362,1.98107,1.99805,0.594361,0.848202,-0.569263,-0.0476036,
						-0.80619,1.95987,-1.22832,1.99944,1.36008,3.53188,1.98999,1.9977,0.615285,0.880985,-0.589096,-0.0472462,
						-0.826624,1.94711,-1.21903,2.02091,1.3773,3.54545,1.99204,1.996,0.627914,0.892092,-0.599341,-0.0468107,
						-0.823165,1.95238,-1.21557,2.0011,1.38144,3.53515,1.98387,2.00611,0.622036,0.881657,-0.594386,-0.0467673,
						-0.737714,1.98472,-1.25664,1.93183,1.30214,3.49721,1.9714,1.98175,0.585885,0.827612,-0.555204,-0.0476452,
						-0.640901,2.03567,-1.29439,1.83752,1.1967,3.43036,1.95332,1.96461,0.547906,0.761859,-0.511566,-0.0475907,
						-0.630885,2.04799,-1.29441,1.80541,1.19415,3.39949,1.94226,1.98011,0.537775,0.752191,-0.50367,-0.0484554,
						-0.651424,2.03615,-1.28532,1.81615,1.21766,3.38351,1.95696,1.99182,0.535978,0.762252,-0.510209,-0.0490981,
						-0.649103,2.05004,-1.28477,1.80574,1.22288,3.37339,1.95342,2.00406,0.530669,0.772641,-0.509485,-0.0510065,
						-0.665917,2.04078,-1.26697,1.80152,1.2611,3.3846,1.93474,2.01909,0.532678,0.766716,-0.513425,-0.0506318};


	float d_a = 0;
	for(int i = 0 ; i < 10 ; i++)
	{
		//cout << tokura_dist( Cis[i] , a_cis[i]) << endl ;
		d_a += tokura_dist( Cis[i] , a_cis[i]);
	}

	float d_e = 0;
	for(int i = 0 ; i < 10 ; i++)
	{
		//cout << tokura_dist( Cis[i] , e_cis[i]) << endl ;
		d_e += tokura_dist( Cis[i] , e_cis[i]);
	}

	float d_i = 0;
	for(int i = 0 ; i < 10 ; i++)
	{
		d_i += tokura_dist( Cis[i] , i_cis[i]);
	}

	float d_o = 0;
	for(int i = 0 ; i < 10 ; i++)
	{
		d_o += tokura_dist( Cis[i] , o_cis[i]);
	}

	float d_u = 0;
	for(int i = 0 ; i < 10 ; i++)
	{
		d_u += tokura_dist( Cis[i] , u_cis[i]);
	}

	vector<float> distance_vector;
	distance_vector.push_back( d_a );
	distance_vector.push_back( d_e );
	distance_vector.push_back( d_i );
	distance_vector.push_back( d_o );
	distance_vector.push_back( d_u );

	// cout << d_a << endl;
	// cout << d_e << endl;
	// cout << d_i << endl;
	// cout << d_o << endl;
	// cout << d_u << endl;

	int argMin = distance(distance_vector.begin(), min_element(distance_vector.begin(), distance_vector.end()));

	switch( argMin )
	{
		case 0:{
			cout << "You spoke 'a' " << endl;
			break;
		}
		case 1:{
			cout << "You spoke 'e' " << endl;
			break;
		}
		case 2:{
			cout << "You spoke 'i' " << endl;
			break;
		}
		case 3:{
			cout << "You spoke 'o' " << endl;
			break;
		}
		case 4:{
			cout << "You spoke 'u' " << endl;
		}
	}
	
	return 0;
}