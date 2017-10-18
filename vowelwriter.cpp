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
    float a[R.size()][R.size()];
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

float tokura_dist(vector<float> v1, vector<float> v2)
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

	ofstream outfile;
	outfile.open ("wow.txt");
	for(int i = 0 ; i < 10 ; i++)
	{
		//outfile << "___" << i << "____" << endl;
		for(int j = 0 ; j < 12 ; j++)
		{
			outfile << Cis[i][j] << ",";
		}
		outfile << endl;
	}
	outfile.close();


}