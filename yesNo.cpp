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

// function computer Short Term Energy as described in class for 100 sample intervals
void get_STE(vector <float> input_vector, vector <float>& A)
{
	float STE_100_sample;

	for(int i = 0; i < input_vector.size(); i+=100)
	{
		STE_100_sample = 0;
		for(int j = i; (j < i + 100)&&(j < input_vector.size()); j++)
		{
			STE_100_sample += ( (input_vector[j]) * (input_vector[j]) );
		}

		A.push_back(STE_100_sample);
	}
}

// function computer Zero Crossing Rate as described in class for 100 sample intervals
void get_ZCR(vector <float> input_vector, vector <float>& B)
{
	float ZCR_100_sample;
	int current_positive;

	for(int i = 0; i < input_vector.size(); i+=100)
	{
		ZCR_100_sample = 0;
		if(input_vector[i] >= 0)
		{
			current_positive = 1;
		}
		else
		{
			current_positive = 0;
		}

		for(int j = i+1; (j < i + 100)&&(j < input_vector.size()); j++)
		{
			if(current_positive == 1)
			{
				if(input_vector[j] < 0)
				{
					ZCR_100_sample++;
					current_positive = 0;
				}
			}
			else
			{
				if(input_vector[j] >= 0)
				{
					ZCR_100_sample++;
					current_positive = 1;
				}
			}
		}

		B.push_back(ZCR_100_sample);
	}
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

// program starts here
int main()
{
	ifstream infile("rr_yes_1.txt");

	vector <float> input_vector;
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
	
	// STE calculation for 100 samples groups stored in A
	vector <float> A;
	get_STE(input_vector, A);

	// processing and comparision on basis of STE
	float max_STE = *max_element(A.begin(), A.end());

	// offset that seperates the front 2/3 of vector from the last 1/3
	int offset = 2 * A.size() / 3;
	int N = A.size() / 3;

	float mean_STE = ( accumulate(A.begin() + offset, A.end() ,0.0) ) / N;

	int is_yes, is_no;

	// this check gives 1 weightage to the sample being a "NO" if the initial 1/3 of the 
	// STE vector had a value greater than 0.06 of the maximum STE (the "ooooooo" sound is captured)
	// The value of 0.06 was determined by testing on multiple samples of my voice.
	if( mean_STE < 0.06 * max_STE)
	{
		is_yes = 1;
		is_no = 0;
	}
	else
	{
		is_yes = 0;
		is_no = 1;	
	}

	// ZCR calculation for 100 samples groups stored in B
	vector <float> B;
	get_ZCR(input_vector, B);
	
	// processing and comparision on basis of ZCR
	offset = 7 * B.size() / 12;
	
	// the mean for initial 7 / 12 of the ZCR values.
	N = 7 * B.size() / 12;
	float m1 =  ( accumulate(B.begin(), B.begin() + offset ,0.0) ) / N;
	
	// the mean for later 5 / 12 of the ZCR values.
	N = 5 * B.size() / 12;
	float m2 = ( accumulate(B.begin() + offset , B.end() ,0.0) ) / N;

	// if the later part has a higher mean ZCR then attribute 1 weightage to the sample being a "NO"
	// the "sssssssss" part of YES is captured.
	if(m2 > m1)
	{
		is_yes++;
	}
	else
	{
		is_no++;
	}

	// Final yes/no decision. If both ZCR and STE point to YES or NO then it is the output.
	// If the word said is ambiguous then "Ambiguous input" is reported.
	if( is_yes == 2)
	{
		cout << "You said YES." << endl;
	}
	else if( is_no == 2)
	{
		cout << "You said NO" << endl;
	}
	else
	{
		cout << "Ambiguous input" << endl;
	}
}