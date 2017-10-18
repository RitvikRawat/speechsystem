#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <cfloat>
#include <algorithm>

using namespace std;

float min(float a, float b){
	if(a <= b){
		return a;
	}
	else{
		return b;
	}
}
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
void print_2D_vector(vector<vector<float> > V){
	for(int i = 0 ; i < V.size() ; i++){
		for(int j = 0 ; j < V[i].size() ; j++){
			cout << V[i][j] << "__" ;
		}
		cout << endl;
	}
}
vector<vector<float> > getuniverse(){
	vector<vector<float> > A;
	vector<float> row;
	ifstream infile("Universe.csv");
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
float getDistortion(vector<vector<float> > in, vector<vector<float> > c_book){
	float global_dist = 0, vector_dist = 0;
	for(int i = 0 ; i < in.size() ; i++){
		vector_dist = FLT_MAX;
		for(int j = 0 ; j < c_book.size() ; j++){
			vector_dist = min( vector_dist , tokura_dist( in[i], c_book[j] ) );
		}
		global_dist += vector_dist;
	}
	return global_dist;
}
int getMinDist(vector<float> v, vector<vector<float> > c_book){
	int min_index = 0;
	float dist = FLT_MAX, min_dist = FLT_MAX;
	for(int i = 0 ; i < c_book.size() ; i++){
		dist = tokura_dist( v, c_book[i] );
		if(dist < min_dist){
			min_index = i;
			min_dist = dist;
		}
	}
	return min_index;
}
vector<vector<float> > K_means(vector<vector<float> > in, vector<vector<float> > c_book, vector<float> &dist){
	vector<int> map( in.size() , 0 );
	vector<float> D;
	//initial distortion
	D.push_back( getDistortion(in, c_book) );
	int epochs = 0;
	while(epochs < 40){
		// populate the mapping 
		for(int i = 0 ; i < in.size() ; i++){
			int cluster_index = getMinDist( in[i], c_book );
			map[i] = cluster_index;
		}
		
		// for(int i = 0 ; i < map.size() ; i++){
		// 	cout << map[i] << endl;
		// }
		// temp matrix
		vector<vector<float> > temp_c_book;
		// cluster_size stores the number of vectors in this bucket
		vector<float> cluster_size( c_book.size(), 0);
		vector<float> row( 12 , 0);
		for(int i = 0 ; i < c_book.size() ; i++){
			temp_c_book.push_back(row);
		}
		// get sum of all vectors in a bucket store in temp_c_book
		for(int i = 0 ; i < in.size() ; i++){
			int cbi = map[i];
			transform( temp_c_book[cbi].begin(), temp_c_book[cbi].end(), in[i].begin(), temp_c_book[cbi].begin(), plus<float>());
			cluster_size[cbi]++;
		}
		// get centroid by dividing with cluster_size
		for(int i = 0 ; i < temp_c_book.size() ; i++){
			for(int j = 0 ; j < 12 ; j++){
				temp_c_book[i][j] /= cluster_size[i]; 
			}
		}
		c_book = temp_c_book;
		temp_c_book.clear();
		cluster_size.clear();
		D.push_back( getDistortion(in, c_book) );
		if( D[ D.size()-1 ] > 0.99 * D[ D.size()-2 ]){
			break;
		}
		epochs++;
	}
	//print_vector(D);
	for(int i = 0 ; i < D.size() ; i++){
		dist.push_back( D[i] );
	}
	return c_book;
}
vector<vector<float> > LBG(vector<vector<float> > in, int M, vector<float> &dist){
	// centroid of universe
	vector<float> centroid(12, 0);
	vector<float> row(12, 0);
	for(int i = 0 ; i < in.size() ; i++){
		transform( centroid.begin(), centroid.end(), in[i].begin(), centroid.begin(), plus<float>());
	}
	for(int i = 0 ; i < 12 ; i++){
		centroid[i] /= (float)in.size();
	}
	vector<vector<float> > c1,c2;
	c2.push_back(centroid);
	int m = 2;
	while(m <= M){
		float eps = 0.03;
		c1.clear();
		for(int i = 0 ; i < c2.size() ; i++){
			vector<float> row_minus(12, 0); vector<float> row_plus(12, 0);
			for(int j = 0 ; j < 12 ; j++){
				row_plus[j] = c2[i][j] + eps;
				row_minus[j] = c2[i][j] - eps;
			}
			c1.push_back(row_plus);
			c1.push_back(row_minus);
		}
		c2.clear();
		c2 = K_means(in, c1, dist);
		m *= 2;
	}
	return c2;
}
int main()
{
    vector<vector<float> > in = getuniverse();
    vector<float> row(12, 0);
    
    // generate 8 size codebook.
    vector<float> dist_8;
    vector<vector<float> > cbook8 = LBG(in, 8, dist_8);
    ofstream myfile;
    myfile.open ("LBG_codebook_8.csv");
	for(int i = 0 ; i < cbook8.size() ; i++){
		for(int j = 0 ; j < 12 ; j++){
			//myfile << cbook8[i][j] << "\t\t";
			myfile << cbook8[i][j] << ",";
		}
		myfile << endl;
	}
	myfile.close();
	myfile.open("LBG_distortion_8.csv");
	for(int i = 0 ; i < dist_8.size() ; i++){
		myfile << dist_8[i] << ",";
	}
	myfile << endl;
	myfile.close();

	// generate 32 size codebook.
	vector<float> dist_32;
    vector<vector<float> > cbook32 = LBG(in, 32, dist_32);
    ofstream myfile2;
    myfile2.open ("LBG_codebook_32.csv");
	for(int i = 0 ; i < cbook32.size() ; i++){
		for(int j = 0 ; j < 12 ; j++){
			//myfile2 << cbook32[i][j] << "\t\t";
			myfile2 << cbook32[i][j] << ",";
		}
		myfile2 << endl;
	}
	myfile2.close();
	myfile2.open("LBG_distortion_32.csv");
	for(int i = 0 ; i < dist_32.size() ; i++){
		myfile2 << dist_32[i] << ",";
	}
	myfile2 << endl;
	myfile2.close();
    
    return 0;
} 