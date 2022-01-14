//Rishabh Rathi
//18HS20027
//gauss seidel method for n equations, n variables

#include<iostream>
#include<vector>
using namespace std;
#define vd vector<double>
#define vvd vector<vector<double>>

double dot_prod(vd arr1, vd arr2){
	int n = arr1.size();
	double p=0;
	for(int i=0;i<n;i++)
	{
          p+=(arr1[i]*arr2[i]);
        }
	return p;
}

void pritntvecd(vd arr){
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

vd solve_gauss_siedel(vvd mat){
	vd temp(mat.size(), 0);
	int n = mat.size();

	int num_iterations=10;
	for(int i=0;i<num_iterations;i++){
		for(int j=0;j<n;j++){
			double dot = dot_prod(temp_arr, mat[j]);
			dot-=(mat[j][j]*temp[j]);
			temp_arr[j] = (mat[j][n] - dot)/mat[j][j];
		}
		pritntvecd(temp);
	}	

	return temp;
}

void printans(vd vec){
	int n = vec.size();
	for(int i=0;i<n;i++) cout<<vec[i]<<" ";
	cout<<endl;
}

int main(){
	int n;
	cout<<"(n): "; cin>>n;
	vector<vector<double>> ma(n , vector<double> (n+1, 0));

	for(int i=0;i<n;i++){
		cout<<"Enter equation "<<i<<" coefficients: "<<endl;
		for(int j=0;j<n+1;j++){
			cin>>ma[i][j];
		}
	}

	vd sol = solve_gauss_siedel(ma);
	cout<<"Basic solutions: ";
	printans(sol);
}
