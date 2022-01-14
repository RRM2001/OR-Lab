
#include<iostream>
#include<vector>
#include<algorithm>

using namespace std;
#define vd vector<double>
#define vvd vector<vector<double>>
#define vi vector<int>
#define num_iterations 50
#define pb push_back



//Find dot product of 2 n dimensional vectors
double dot_prod(vd arr1, vd arr2){
	int n = arr1.size();
	double p=0;
	for(int i=0;i<n;i++)
		p+=(arr1[i]*arr2[i]);
	return p;
}

//Solve a system of n equations, n variables using Gauss-Siedel method using 0 initialisation
vd solve_gauss_siedel(vvd mat){
	vd temp(mat.size(), 0);
	int n = mat.size();

	for(int i=0;i<num_iterations;i++){
		for(int j=0;j<n;j++){
			double dot = dot_prod(temp, mat[j]);
			dot-=(mat[j][j]*temp[j]);
			temp[j] = (mat[j][n] - dot)/mat[j][j];
		}
	}	

	return temp;
}

//Print the obtained basic solution
void printans(vd vec, vi index_set){
	int n = index_set.size();
	int count=0;
	for(int i=0;i<n;i++){
		if(index_set[i]==0) cout<<"0(non-basic var) ";
		else cout<<vec[count++]<<" ";
	}
	cout<<endl<<endl;
}

//calculate factorial of n
int fact(int n, vi &dp){
	if(n==1) return 1;
	if(n==2) return 2;
	if(dp[n]!=1) return dp[n];

	dp[n] = n*fact(n-1,dp);
	return dp[n];
}

//Calculate n combination m
int nCm(int n, int m){
	vi dp(n+1,1);
	dp[2]=2;

	int temp = fact(n,dp);
	int temp2 = fact(m,dp);
	int temp3 = fact(n-m,dp);

	return (temp/(temp2*temp3));
}

//Construct a matrix(from the input set of equations) to be used for the Gauss-Siedel method
vvd form_mat(vi index_set, vvd mat){
	int m = mat.size();
	int n = (mat[0].size())-1;
	vvd ret;
	for(int i=0;i<m;i++){
		vd temp;
		for(int j=0;j<n;j++){
			if(index_set[j]==0) continue;
			temp.pb(mat[i][j]);
		}
		temp.pb(mat[i][n]);
		ret.pb(temp);
	}

	return ret;
}
double printz(vd &vec, vi &index_set, vd &z)
{
    int n = index_set.size();
	int count=0,cz=0;
	double Z=0;
	for(int i=0;i<n;i++){
		if(index_set[i]==0)
		{
		    cz++;
		    continue;
		}
		else Z+=z[cz]*vec[count++];
		cz++;
	}
	cout<<"Objective Function Value : "<<Z<<endl;
	return Z;
}

int main(){
	int n, m;
	cin>>m;//eqns
	cin>>n;//var
	vector<vector<double>> mat(m , vector<double> (n+1, 0)); 
    vd z(n);
    vd obj;
    for(int i=0;i<n;i++)
    {
        cin>>z[i];
    }
	for(int i=0;i<m;i++){
		for(int j=0;j<n+1;j++){
			cin>>mat[i][j];
		}
	}

	vi index_set(n,1);
	for(int i=0;i<n-m;i++) index_set[i]=0;

	cout<<"Basic solutions: "<<endl;
	vvd mat_temp;

	int num_cases = nCm(n,m);
	for(int i=0;i<num_cases;i++){
		mat_temp = form_mat(index_set,mat);
		
		vd ans = solve_gauss_siedel(mat_temp);
		printans(ans,index_set);
        obj.pb(printz(ans, index_set, z));
		next_permutation(index_set.begin(), index_set.end());
	}
	cout<<endl<<"Optimal Value : "<<*max_element(obj.begin(), obj.end());
}
