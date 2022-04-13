//Rishabh Rathi
//18HS20027
//Big M method is used

#include <iostream>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <iomanip>
#include <cmath>

using namespace std;

#define vd vector<double>
#define pb push_back

//to add appropriate sign to the objective function coefficients
void construct_obj(vd &z,int min_max)
{
	if(min_max==1)return;
	for(int i=0;i<z.size();i++)
	{
		if(z[i]!=0)
		z[i]*=-1;
	}
	return;
}

//adding all the slack/artificial vairables to the constraints
void construct_con(vector<int>&art, vector<int>&basis,vector<vd>&a,vector<int>&con, int n)
{
	for(int i=0;i<n;i++)
	{
		if(con[i]==0||con[i]==1)
		{
			a[i].pb(1);
			if(con[i]) art.pb(a[i].size());
			for(int j=0;j<n;j++)
			{
				if(j!=i)
				{
					a[j].pb(0);
				}
			}
		}
		else
		{
			a[i].pb(-1);
			a[i].pb(1);
			art.pb(a[i].size());
			for(int j=0;j<n;j++)
			{
				if(j!=i)
				{
					a[j].pb(0);
					a[j].pb(0);
				}
			}
		}
		basis[i]=a[i].size();
	}
	return;
}

//printing a simplex iteration
void print(vector<vd> &table, vd &vm)
{
	for(int i=0;i<table[0].size();i++)
	{
		if(i==0)cout<<' '<<'\t';
		else if(i==table[0].size()-1) cout<<"RHS";
		else cout<<'x'<<i<<'\t';
	}
	cout<<endl;
	cout<<'\t';
	for(int i=1;i<table[1].size();i++)
	{
		if(vm[i]==0)cout<<'\t';
		else cout<<vm[i]<<"M\t";
	}
	cout<<endl;
	for(int i=0;i<table.size();i++)
	{
		for(int j=0;j<table[1].size();j++)
		{
			if(i==0&&j==0) cout<<'Z'<<'\t';
			else if(i>0&&j==0) cout<<fixed<<setprecision(0)<<'x'<<table[i][0]<<'\t';
			else cout<<fixed<<setprecision(2)<<table[i][j]<<'\t';
		}
		cout<<endl;
	}
	return;
}

//checking for optimality or infeasiblility
int is_solved(vd &vm, vector<vd>&v, vector<int>&art)
{
	if(*min_element(vm.begin()+1,vm.end()-1)<-0.001) return 0;
	double mini=DBL_MAX;
	for(int i=1;i<v[0].size()-1;i++)
	{
		if(vm[i]==0)
		{
			mini=min(mini,v[0][i]);
		}
	}
	if(mini<-0.001)return 0;
	else
	{
		for(int i=0;i<art.size();i++)
		{
			for(int j=1;j<v.size();j++)
			{
				if(art[i]==v[j][0]&&v[j][v[j].size()-1]>0)
				{
					return 1;
				}
			}
		}
	}
	return 2;
}

//checking for unboundedness
int is_unbounded(vd &vm, vector<vd> &v)
{
	for(int i=1;i<v[0].size()-1;i++)
	{
		if((vm[i]<0)||(vm[i]==0&&v[0][i]<0))
		{
			int ch=0;
			for(int j=1;j<v.size();j++)
			{
				if(v[j][i]>0)
				{
					ch=1;
					break;
				}
			}
			if(ch==0)
			{
				return 1;
			}
		}
	}
	return 0;
}

//making changes to the original simplex table so that the coefficients of
//the basic variables are made 0
void resolve(vd &M, vector<vd> &table)
{
	for(int i=1;i<M.size()-1;i++)
	{
		if(M[i])
		{
			for(int j=1;j<table.size();j++)
			{
				if(table[j][0]==i)
				{
					for(int k=1;k<table[j].size();k++)
					{
						M[k]+=-1*table[j][k];
					}
					break;
				}
			}
		}
	}
	return;
}

//solving a simplex iteration by determining entering and leaving variables
void solve_simplex(vector<vd>&M, vd vm, vector<vector<vd> > &table, vector<vd>v)
{
	int entering_id;
	double min_m=*min_element(vm.begin()+1,vm.end()-1);
	if(min_m<-0.001)
		entering_id=find(vm.begin()+1,vm.end()-1,min_m)-vm.begin();
	else
	{
		double mini=DBL_MAX;
		for(int i=1;i<v[0].size()-1;i++)
		{
			if(vm[i]==0)
			{
				if(v[0][i]<mini)
				{
					mini=v[0][i];
					entering_id=i;
				}
			}
		}
	}
	vector<double>ratio(v.size()-1);
	for(int i=0;i<v.size()-1;i++)
	{
		if(v[i+1][entering_id]<=0||*(v[i+1].end()-1)<=0)
		ratio[i]=DBL_MAX;
		else
		ratio[i]=(*(v[i+1].end()-1))/v[i+1][entering_id];
		
	}
	int leaving_id = (find(ratio.begin(),ratio.end(),*min_element(ratio.begin(),ratio.end()))-ratio.begin())+1;
	v[leaving_id][0]=entering_id;
	double div=v[leaving_id][entering_id];
	for(int i=1;i<v[leaving_id].size();i++)
	{
		v[leaving_id][i]/=div;
	}
	double mul=v[0][entering_id];
	double mulm=vm[entering_id];
	for(int j=1;j<v[0].size();j++)
	{
		v[0][j]+=((double)round(v[leaving_id][j]*mul*(-1)*10000))/10000.0;
		vm[j]+=((double)round(v[leaving_id][j]*mulm*(-1)*10000))/10000.0;
	}
	for(int i=1;i<v.size();i++)
	{
		if(i!=leaving_id)
		{	
			mul=v[i][entering_id];
			for(int j=1;j<v[i].size();j++)
			{
				v[i][j]+=((double)round(v[leaving_id][j]*mul*(-1)*10000))/10000.0;
			}
		}
	}
	table.pb(v);
	M.pb(vm);
	return;
}

//printing non basic variables and the evaluation in an iteration
void non_basic(vector<vector<vd> > &table, int iter, int min_max)
{
	cout<<"Non-basic variables :"<<endl;
	vector<int>v;
	for(int i=1;i<table[iter].size();i++)
		v.pb(table[iter][i][0]);
	sort(v.begin(),v.end());
	int x=0;
	for(int i=1;i<table[iter][0].size()-1;i++)
	{
		if(i==v[x])
		{
			x++;
			continue;
		}
		else
		{
			cout<<'x'<<i<<'\t';
		}
	}
	cout<<endl;
	cout<<"Evaluation at this iteration : "<<min_max*(*(table[iter][0].end()-1))<<endl;
	return;
}

//printing basic variables and their ratios in an iteration
void basic(vector<vd> &v)
{
	cout<<"Basic variables with their ratios :"<<endl;
	for(int i=1;i<v.size();i++)
		cout<<fixed<<setprecision(0)<<'x'<<v[i][0]<<'\t';
	cout<<endl;
	int entering_id=find(v[0].begin()+1,v[0].end()-1,*min_element(v[0].begin()+1,v[0].end()-1))-v[0].begin();
	vector<double>ratio(v.size()-1);
	for(int i=0;i<v.size()-1;i++)
	{
		if(v[i+1][entering_id]<=0)
		ratio[i]=DBL_MAX;
		else
		ratio[i]=(*(v[i+1].end()-1))/v[i+1][entering_id];
	}
	for(auto x :ratio)
	{
		if(x==DBL_MAX)
		{
			cout<<"NA\t";
		}
		else cout<<fixed<<setprecision(2)<<x<<'\t';
	}
	cout<<endl;
	return;
}

//using gauss seidel method to solve a set of equations
vd solve_gauss_seidel(vector<vd>&V,vd &b,vector<int>index)
{
	vector<vd>A;
	if(V.size()!=V[0].size())
	{
		for(int i=0;i<V.size();i++)
		{
			vd v;
			for(int j=0;j<V[i].size();j++)
			{
				if(index[j])
				{
					v.push_back(V[i][j]);
				}
			}
			A.push_back(v);
		}
	}
	else A=V;
	double err=0.00001;
	int ch=0;
	int iter=100000;
	vd x1(b.size(),0);
	vd x2(b.size(),0);
	while(--iter)
	{
		for(int i=0;i<A.size();i++)
		{
			x2[i]=b[i]/A[i][i];
			for(int j=0;j<b.size();j++)
			{
				if(i!=j)
				{
					x2[i]=x2[i]-((A[i][j]/A[i][i])*x1[j]);
					x1[i]=x2[i];
				}
			}
		}
		for(int i=0;i<A.size();i++)
		{
			double dot=0.0;
			for(int c=0;c<x2.size();c++)
			dot+=x2[c]*A[i][c];
			if(abs(dot-b[i])>err)
			{
				ch=1;
				break;
			}
		}
		if(ch)
		{
			ch=0;
			continue;
		}
		else
		{
			return x2;
		}
	}
	return x2;
}

//calculating factorial using dynamic programming
long long int fact(vector<long long int>&dp,int n)
{
	if(dp[n]!=0)return dp[n];
	for(int i=3;i<=n;i++)
	{
		dp[i]=dp[i-1]*i;
	}
	return dp[n];
}

//calculating number of combinations
long long int comb(int n,int m,vector<long long int>&dp)
{
	return fact(dp,n)/(fact(dp,m)*fact(dp,n-m));
}

void print_bfs(vd x,vector<int>&index,long long int i)
{
	if(*min_element(x.begin(),x.end())<0)return;
	for(int i=0;i<x.size();i++)
	{
		if(isnan(x[i]))return;
	}
	cout<<"Solution "<<i+1<<endl;
	int xp=0;
	for(int i=0;i<index.size();i++)
	{
		if(index[i])
		cout<<'x'<<i+1<<" : "<<x[xp++]<<'\t';
		else cout<<'x'<<i+1<<" : "<<'0'<<'\t';
	}
	cout<<endl;
	return;
}

//finding and printing bfs for a given problem
void bfs(vector<vd>a, vd rhs)
{
	int m=a[0].size();
	int n=a.size();
	vector<int>index(m,1);
	vector<long long int> dp(m,0);
	dp[0]=1;
	dp[1]=1;
	dp[2]=2;
	long long int mCn=comb(m,n,dp);
	for(int i=0;i<m-n;i++)index[i]=0;
	cout<<"Basic Feasible Solutions : "<<endl;
	for(long long int i=0;i<mCn;i++)
	{
		vd x=solve_gauss_seidel(a,rhs,index);
		print_bfs(x,index,i);
		next_permutation(index.begin(),index.end());
	}
	return;
}

int main()
{
	cout<<fixed<<setprecision(2);
	cout<<"Enter type of problem :\nMinimization : -1\nMaximization : 1"<<endl;
	int min_max;
	cin>>min_max;
	int m;
	cout<<"Enter total no. of variables :"<<endl;
	cin>>m;
	cout<<"Enter objective function coeffecients (variable wise) :"<<endl;
	vd z(m);
	vd zm(m,0);
	for(int i=0;i<m;i++)
	{
		cin>>z[i];
	}
	construct_obj(z,min_max);
	cout<<"Enter total no. of constraints :"<<endl;
	int n;
	cin>>n;
	cout<<"Enter the LHS of the constraints :"<<endl;
	vector<vd>a(n,vd(m));
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<m;j++)
		{
			cin>>a[i][j];
		}
	}
	cout<<"Enter type of constraints :"<<endl;
	cout<<"<= : 0\n = : 1\n>= : 2"<<endl;
	vector<int>con(n);
	for(int i=0;i<n;i++)
	{
		cin>>con[i];
	}
	cout<<"Enter RHS of each constraint :"<<endl;
	vd rhs(n);
	for(int i=0;i<n;i++)
	{
		cin>>rhs[i];
		if(rhs[i]<0)
		{
			rhs[i]=-rhs[i];
			if(con[i]==0)con[i]=2;
			else if(con[i]==2)con[i]=0;
			for(int k=0;k<a[i].size();k++)
				a[i][k]=-a[i][k];
		}
	}
	vector<int> basis(n);
	vector<int>art;
	construct_con(art,basis,a,con,n);
	while(z.size()<a[0].size())
	{
		z.pb(0);
		zm.pb(0);
	}
	for(auto x:art)
		zm[x-1]=1;
	vector<vector<vd> >simplex_iterations;
	vector<vd>M;
	vector<vd> v(n+1,vd(a[0].size()+2,0));
	vd vm(a[0].size()+2,0);
	for(int i=1;i<=a[0].size();i++)
	{
		if(z[i-1])
			v[0][i]=-1*z[i-1];
		else v[0][i]=0;
		if(zm[i-1])
			vm[i]=1;
		else vm[i]=0;
	}
	for(int i=1;i<n+1;i++)
	{
		for(int j=0;j<a[i-1].size()+2;j++)
		{
			if(j==0)
			{
				v[i][j]=basis[i-1];
			}
			else if(j==a[i-1].size()+1)
			{
				v[i][j]=rhs[i-1];
			}
			else
			{
				v[i][j]=a[i-1][j-1];
			}
			
		}
	}
	simplex_iterations.pb(v);
	M.pb(vm);
	resolve(M[0],simplex_iterations[0]);
	int sol,ub=0;
	//solving the problem using simplex method
	while(1)
	{
		sol = is_solved(M[M.size()-1],simplex_iterations[simplex_iterations.size()-1],art);
		if(sol)break;
		ub=is_unbounded(M[M.size()-1],simplex_iterations[simplex_iterations.size()-1]);
		if(ub==1)break;
		solve_simplex(M, M[M.size()-1], simplex_iterations,simplex_iterations[simplex_iterations.size()-1]);
	}
	//menu driven answers
	while(1)
	{
		cout<<"MENU"<<endl;
		cout<<"1. List of all basic feasible solutions"<<endl;
		cout<<"2. Number of iterations to solve the problem"<<endl;
		cout<<"3. List of all non basic variables with net evaluation on the ith iteration"<<endl;
		cout<<"4. List of basic variables with minimum ratios on the ith iteration"<<endl;
		cout<<"5. Simplex table on ith iteration"<<endl;
		cout<<"6. Solution"<<endl;
		cout<<"7. EXIT"<<endl;
		cout<<"Enter your choice : ";
		int ch;
		cin>>ch;
		if(ch==1)
		{
			bfs(a,rhs);
		}
		else if(ch==2) cout<<"Number of iterations : "<<simplex_iterations.size()-1<<endl;
		else if(ch==3)
		{
			int i;
			cout<<"Enter iteration number : ";
			cin>>i;
			non_basic(simplex_iterations,i,min_max);
		}
		else if(ch==4)
		{
			int i;
			cout<<"Enter iteration number : ";
			cin>>i;
			basic(simplex_iterations[i]);
		}
		else if(ch==5)
		{
			int i;
			cout<<"Enter iteration number : ";
			cin>>i;
			print(simplex_iterations[i],M[i]);
		}
		else if(ch==6)
		{
			if(ub) cout<<"Unbounded"<<endl;
			else if(sol==1) cout<<"Infeasible"<<endl;
			else cout<<"Optimal solution reached : "<<min_max*(*(simplex_iterations[simplex_iterations.size()-1][0].end()-1))<<endl;
		}
		else break;
	}
	return 0;
}