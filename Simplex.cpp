//Rishabh Rathi
//18HS20027

#include <iostream>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <iomanip>
#include <cmath>

using namespace std;

#define vd vector<double>
#define pb push_back

void construct_obj(vd &z,int min_max)
{
	if(min_max)return;
	for(int i=0;i<z.size();i++)
	{
		if(z[i]!=0)
		z[i]*=-1;
	}
	return;
}

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

void print(vector<vd> &table)
{
	for(int i=0;i<table[0].size();i++)
	{
		if(i==0)cout<<' '<<'\t';
		else if(i==table[0].size()-1) cout<<"RHS";
		else cout<<'x'<<i<<'\t';
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

int is_solved(vector<vd>&v, vector<int>&art)
{
	if(*min_element(v[0].begin()+1,v[0].end()-1)<0) return 0;
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

int is_unbounded(vector<vd>&v)
{
	for(int i=1;i<v[0].size()-1;i++)
	{
		if(v[0][i]<0)
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

void solve_simplex(vector<vector<vd> > &table, vector<vd>v)
{
	int entering_id=find(v[0].begin()+1,v[0].end()-1,*min_element(v[0].begin()+1,v[0].end()-1))-v[0].begin();
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
	for(int i=0;i<v.size();i++)
	{
		if(i!=leaving_id)
		{	
			double mul=v[i][entering_id];
			for(int j=1;j<v[i].size();j++)
			{
				v[i][j]+=(v[leaving_id][j]*mul*(-1));
			}
		}
	}
	table.pb(v);
	return;
}

void non_basic(vector<vector<vd> > &table, int iter)
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
	cout<<"Evaluation at this iteration : "<<(*(table[iter][0].end()-1))<<endl;
	return;
}

void basic(vector<vd> &v)
{
	cout<<"Basic variables :"<<endl;
	for(int i=1;i<v.size();i++)
		cout<<'x'<<v[i][0]<<'\t';
	cout<<"Minimum Ratio : ";
	int entering_id=find(v[0].begin()+1,v[0].end()-1,*min_element(v[0].begin()+1,v[0].end()-1))-v[0].begin();
	vector<double>ratio(v.size()-1);
	for(int i=0;i<v.size()-1;i++)
	{
		if(v[i+1][entering_id]<=0)
		ratio[i]=DBL_MAX;
		else
		ratio[i]=(*(v[i+1].end()-1))/v[i+1][entering_id];
	}
	cout<<*min_element(ratio.begin(),ratio.end())<<endl;
	return;
}

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

int fact(vector<int>&dp,int n)
{
	if(dp[n]!=0)return dp[n];
	for(int i=3;i<=n;i++)
	{
		dp[i]=dp[i-1]*i;
	}
	return dp[n];
}

int comb(int n,int m,vector<int>&dp)
{
	return fact(dp,n)/(fact(dp,m)*fact(dp,n-m));
}

void print_bfs(vd x,vector<int>&index,int i)
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

void bfs(vector<vd>a, vd rhs)
{
	int m=a[0].size();
	int n=a.size();
	vector<int>index(m,1);
	vector<int> dp(m,0);
	dp[0]=1;
	dp[1]=1;
	dp[2]=2;
	int mCn=comb(m,n,dp);
	for(int i=0;i<m-n;i++)index[i]=0;
	cout<<"Basic Feasible Solutions : "<<endl;
	for(int i=0;i<mCn;i++)
	{
		vd x=solve_gauss_seidel(a,rhs,index);
		print_bfs(x,index,i);
		next_permutation(index.begin(),index.end());
		//x.clear();
	}
	return;
}

int main()
{
	cout<<"Enter type of problem :\nMinimization : 0\nMaximization : 1"<<endl;
	int min_max;
	cin>>min_max;
	int m;
	cout<<"Enter total no. of variables :"<<endl;
	cin>>m;
	cout<<"Enter objective function coeffecients (variable wise) :"<<endl;
	vd z(m);
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
	}
	vector<int> basis(n);
	vector<int>art;
	construct_con(art,basis,a,con,n);
	while(z.size()<a[0].size())
	{
		z.pb(0);
	}
	vector<vector<vd> >simplex_iterations;
	vector<vd> v(n+1,vd(a[0].size()+2,0));
	for(int i=1;i<a[0].size();i++)
	{
		if(z[i-1])
		v[0][i]=-1*z[i-1];
		else v[0][i]=0;
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
	int sol,ub=0;
	while(1)
	{
		sol = is_solved(simplex_iterations[simplex_iterations.size()-1],art);
		if(sol)break;
		ub=is_unbounded(simplex_iterations[simplex_iterations.size()-1]);
		if(ub==1)break;
		solve_simplex(simplex_iterations,simplex_iterations[simplex_iterations.size()-1]);
	}
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
			non_basic(simplex_iterations,i);
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
			print(simplex_iterations[i]);
		}
		else if(ch==6)
		{
			if(ub) cout<<"Unbounded"<<endl;
			else if(sol==1) cout<<"Infeasible"<<endl;
			else cout<<"Optimal solution reached : "<<*(simplex_iterations[simplex_iterations.size()-1][0].end()-1)<<endl;
		}
		else break;
	}
	return 0;
}