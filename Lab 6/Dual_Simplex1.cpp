//Rishabh Rathi
//18HS20027
//Dual Simplex

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
	if(min_max==1)return;
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

int is_solved(vector<vd>&v)
{
	for(int i=1;i<v.size();i++)
	{
		if(*(v[i].end()-1)<0)return 0;
	}
	return 2;
}
int is_solved_simplex(vector<vd>&v, vector<int>&art)
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


int solve_simplex(vector<vector<vd> > &table, vector<vd>v,int xx, int iter)
{
	int leaving_id,entering_id;
	if(xx==1)
	{
	double min=DBL_MAX;
	for(int i=1;i<v.size();i++)
	{
		if(*(v[i].end()-1)<min)
		{
			min=*(v[i].end()-1);
			leaving_id=i;
		}
	}
	vector<double>ratio(v[0].size()-2);
	for(int i=1;i<v[0].size()-1;i++)
	{
		if(v[leaving_id][i]>=0||v[0][i]==0)
		ratio[i-1]=DBL_MAX;
		else
		ratio[i-1]=abs(v[0][i]/v[leaving_id][i]);
		
	}
	entering_id = (find(ratio.begin(),ratio.end(),*min_element(ratio.begin(),ratio.end()))-ratio.begin())+1;
	}
	else
	{
	entering_id=find(v[0].begin()+1,v[0].end()-1,*min_element(v[0].begin()+1,v[0].end()-1))-v[0].begin();
	vector<double>ratio(v.size()-1);
	for(int i=0;i<v.size()-1;i++)
	{
		if(v[i+1][entering_id]<=0||*(v[i+1].end()-1)<=0)
		ratio[i]=DBL_MAX;
		else
		ratio[i]=(*(v[i+1].end()-1))/v[i+1][entering_id];
		
	}
	double d=*min_element(ratio.begin(),ratio.end());
	if(d==DBL_MAX)
		return 1;
	leaving_id = (find(ratio.begin(),ratio.end(),*min_element(ratio.begin(),ratio.end()))-ratio.begin())+1;
	}
	cout<<"Iteration "<<iter<<endl<<endl;
	cout<<"Entering variable : x"<<entering_id<<endl;
	cout<<fixed<<setprecision(0)<<"Leaving variable : x"<<v[leaving_id][0]<<endl;
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
	cout<<"Objective value : "<<*(v[0].end()-1)<<endl<<endl;
	//print(v);
	table.pb(v);
	return 0;
}

int main()
{
	cout<<"Enter type of problem :\nMinimization : 1\nMaximization : -1"<<endl;
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
		if(rhs[i]>0)
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
	}
	vector<vector<vd> >simplex_iterations;
	vector<vd> v(n+1,vd(a[0].size()+2,0));
	for(int i=1;i<a[0].size();i++)
	{
		if(z[i-1])
		v[0][i]=1*z[i-1];
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
	cout<<"Initial table :"<<endl;
	print(simplex_iterations[0]);
	cout<<endl<<endl<<endl;
	int sol=0,sol2,ub=0,iter=1;
	while(iter++)
	{
		int s=0;
		if(sol==0)
		sol = is_solved(simplex_iterations[simplex_iterations.size()-1]);
		if(sol)
		{
			sol2=is_solved_simplex(simplex_iterations[simplex_iterations.size()-1],art);
			if(sol2)
			break;
			else
			s=solve_simplex(simplex_iterations,simplex_iterations[simplex_iterations.size()-1],2,iter-1);
			if(s)break;
			continue;
		}
		//ub=is_unbounded(simplex_iterations[simplex_iterations.size()-1]);
		//if(ub==1)break;
		s=solve_simplex(simplex_iterations,simplex_iterations[simplex_iterations.size()-1],1,iter-1);
	}
	cout<<endl<<endl<<endl;
	cout<<"FINAL SOLUTION"<<endl;
	cout<<"Optimal value : "<<*(simplex_iterations[simplex_iterations.size()-1][0].end()-1);
	return 0;
}