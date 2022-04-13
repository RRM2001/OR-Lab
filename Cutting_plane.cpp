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
#define vvd vector<vector<double> >

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
		solve_simplex(simplex_iterations,simplex_iterations[simplex_iterations.size()-1]);
	}
	while(1)
	{
		if(int_sol(simplex_iterations[simplex_iterations.size()-1])) break;
		solve_int(simplex_iterations[simplex_iterations.size()-1], simplex_iterations)
	}
	while(1)
	{
		cout<<"MENU"<<endl;
		cout<<"1. Equation of cutting plane in each iteration"<<endl;
		cout<<"2. Incoming variable in each iteration"<<endl;
		cout<<"3. Outgoing variable in each iteration"<<endl;
		cout<<"4. Print the pivot element"<<endl;
		cout<<"5. EXIT"<<endl;
		cout<<"Enter your choice : ";
		int ch;
		cin>>ch;
		if(ch==1)
		else if(ch==2)
		else if(ch==3)
		else if(ch==4)
		else break;
	}
	return 0;
}