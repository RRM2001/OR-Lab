//Rishabh Rathi
//18HS20027
//Big M method is used for each subproblem

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

void so(int m, int n, vd z, vd zm, vector<vd>a, vector<int>con, vd rhs)
{
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
	if(sol==1)return;
	if(ub==1)return;
	int b=0;
	for(int i=1;i<simplex_iterations[simplex_iterations.size()-1].size();i++)
	{
		if(*(simplex_iterations[simplex_iterations.size()-1][i].end()-1)!=int(*(simplex_iterations[simplex_iterations.size()-1][i].end()-1)))
		{
			b=i;
			break;
		}
	}
	if(b>0)
	{
		vd v;
		for(int i=1;i<m+1;i++)
		{
			if(b!=i)v.pb(0);
			else v.pb(1);
		}
		for(auto x:a)
		{
			x.erase(x.begin()+m+1,x.end());
		}
		a.pb(v);
		con.pb(1);
		rhs.pb(floor(*(simplex_iterations[simplex_iterations.size()-1][b].end()-1)));
		so(m,n+1,z,zm,a,con,rhs);
		for(auto x:a)
		{
			x.erase(x.begin()+m+1,x.end());
		}
		con.erase(con.end()-1,con.end());
		rhs[m-1]=ceil(*(simplex_iterations[simplex_iterations.size()-1][b].end()-1));
		a.pb(v);
		con.pb(2);
		so(m,n+1,z,zm,a,con,rhs);
	}
	else
	{
		vector<int>x(m,0);
		for(int i=1;i<simplex_iterations[simplex_iterations.size()-1].size();i++)
		{
			x[simplex_iterations[simplex_iterations.size()-1][i][0]-1]=*(simplex_iterations[simplex_iterations.size()-1][i].end()-1);
		}
		for(int i=0;i<x.size();i++)
		{
			cout<<"x"<<i+1<<" = "<<x[i]<<endl;
		}
		cout<<"Z = "<<*(simplex_iterations[simplex_iterations.size()-1][0].end()-1)<<endl;
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
	so(m,n,z,zm,a,con,rhs);
	return 0;
}