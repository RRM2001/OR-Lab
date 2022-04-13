//Rishabh Rathi
//18HS20027
//2-phase method

#include <iostream>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <iomanip>
#include <cmath>

using namespace std;

#define vd vector<double>
#define pb push_back

//constructing objective function
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

//adding slack/surplus/artificial variables
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

//printing a simplex table
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

void resolve_p1(vector<vd> &table)
{
	for(int i=1;i<table[0].size()-1;i++)
	{
		if(table[0][i]==1)
		{
			for(int j=1;j<table.size();j++)
			{
				if(table[j][0]==i)
				{
					for(int k=1;k<table[j].size();k++)
					{
						table[0][k]-=table[j][k];
					}
					break;
				}
			}
		}
	}
	return;
}

//solving via simplex method
void solve_simplex(vector<vector<vd> > &table, vector<vd>v,int *art_basis,vector<int>&art)
{
	int entering_id=find(v[0].begin()+1,v[0].end()-1,*min_element(v[0].begin()+1,v[0].end()-1))-v[0].begin();
	int leaving_id=0;
	if(*art_basis==0)
	{
		vector<double>ratio(v.size()-1);
		for(int i=0;i<v.size()-1;i++)
		{
			if((v[i+1][entering_id]<=0&&*(v[i+1].end()-1)>=0)||(v[i+1][entering_id]>=0&&*(v[i+1].end()-1)<=0))
			ratio[i]=DBL_MAX;
			else
			ratio[i]=(*(v[i+1].end()-1))/v[i+1][entering_id];
		}
		leaving_id = (find(ratio.begin(),ratio.end(),*min_element(ratio.begin(),ratio.end()))-ratio.begin())+1;
	}
	else
	{
		for(int i=1;i<v.size();i++)
		{
			for(auto x:art)
			{
				if(x==v[i][0])
				{
					if(v[i][entering_id]<0)
					{
						leaving_id=i;
						break;
					}
				}
			}
			if(leaving_id!=0)break;
		}
	}
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
	*art_basis=0;
	for(int i=1;i<v.size();i++)
	{
		for(auto x:art)
		{
			if(x==v[i][0]) *art_basis=1;
		}
	}
	return;
}

//checking for solution/feasibility in phase 1
int is_solved_p1(vector<vd> &v, vector<int> &art)
{
	int ch=0,c=0;
	for(int i=1;i<v.size();i++)
	{
		if(find(art.begin(),art.end(),v[i][0])!=art.end())
		{
			if(*(v[i].end()-1)!=0) ch=2;//infeasible
		}
	}
	if(*(v[0].end()-1)==0&&ch==2)return 2;
	if(*(v[0].end()-1)==0&&ch==0)return 1;
	for(int i=1;i<v[0].size()-1;i++)
	{
		if(v[0][i]<0)c=1;
	}
	if(c==0&&ch==0)return 1;
	if(c==0&&ch==2)return 2;
	if(c==1)return 0;
}

//checking for solution in phase 2
int is_solved_p2(vector<vd>&v, vector<int>&art)
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

int resolve(vector<vd>&p2,vector<int> &art, vd &zp2, vector<int> &var)
{
	vector<int> drop;
	int ch=0,ar=0;
	for(auto x:art)
	{
		ch=0;
		for(int i=1;i<p2.size();i++)
		{
			if(x==p2[i][0])
			{
				ch=1;
				ar=1;
			}
		}
		if(ch==0)drop.pb(x);
	}
	//dropping artificial variable columns
	if(drop.size())
	{
		if(drop.size()>1)sort(drop.begin(),drop.end());
		while(drop.size()>=1)
		{
			var.erase(var.begin()+drop[0],var.begin()+drop[0]+1);
			for(int i=0;i<p2.size();i++)
			{
				p2[i].erase(p2[i].begin()+drop[0],p2[i].begin()+drop[0]+1);
			}
			for(int i=0;i<drop.size();i++)
				drop[i]--;
			if(drop.size()==1) break;
			drop.erase(drop.begin(),drop.begin()+1);
		}
	}
	//replacing coefficients
	for(int i=1;i<p2[0].size()-1;i++)
		p2[0][i]=-zp2[var[i]-1];
	//calculating indicator row values
	for(int i=1;i<p2.size();i++)
	{
		int id=find(var.begin(),var.end(),p2[i][0])-var.begin();
		if(p2[0][id]!=0)
		{
			double div=p2[i][id];
			for(int j=1;j<p2[i].size();j++)
			{
				p2[i][j]/=div;
			}
			double mul=-p2[0][id];
			for(int j=1;j<p2[i].size();j++)
			{
				p2[0][j]+=(mul*p2[i][j]);
			}
		}
	}
	return ar;
}

int main()
{
	cout<<"Enter type of problem :\nMinimization : -1\nMaximization : 1"<<endl;
	int min_max;
	cin>>min_max;
	int m;
	cout<<"Enter total no. of variables :"<<endl;
	cin>>m;
	cout<<"Enter objective function coeffecients (variable wise) :"<<endl;
	vd zp2(m);
	for(int i=0;i<m;i++)
	{
		cin>>zp2[i];
	}
	construct_obj(zp2,min_max);
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
	vector<int> variables(a[0].size()+1,0);
	for(int i=0;i<variables.size();i++)
		variables[i]=i;
	while(zp2.size()<a[0].size())
	{
		zp2.pb(0);
	}
	vd zp1(zp2.size(),0);
	for(auto x:art)
	{
		zp1[x-1]=1;
	}
	vector<vector<vd> >simplex_phase1;
	vector<vd> v(n+1,vd(a[0].size()+2,0));
	for(int i=1;i<v[0].size()-1;i++)
	{
		v[0][i]=zp1[i-1];
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
	simplex_phase1.pb(v);
	resolve_p1(simplex_phase1[0]);
	int sol,ub=0;
	while(1)
	{
		if(art.size()==0)break;
		sol=is_solved_p1(simplex_phase1[simplex_phase1.size()-1],art);
		if(sol==2)
		{
			cout<<"The problem is infeasible.";
			return 0;
		}
		if(sol==1)
		{
			break;
		}
		int x=0;
		solve_simplex(simplex_phase1,simplex_phase1[simplex_phase1.size()-1],&x,art);
	}
	sol=0;
	vector<vector<vd> >simplex_phase2;
	simplex_phase2.pb(simplex_phase1[simplex_phase1.size()-1]);
	int art_basis=resolve(simplex_phase2[0],art,zp2,variables);
	while(1)
	{
		sol = is_solved_p2(simplex_phase2[simplex_phase2.size()-1],art);
		if(sol)break;
		ub=is_unbounded(simplex_phase2[simplex_phase2.size()-1]);
		if(ub==1)break;
		solve_simplex(simplex_phase2,simplex_phase2[simplex_phase2.size()-1],&art_basis,art);
	}
	//menu driven answers
	while(1)
	{
		cout<<"MENU"<<endl;
		cout<<"1. Initial table for phase 1"<<endl;
		cout<<"2. Initial table for phase 2"<<endl;
		cout<<"3. Solution"<<endl;
		cout<<"4. EXIT"<<endl;
		cout<<"Enter your choice : ";
		int ch;
		cin>>ch;
		if(ch==1)
		{
			print(simplex_phase1[0]);
		}
		else if(ch==2)
		{
			print(simplex_phase2[0]);
		}
		else if(ch==3)
		{
			if(ub==1)
				cout<<"The problem is unbounded."<<endl;
			else
			{
				cout<<"Optimal value of objective function : "<<*(simplex_phase2[simplex_phase2.size()-1][0].end()-1)*min_max<<endl;
				for(int i=1;i<simplex_phase2[simplex_phase2.size()-1].size();i++)
				{
					cout<<fixed<<setprecision(0)<<"x"<<simplex_phase2[simplex_phase2.size()-1][i][0];
					cout<<fixed<<setprecision(2)<<" = "<<*(simplex_phase2[simplex_phase2.size()-1][i].end()-1)<<endl;
				}
			}
		}
		else break;
	}
	return 0;
}