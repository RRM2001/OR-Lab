//Rishabh Rathi
//18HS20027

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <numeric>
#include <unordered_set>


#define vi vector<int>
#define vvi vector<vector<int> >
#define pii pair<int,int>
#define ppi pair<pair<int,int>,pair<int,int> >
#define pb push_back
#define mp make_pair
#define ps prob.size()

using namespace std;

bool check_z(vector<ppi >&v, vector<pii >z)
{
	for(auto y:z)
	{
		int f=0;
		for(auto x:v)
		{
			if((x.first.first==y.first && x.second.first==y.first) ||(x.first.second==y.second && x.second.second==y.second))
			{
				f=1;
				break;
			}
		}
		if(f==0)return false;
	}
	return true;
}

bool move(int m,int vsize,vector<ppi>&lines,vector<ppi>v,vvi &prob, vector<pii>&z,char l)
{
	cout<<endl;
	for(auto x:v)
	{
		cout<<x.first.first<<", "<<x.first.second<<" --> "<<x.second.first<<", "<<x.second.second<<endl;
	}
	cout<<endl<<endl;
	if(m==vsize)
	{
		if(check_z(v,z))
		{
			for(auto x:v)
				lines.pb(x);
			return true;
		}
		return false;
	}
	else if(m<vsize)
	{
		if(v.size()==0)
			v.pb(make_pair(make_pair(0,0),make_pair(prob.size()-1,0)));
		else
		{
			v.pb(v[m-1]);
			if(l=='v'&& v[v.size()-1].first.second<prob.size()-1)
			{
				v[v.size()-1].first.second += 1;
				v[v.size()-1].second.second += 1;
			}
			else if(l=='h' && v[v.size()-1].first.first<prob.size()-1)
			{
				v[v.size()-1].first.first += 1;
				v[v.size()-1].second.first += 1;
			}
			else if(l=='v'&& v[v.size()-1].first.second==prob.size()-1)
			{
				v[v.size()-1].first.first = 0;
				v[v.size()-1].first.second = 0;
				v[v.size()-1].second.first = 0;
				v[v.size()-1].second.second = prob.size()-1;
			}
		}
	}
	while(1)
	{
		if(v[v.size()-1].second.first == v[v.size()-1].first.first) l='h';
		else l='v';
		if(move(m+1,vsize,lines,v,prob,z,l))
			return true;
		if(l=='v'&& v[v.size()-1].first.second<prob.size()-1)
		{
			v[v.size()-1].first.second += 1;
			v[v.size()-1].second.second += 1;
		}
		else if(l=='h' && v[v.size()-1].first.first<=prob.size()-1)
		{
			v[v.size()-1].first.first += 1;
			v[v.size()-1].second.first += 1;
		}
		else if(l=='v'&& v[v.size()-1].first.second==prob.size()-1)
		{
			v[v.size()-1].first.first = 0;
			v[v.size()-1].first.second = 0;
			v[v.size()-1].second.first = 0;
			v[v.size()-1].second.second = prob.size()-1;
			l='h';
		}
		if(l=='h')
		{
			if(v[v.size()-1].first.first>prob.size()-vsize+m) return false;
		}
	}
}

void make_lines(vector<ppi>&v,vvi &prob, vector<pii>&z)
{
	for(int i=1;i<=ps;i++)
	{
		v.clear();
		vector<ppi>lines;
		if(move(0,i,v,lines,prob,z,'v')) return;
	}
}

void solve_assg(vvi &prob)
{
	vector<pii>z;
	for(int j=0;j<ps;j++)
	{
		int mi=*min_element(prob[j].begin(),prob[j].end());
		for(int i=0;i<ps;i++)
		{
			prob[j][i]-=mi;
			if(prob[j][i]==0)
				z.pb(make_pair(j,i));
		}
	}
	for(auto x:prob)
	{
		for(auto y:x)
			cout<<y<<' ';
		cout<<endl;
	}
	cout<<endl<<endl;
	for(int i=0;i<ps;i++)
	{
		int mi=INT_MAX;
		for(int j=0;j<ps;j++)
		{
			mi=min(mi,prob[j][i]);
		}
		if(mi!=0)
		{
			for(int j=0;j<ps;j++)
			{
				prob[j][i]-=mi;
				if(prob[j][i]==0)
					z.pb(make_pair(j,i));
			}
		}
	}
	while(1)
	{
		for(auto x:prob)
		{
			for(auto y:x)
				cout<<y<<' ';
			cout<<endl;
		}
		cout<<endl<<endl;
		vector<ppi> lines;
		make_lines(lines,prob,z);
		if(lines.size()==ps)
		{
			return;
		}
		
		for(auto x:lines)
		{
			cout<<x.first.first<<", "<<x.first.second<<" --> "<<x.second.first<<", "<<x.second.second<<endl;
		}
		vvi covered(ps,vi(ps,0));
		for(auto x:lines)
		{
			if(x.first.first==x.second.first)
			{
				for(int i=0;i<ps;i++)
					covered[x.first.first][i]++;
			}
			else
			{
				for(int i=0;i<ps;i++)
					covered[i][x.first.second]++;
			}
		}
		int mi=INT_MAX;
		for(int i=0;i<ps;i++)
		{
			for(int j=0;j<ps;j++)
			{
				if(covered[i][j]==0)
				{
					mi=min(mi,prob[i][j]);
				}
			}
		}
		for(int i=0;i<ps;i++)
		{
			for(int j=0;j<ps;j++)
			{
				if(covered[i][j]==0)
				{
					prob[i][j]-=mi;
				}
				else if(covered[i][j]==2)
				{
					prob[i][j]+=mi;
				}
			}
		}
		z.clear();
		for(int i=0;i<ps;i++)
		{
			for(int j=0;j<ps;j++)
			{
				if(prob[i][j]==0)
				{
					z.pb(mp(i,j));
				}
			}
		}
	}
	return;
}

int allot(vector<pii> &v, vvi &prob, vi &c)
{
	if(v.size()==ps) return 1;
	for(int i=0;i<ps;i++)
	{
		if(prob[v.size()][i]==0 && c[i]==0)
		{
			v.pb(mp(v.size(),i));
			c[i]=1;
			if(allot(v,prob,c))return 1;
			c[i]=0;
			v.erase(v.end()-1,v.end());
		}
	}
	return 0;
}

int main()
{
	int rows,col;
	cout<<"Enter number of rows :";
	cin>>rows;
	cout<<"Enter number of columns :";
	cin>>col;
	cout<<"Enter the problem matrix :"<<endl;
	vvi prob(rows,vi(col));
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<col;j++)
			cin>>prob[i][j];
	}
	while(rows<col)
	{
		vi v(col,0);
		prob.pb(v);
		rows++;
	}
	while(rows>col)
	{
		for(auto x:prob)
			x.pb(0);
		col++;
	}
	vvi ori=prob;
	solve_assg(prob);
	vector<pii> v;
	vi c(ps,0);
	allot(v,prob,c);
	int Z=0;
	for(auto x:v)
	{
		cout<<x.first<<", "<<x.second<<endl;
		Z+=ori[x.first][x.second];
	}
	cout<<"Optimal cost : "<<Z<<endl;
	return 0;
}
/*
8 7 9 9
5 2 7 8
6 1 4 9
2 3 2 6
*/