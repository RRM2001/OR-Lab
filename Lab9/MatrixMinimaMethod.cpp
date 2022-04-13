#include <iostream>
#include <vector>
#include <stack>
#include <utility>
#include <cmath>
#include <cfloat>
#include <numeric>

#define vvd vector<vector<double> >
#define vd vector<double>
#define pb push_back
#define pii pair<int,int>

using namespace std;

bool make_loop(stack<pii>&s, vvd &table, pii entering, pii current, pii prev)
{
	//cout<<current.first<<", "<<current.second<<endl;
	if(current==entering) return true;
	vector<pii>v;
	while(!s.empty())
	{
		//if(s.top()==current)return false;
		v.pb(s.top());
		//cout<<s.top().first<<", "<<s.top().second<<'\t';
		s.pop();
	}
	//cout<<endl;
	//if(current==entering) return true;
	for(int i=v.size()-1;i>=0;i--)
		s.push(v[i]);
	if(current.first==-1)
	{
		bool ret;
		if(entering.first==0)
		{
			ret = make_loop(s,table,entering,make_pair(1,entering.second),entering);
			if(ret) return ret;
			if(entering.second==0)
				ret=ret||make_loop(s,table,entering, make_pair(0,entering.second+1),entering);
			else if(entering.second==table[0].size()-1)
				ret=ret||make_loop(s,table,entering, make_pair(0,entering.second-1),entering);
			else
			{
				ret=ret||make_loop(s,table,entering, make_pair(0,entering.second+1),entering);
				//if(ret) return ret;
				ret=ret||make_loop(s,table,entering, make_pair(0,entering.second-1),entering);
			}
		}
		else if(entering.first==table.size()-1)
		{
			ret = make_loop(s,table,entering, make_pair(table.size()-2,entering.second),entering);
			if(ret) return ret;
			if(entering.second==0)
				ret=ret||make_loop(s,table,entering,make_pair(table.size()-1,entering.second+1),entering);
			else if(entering.second==table[0].size()-1)
				ret=ret||make_loop(s,table,entering,make_pair(table.size()-1,entering.second-1),entering);
			else
			{
				ret=ret||make_loop(s,table,entering,make_pair(table.size()-1,entering.second+1),entering);
				//if(ret) return ret;
				ret=ret||make_loop(s,table,entering,make_pair(table.size()-1,entering.second-1),entering);
			}
		}
		else
		{
			ret = make_loop(s, table, entering, make_pair(entering.first+1,entering.second),entering);
			if(ret) return ret;
			ret=ret||make_loop(s, table, entering, make_pair(entering.first-1,entering.second),entering);
			if(ret) return ret;
			ret=ret||make_loop(s, table, entering, make_pair(entering.first,entering.second-1),entering);
			if(ret) return ret;
			ret=ret||make_loop(s, table, entering, make_pair(entering.first,entering.second+1),entering);
		}
		return ret;
	}
	if(table[current.first][current.second]==DBL_MIN)
	{
		if(current.first==table.size()-1 && current.first!=prev.first) return false;
		if(current.first==0 && current.first!=prev.first) return false;
		if(current.second==0 && current.second!=prev.second) return false;
		if(current.second==table[0].size()-1 && current.second!=prev.second) return false;
		pii new_curr;
		bool ret=false;
		if(current.first==prev.first)
		{
			if(current.second==prev.second+1)
				new_curr=make_pair(current.first,current.second+1);
			else
				new_curr=make_pair(current.first,current.second-1);
		}
		else
		{
			if(current.first==prev.first+1)
				new_curr=make_pair(current.first+1,current.second);
			else
				new_curr=make_pair(current.first-1,current.second);
		}
		s.push(current);
		ret=make_loop(s,table,entering,new_curr,current);
		if(!ret)s.pop();
		return ret;
	}
	else
	{
		s.push(current);
		bool ret=false;
		if(current.first==prev.first)
		{
			if(current.second==prev.second+1)
			{
				if(current.second==table[0].size()-1)
				{
					if(current.first==0)
						ret=ret||make_loop(s,table,entering,make_pair(1,current.second),current);
					else if(current.first==table.size()-1)
						ret=ret||make_loop(s,table,entering,make_pair(table.size()-2,current.second),current);
					else
					{
						ret=ret||make_loop(s,table,entering,make_pair(current.first+1,current.second),current);
						if(ret) return ret;
						ret=ret||make_loop(s,table,entering,make_pair(current.first-1,current.second),current);
					}
				}
				else
				{
					ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second+1),current);
					if(ret) return ret;
					if(current.first==0)
						ret=ret||make_loop(s,table,entering,make_pair(current.first+1,current.second),current);
					else if(current.first==table.size()-1)
						ret=ret||make_loop(s,table,entering,make_pair(current.first-1,current.second),current);
					else
					{
						ret=ret||make_loop(s,table,entering,make_pair(current.first-1,current.second),current);
						//if(ret) return ret;
						ret=ret||make_loop(s,table,entering,make_pair(current.first+1,current.second),current);
					}
				}
			}
			else
			{
				if(current.second==0)
				{
					if(current.first==0)
						ret=ret||make_loop(s,table,entering,make_pair(1,0),current);
					else if(current.first==table.size()-1)
						ret=ret||make_loop(s,table,entering,make_pair(table.size()-2,0),current);
					else
					{
						ret=ret||make_loop(s,table,entering,make_pair(current.first+1,0),current);
						//if(ret) return ret;
						ret=ret||make_loop(s,table,entering,make_pair(current.first-1,0),current);
					}
				}
				else
				{
					ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second-1),current);
					if(ret) return ret;
					if(current.first==0)
						ret=ret||make_loop(s,table,entering,make_pair(current.first+1,current.second),current);
					else if(current.first==table.size()-1)
						ret=ret||make_loop(s,table,entering,make_pair(current.first-1,current.second),current);
					else
					{
						ret=ret||make_loop(s,table,entering,make_pair(current.first-1,current.second),current);
						//if(ret) return ret;
						ret=ret||make_loop(s,table,entering,make_pair(current.first+1,current.second),current);
					}
				}
			}
		}
		else
		{
			if(current.first==prev.first+1)
			{
				if(current.first==table.size()-1)
				{
					if(current.second==0)
						ret=ret||make_loop(s,table,entering,make_pair(current.first,1),current);
					else if(current.second==table[0].size()-1)
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second-1),current);
					else
					{
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second+1),current);
						//if(ret) return ret;
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second-1),current);
					}
				}
				else
				{
					ret=ret||make_loop(s,table,entering,make_pair(current.first+1,current.second),current);
					if(ret) return ret;
					if(current.second==0)
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second+1),current);
					else if(current.second==table[0].size()-1)
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second-1),current);
					else
					{
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second+1),current);
						//if(ret) return ret;
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second-1),current);
					}
				}
			}
			else
			{
				if(current.first==0)
				{
					if(current.second==0)
						ret=ret||make_loop(s,table,entering,make_pair(0,1),current);
					else if(current.second==table[0].size()-1)
						ret=ret||make_loop(s,table,entering,make_pair(0,table[0].size()-2),current);
					else
					{
						ret=ret||make_loop(s,table,entering,make_pair(0,current.second-1),current);
						//if(ret) return ret;
						ret=ret||make_loop(s,table,entering,make_pair(0,current.second+1),current);
					}
				}
				else
				{
					ret=ret||make_loop(s,table,entering,make_pair(current.first-1,current.second),current);
					if(ret) return ret;
					if(current.second==0)
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second+1),current);
					else if(current.second==table[0].size()-1)
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second-1),current);
					else
					{
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second+1),current);
						//if(ret) return ret;
						ret=ret||make_loop(s,table,entering,make_pair(current.first,current.second-1),current);
					}
				}
			}
		}
		if(!ret)
		{
			s.pop();
		}
		return ret;
	}
}

int main()
{
	double theta=0;
	int dds, sss;
	cout<<"Enter number of demand points :";
	cin>>dds;
	vd demand(dds,0);//{3,3,2,2};
	cout<<"Enter number of supply sources:";
	cin>>sss;
	vd supply(sss,0);//{5,2,3};
	cout<<"Enter demand :"<<endl;
	for(int i=0;i<dds;i++)
		cin>>demand[i];
	cout<<"Enter supply :"<<endl;
	for(int i=0;i<sss;i++)
		cin>>supply[i];
	vvd costs(supply.size(),vd(demand.size()));
	cout<<"Enter costs :"<<endl;
	for(int i=0;i<sss;i++)
	{
		for(int j=0;j<dds;j++)
			cin>>costs[i][j];
	}
	double sup=accumulate(supply.begin(),supply.end(),0);
	double dd=accumulate(demand.begin(),demand.end(),0);
	if(sup>dd)
	{
		demand.pb(sup-dd);
		for(auto x:costs)
			x.pb(DBL_MIN);
	}
	else if(dd>sup)
	{
		supply.pb(dd-sup);
		vd v;
		for(auto x:costs[0])
			v.pb(DBL_MIN);
		costs.pb(v);
	}
	vvd table(supply.size(),vd(demand.size(),DBL_MIN));
	double min_cost=DBL_MAX;
	double alloc=-1;
	pii loc;
	while(*max_element(supply.begin(),supply.end())!=0 && *max_element(demand.begin(),demand.end())!=0)
	{
		for(int i=0;i<costs.size();i++)
		{
			for(int j=0;j<costs[i].size();j++)
			{
				if(supply[i]!=0 && demand[j]!=0)
				{
					
					if(min_cost==costs[i][j])
					{
						if(alloc<min(supply[i],demand[j]))
						{
							min_cost=costs[i][j];
							alloc=min(supply[i],demand[j]);
							loc=make_pair(i,j);
						}
					}
					else if(min_cost>costs[i][j])
					{
						min_cost=costs[i][j];
						alloc=min(supply[i],demand[j]);
						loc=make_pair(i,j);
					}
				}
			}
		}
		table[loc.first][loc.second]=alloc;
		demand[loc.second]-=alloc;
		supply[loc.first]-=alloc;
		min_cost=DBL_MAX;
		alloc=-1;
	}
	cout<<endl<<endl<<endl<<"Matrix Minima Method"<<endl;
	while(1){
	cout<<"BFS :"<<endl;
	for(auto x:table)
	{
		for(auto y:x)
		{
			if(y==DBL_MIN)cout<<"0\t";
			else cout<<y<<'\t';
		}
		cout<<endl;
	}
	cout<<endl<<endl<<endl;
	vd u(supply.size(),-100000000);
	vd v(demand.size(),-100000000);
	u[0]=0;
	while(*min_element(u.begin(),u.end())==-100000000||*min_element(v.begin(),v.end())==-100000000)
	{
		for(int i=0;i<table.size();i++)
		{
			for(int j=0;j<table[0].size();j++)
			{
				if(table[i][j]!=DBL_MIN)
				{
					if(u[i]!=-100000000 && v[j]==-100000000)
					{
						v[j]=costs[i][j]-u[i];
					}
					else if(u[i]==-100000000 && v[j]!=-100000000)
					{
						u[i]=costs[i][j]-v[j];
					}
				}
			}
		}
	}
	vvd enter(table.size(),vd(table[0].size(),DBL_MIN));
	double ma=-100000000;
	pii entering;
	for(int i=0;i<table.size();i++)
	{
		for(int j=0;j<table[0].size();j++)
		{
			if(table[i][j]==DBL_MIN)
			{
				enter[i][j]=u[i]+v[j]-costs[i][j];
				if(ma<enter[i][j])
				{
					ma=enter[i][j];
					entering=make_pair(i,j);
				}
			}
		}
	}
	cout<<"ui+vj-cij for various locations :"<<endl;
	for(auto x:enter)
	{
		for(auto y:x)
		{
			if(y==DBL_MIN)cout<<"NA\t";
			else cout<<y<<'\t';
		}
		cout<<endl;
	}
	cout<<endl<<endl;
	cout<<"Maximum of ui+vj-cij = "<<ma<<endl;
	if(ma>0) cout<<"Entering location "<<entering.first<<", "<<entering.second<<endl<<endl;
	if(ma<0)
	{
		//optimal
		cout<<"Optimal"<<endl;
		double Z=0;
		for(int i=0;i<table.size();i++)
		{
			for(int j=0;j<table[0].size();j++)
			{
				if(table[i][j]!=DBL_MIN)
					Z+=costs[i][j]*table[i][j];
			}
		}
		cout<<endl<<"Minimum cost = "<<Z<<endl;
		return 0;
	}
	stack<pii>s;
	s.push(entering);
	bool ret=make_loop(s,table,entering,make_pair(-1,-1),make_pair(-1,-1));
	//cout<<"ret "<<ret<<endl<<endl;
	vector<pii>vt;
	vector<pii>loop;
	while(!s.empty())
	{
		//cout<<s.top().first<<", "<<s.top().second<<endl;
		if(table[s.top().first][s.top().second]!=DBL_MIN)vt.pb(s.top());
		s.pop();
	}
	//cout<<endl<<endl;
	for(int i=0;i<vt.size();i++)
	{
		if(i==0)
		{
			if((entering.first==vt[i].first&&vt[i].first==vt[i+1].first)||(entering.second==vt[i].second&&vt[i].second==vt[i+1].second))
			continue;
		}
		else if(i==vt.size()-1)
		{
			if((entering.first==vt[i].first&&vt[i].first==vt[i-1].first)||(entering.second==vt[i].second&&vt[i].second==vt[i-1].second))
			continue;
		}
		else
		{
			if((vt[i-1].first==vt[i].first&&vt[i].first==vt[i+1].first)||(vt[i-1].second==vt[i].second&&vt[i].second==vt[i+1].second))
			continue;
		}
		loop.pb(vt[i]);
	}
	cout<<endl<<"loop"<<endl;
	for(auto x:loop)
		cout<<x.first<<", "<<x.second<<endl;
	cout<<endl;
	theta=min(table[loop[0].first][loop[0].second],table[loop[loop.size()-1].first][loop[loop.size()-1].second]);
	cout<<"Theta = "<<theta<<endl<<endl;
	table[entering.first][entering.second]=theta;
	if(table[loop[0].first][loop[0].second]==table[loop[loop.size()-1].first][loop[loop.size()-1].second])
	{
		for(int i=0;i<loop.size();i++)
		{
			table[loop[i].first][loop[i].second]+=theta*pow(-1,i+1);
		}
		if(costs[loop[0].first][loop[0].second]<costs[loop[loop.size()-1].first][loop[loop.size()-1].second])
		table[loop[loop.size()-1].first][loop[loop.size()-1].second]=DBL_MIN;
		else
		table[loop[0].first][loop[0].second]=DBL_MIN;
	}
	else
	{
		for(int i=0;i<loop.size();i++)
		{
			table[loop[i].first][loop[i].second]+=theta*pow(-1,i+1);
		}
		if(table[loop[loop.size()-1].first][loop[loop.size()-1].second]==0)
			table[loop[loop.size()-1].first][loop[loop.size()-1].second]=DBL_MIN;
		else
			table[loop[0].first][loop[0].second]=DBL_MIN;
	}
	}
	return 0;
}