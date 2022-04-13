//Rishabh Rathi
//18HS20027
//Revised Simplex

#include <iostream>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <iomanip>
#include <cmath>

using namespace std;

#define vd vector<double>
#define vi vector<int>
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

void print_b(vector<vd>&b)
{
	for(auto x:b)
	{
		for(auto y:x)
			cout<<y<<' ';
		cout<<endl;
	}
	cout<<endl;
	return;
}
void print(vd m)
{
	for(auto x:m)
		if(x==DBL_MAX)
		cout<<"NA ";
		else
		cout<<x<<' ';
	cout<<endl;
}
void printi(vi m)
{
	for(auto x:m)
		cout<<x<<' ';
	cout<<endl;
}

void construct_con(vi &art, vi &basis,vector<vd>&a,vi &con, int n)
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

vector<vd> inverse(vector<vd>b, vector<vd>b_1, vector<vd> &binv_1)
{
	int col=0;
	for(int i=0;i<b.size();i++)
	{
		for(int j=0;j<b.size();j++)
		{
			if(b[j][i]!=b_1[j][i])
			{
				col=i;
				break;
			}
		}
		if(col) break;
	}
	vd e(b.size(), 0);
	for(int i=0;i<e.size();i++)
	{
		for(int j=0;j<b.size();j++)
		{
			e[i]+=binv_1[i][j]*b[j][col];
		}
	}
	double ecol=e[col];
	e[col]=1/e[col];
	for(int i=0;i<e.size();i++)
	{
		if(i!=col)
		{
			e[i]/=(-ecol);
		}
	}
	vector<vd>Er(e.size(),vd(e.size(),0));
	for(int i=0;i<Er.size();i++)
	{
		Er[i][i]=1;
		Er[i][col]=e[i];
	}
	vector<vd>binv(b.size(),vd(b.size(),0));
	for(int i=0;i<binv.size();i++)
	{
		for(int j=0;j<binv.size();j++)
		{
			for(int k=0;k<binv_1.size();k++)
			{
				binv[i][j]+=Er[i][k]*binv_1[k][j];
			}
		}
	}
	return binv;
}

int solve_revised_simplex(vector<vd>&P, vd &rhs, vector<vd>b,vector<vector<vd> >&B, vector<vd>binv, vector<vector<vd> >&Binv, vd &z, vd &zm, vd cb, vector<vd> &CB, vd cbm, vector<vd> &CBM, vi basis, vector<vi> &BASIS)
{
	vd y(binv.size(),0);
	vd ym(binv.size(),0);
	for(int i=0;i<y.size();i++)
	{
		for(int j=0;j<y.size();j++)
		{
			y[i]+=cb[j]*binv[j][i];
			ym[i]+=cbm[j]*binv[j][i];
		}
	}
	vd nbz(P.size(),DBL_MAX);
	vd nbzm(P.size(),DBL_MAX);
	for(int i=0;i<nbz.size();i++)
	{
		auto it=find(basis.begin(),basis.end(),i+1);
		if(it==basis.end())
		{
			int id=it-basis.begin();
			nbz[i]=0;
			nbzm[i]=0;
			for(int j=0;j<y.size();j++)
			{
				nbz[i]+=y[j]*P[i][j];
				nbzm[i]+=ym[j]*P[i][j];
			}
			nbz[i]-=z[i];
			nbzm[i]-=zm[i];
		}
	}
	print(nbz);
	int id= min_element(nbzm.begin(),nbzm.end())-nbzm.begin();
	int entering;
	if(nbzm[id]>=0)
	{
		id=min_element(nbz.begin(),nbz.end())-nbz.begin();
		if(nbz[id]>=0)
			return 0;//solved
		else entering=id;
	}
	else
	{
		entering=id;
	}
	cout<<"Entering id = "<<entering<<endl;
	vd xb(basis.size(),0);
	for(int i=0;i<xb.size();i++)
	{
		for(int j=0;j<xb.size();j++)
		{
			xb[i]+=binv[i][j]*rhs[j];
		}
	}
	cout<<"Xb"<<endl;
	print(xb);
	vd alpha(basis.size(),0);
	for(int i=0;i<alpha.size();i++)
	{
		for(int j=0;j<alpha.size();j++)
		{
			alpha[i]+=binv[i][j]*P[id][j];
		}
	}
	cout<<endl<<"Alpha"<<endl;
	print(alpha);
	vd theta(basis.size(),DBL_MAX);
	int ch=1;
	for(int i=0;i<theta.size();i++)
	{
		if(xb[i]<=0) continue;
		if(alpha[i]<=0) continue;
		theta[i]=xb[i]/alpha[i];
		ch=0;
	}
	cout<<endl<<"Theta"<<endl;
	print(theta);
	if(ch)
	{
		return 2;//unbounded
	}
	int leaving=min_element(theta.begin(),theta.end())-theta.begin();
	cout<<endl<<"leaving id : "<<leaving<<endl;
	basis[leaving]=entering+1;
	for(int i=0;i<P[entering].size();i++)
	{
		b[i][leaving]=P[entering][i];
	}
	cb[leaving]=z[basis[leaving]-1];
	cbm[leaving]=zm[basis[leaving]-1];
	binv=inverse(b,B[B.size()-1],binv);
	B.pb(b);
	Binv.pb(binv);
	BASIS.pb(basis);
	CB.pb(cb);
	CBM.pb(cbm);
	return 1;//unsolved
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
	vi con(n);
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
	vector<vi> basis(1,vi(n,0));
	vi art;
	construct_con(art,basis[0],a,con,n);
	while(z.size()<a[0].size())
	{
		z.pb(0);
		zm.pb(0);
	}
	for(auto x:art)
		zm[x-1]=-1;
	vector<vd> cb;
	vector<vd> cbm;
	vd v1(n);
	vd v2(n);
	for(int i=0;i<n;i++)
	{
		vd v1(n);
		vd v2(n);
		v1[i]=zm[basis[0][i]-1];
		v2[i]=z[basis[0][i]-1];
	}
	cb.pb(v2);
	cbm.pb(v1);
	vector<vector<vd> >B;
	vector<vector<vd> >Binv;
	vector<vd>w(n,vd(n,0));
	for(int i=0;i<n;i++)
	{
		w[i][i]=1;
	}
	B.pb(w);
	Binv.pb(w);
	vector<vd>P;
	for(int i=0;i<a[0].size();i++)
	{
		vd v(a.size(),0);
		for(int j=0;j<a.size();j++)
		{
			v[j]=a[j][i];
		}
		P.pb(v);
	}
	int sol=1;
	for(auto x:P)
	{
		for(auto y:x)
			cout<<y<<' ';
		cout<<endl;
	}
	cout<<endl;
	while(sol==1)
	{
		/*print_b(B[B.size()-1]);
		cout<<endl;
		print_b(Binv[Binv.size()-1]);
		cout<<endl;
		print(cb[cb.size()-1]);
		cout<<endl;
		printi(basis[basis.size()-1]);
		cout<<endl;*/
		sol=solve_revised_simplex(P,rhs,B[B.size()-1],B,Binv[Binv.size()-1],Binv,z,zm,cb[cb.size()-1],cb,cbm[cbm.size()-1],cbm, basis[basis.size()-1],basis);
	}
	if(sol==0)
	{
		//calculate xb and z here
		vd xb(Binv[0].size(),0);
		for(int i=0;i<xb.size();i++)
		{
			for(int j=0;j<Binv[Binv.size()-1].size();j++)
				xb[i]+=Binv[Binv.size()-1][i][j]*rhs[j];
		}
		double op=0;
		for(int i=0;i<xb.size();i++)
			op+=(cb[cb.size()-1][i]*xb[i]);
		cout<<"Solved : "<<op<<endl;
		int c=0;
		for(int i=0;i<xb.size();i++)
		{
			if(basis[basis.size()-1][i]>m)continue;
			cout<<fixed<<setprecision(0)<<"x"<<basis[basis.size()-1][i];
			cout<<fixed<<setprecision(2)<<" = "<<xb[i]<<endl;
			c++;
		}
		if(c<xb.size())
			cout<<"The rest of the original variables are 0.";
	}
	else cout<<"Unbounded";
	return 0;
}