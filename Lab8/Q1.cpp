//Rishabh Rathi
//18HS20027
//Cutting plane method

#include<iostream>
#include<vector>
#include<algorithm>
#include<cmath>

using namespace std;
#define vd vector<double>
#define vvd vector<vector<double> >
#define vi vector<int>
#define pb push_back
#define epsilon_threshold 0.00001
#define double_MAXX 1000000
#define zero_threshold 0.0001
#define counter_overflow_threshold 15
#define M_value 1000000

bool is_infeasible=0;

void print_v(vi arr)
{
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

void print_vd(vd arr)
{
	for(int i=0;i<arr.size();i++)
		cout<<arr[i]<<" ";
	cout<<endl;
}

void print_mat(vvd mat)
{
	for(int i=0;i<mat.size();i++){
		for(int j=0;j<mat[i].size();j++)
			cout<<mat[i][j]<<" ";
		cout<<endl;
	}
}

//Find dot product of 2 n dimensional vectors
double dot_prod(vd arr1, vd arr2)
{
	int n = arr1.size();
	double ret=0;
	for(int i=0;i<n;i++)
		ret+=(arr1[i]*arr2[i]);
	return ret;
}


//Modify the objective function to take the slack variables into account
void create_max_obj(vd &objective_fn, bool is_max, int n, vi slack_list, vi surplus_list, vi artificial_list)
{
	for(int i=objective_fn.size();i<n;i++) objective_fn.pb(0);
	cout<<"objective_fn_size: "<<objective_fn.size()<<endl;
	cout<<"n: "<<n<<endl;
	if(!is_max){
		for(int i=0;i<objective_fn.size();i++) objective_fn[i] = -1*objective_fn[i];
	}
	
	//Use artificial variables to modify the objective function to include the bigM terms
	for(int i=0;i<artificial_list.size();i++)
		objective_fn[artificial_list[i]-1] = -1*M_value;
}

//Create the full coefficient matrix which includes the slack variables as well
void create_mat(vi eqn_types, vvd &mat, vi &slack_list, vi &surplus_list, vi &artificial_list)
{
	int m = mat.size();
	int n = mat[0].size();

	vd last_column(m,0);
	for(int i=0;i<m;i++){
		last_column[i] = mat[i][n-1];
		mat[i].pop_back();
	}
	n--;//n is the number of variables till now

	//Add the slack, surplus, artificial variables to the matrix
	for(int i=0;i<m;i++){
		if(eqn_types[i]==1){
			slack_list.pb(n+1); n++;
			for(int j=0;j<m;j++){
				if(j==i) mat[j].pb(1);
				else mat[j].pb(0);
			}
		}

		else if(eqn_types[i]==2){
			surplus_list.pb(n+1);
			artificial_list.pb(n+2);
			n=n+2;//add surplus and artificial variable
			for(int j=0;j<m;j++){
				if(j==i){mat[j].pb(-1); mat[j].pb(1);}//subtract surplus and add artificial variable 
				else {mat[j].pb(0); mat[j].pb(0);}
			}
		}

		else if(eqn_types[i]==3){
			artificial_list.pb(n+1); n++;//add artificial variable
			for(int j=0;j<m;j++){
				if(j==i) mat[j].pb(1);
				else mat[j].pb(0);
			}
		}
	}

	//Add the RHS of the contraints back again into the matrix
	for(int i=0;i<m;i++){
		mat[i].pb(last_column[i]);
	}
}

bool can_simplex_be_used(vi eqn_types)
{
	for(int i=0;i<eqn_types.size();i++){
		if(eqn_types[i]!=1){cout<<endl<<"Normal Simplex method cant be used to solve this system, using BIG_M method"<<endl; return 0;}
	}
	cout<<"Using simplex method to solve optimisation question"<<endl;
	return 1;
}

vvd construct_bigM_tableau(vvd mat, vd objective_fn, vi slack_list, vi surplus_list, vi artificial_list)
{
	int num_rows = mat.size();
	int num_cols = mat[0].size();
	vector<vector<double> > tableau(num_rows , vector<double> (num_cols+2, 0));//2 extra columns for C0, Basis
	
	int curr_row=0;
	int slack_index=0;
	int artificial_index=0;
	while( (slack_index < slack_list.size()) && (artificial_index < artificial_list.size()) ){
		if(slack_list[slack_index] < artificial_list[artificial_index]){
			tableau[curr_row][1] = slack_list[slack_index++] + 1;//column index in tableau
			tableau[curr_row++][0] = 0;
		}
		else{
			tableau[curr_row][1] = artificial_list[artificial_index++] + 1;
			tableau[curr_row++][0] = -1*M_value;
		}
	}

	while(slack_index < slack_list.size()){
		tableau[curr_row][1] = slack_list[slack_index++] + 1;//column index in tableau
		tableau[curr_row++][0] = 0;
	}

	while(artificial_index < artificial_list.size()){
		tableau[curr_row][1] = artificial_list[artificial_index++] + 1;
		tableau[curr_row++][0] = -1*M_value;
	}


	//construct rest of the tableau
	for(int i=0;i<num_rows;i++){
		for(int j=0;j<num_cols;j++) tableau[i][j+2] = mat[i][j];
	}
	return tableau;
}


vd calculate_deviations(vvd tableau, vd objective_fn)
{
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	vd ret(num_cols-3,0);
	
	vd col(num_rows,0);
	for(int i=0;i<num_rows;i++) col[i] = tableau[i][0];

	for(int col_index = 2; col_index < num_cols-1; col_index++){
		vd temp(num_rows,0);		
		for(int i=0;i<num_rows;i++) temp[i] = tableau[i][col_index];

		ret[col_index-2] = objective_fn[col_index-2] - dot_prod(temp, col);
		temp.clear();
	}
	return ret;
}

bool is_termination_reached_dual_simplex(vvd tableau)
{
	int m = tableau.size();
	int n = tableau[0].size();
	for(int i=0;i<m;i++){
		if(tableau[i][n-1] < (-1*zero_threshold)) return 0;
	}
	return 1;
}

bool is_termination_reached(vd deviations)
{
	for(int i=0;i<deviations.size();i++){
		if(deviations[i]>zero_threshold) return 0;
	}
	return 1;
}

int find_max_elem_index(vd arr)
{
	double max = 0;
	int max_index=-1;
	for(int i=0;i<arr.size();i++){
		if(arr[i]>max){
			max = arr[i];
			max_index = i;
		}
	}
	return max_index;
}

// R1 -> R1 - dR2
void subtract_rows(vvd &tableau, double d, int r1, int r2)
{
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for(int i=2;i<num_cols;i++) tableau[r1][i] = tableau[r1][i] - d*tableau[r2][i];
}

// R1 -> R1/d
void divide_row(vvd &tableau, double d, int r)
{
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	for(int i=2;i<num_cols;i++) tableau[r][i] = tableau[r][i]/d;
}

bool is_in(int d, vi arr)
{
	for(int i=0;i<arr.size();i++){if(arr[i] == d) return 1;}
	return 0;
}

bool are_all_vars_int(vd variable_vals)
{
	for(int i=0;i<variable_vals.size();i++){
		int int_part = (int)(variable_vals[i]);
		if(variable_vals[i]-int_part > zero_threshold) return 0;
	}
	return 1;
}

int find_eqn_index_of_var_with_max_fractional(vd variable_vals, vvd tableau)
{
	int var_index=-1;
	double max_frac = -0.001;
	for(int i=0;i<variable_vals.size();i++){
		double curr_frac = variable_vals[i]-((int)(variable_vals[i]));
		if(curr_frac > (max_frac + epsilon_threshold)){
			var_index = i+1;
			max_frac = curr_frac;
		}
	}

	int eqn_index=-1;
	for(int i=0;i<tableau.size();i++){
		if(tableau[i][1]-1 == var_index){
			eqn_index = i; 
			break;
		}
	}
	return eqn_index;
}

bool is_int(double d)
{
	int gif_d = (int)d;
	double left = d - (double)gif_d;
	if( abs(left) < zero_threshold) return 1;
	if( 1 - abs(left) < zero_threshold) return 1;
	return 0;
}

vvd make_new_tableau(int eqn_index, vvd tableau)
{
	int rows = tableau.size();
	int cols = tableau[0].size();

	vvd new_tableau(rows+1, vd(cols+1,0));
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols-1;j++){
			new_tableau[i][j] = tableau[i][j];
		}
		new_tableau[i][cols-1] = 0;
		new_tableau[i][cols] = tableau[i][cols-1];
	}

	//Construct last row of the new tableau
	new_tableau[rows][0]=0;
	new_tableau[rows][1]=cols-1;
	new_tableau[rows][cols-1]=1;

	for(int i=2;i<cols-1;i++){
		double curr = tableau[eqn_index][i];
		if(is_int(curr)){
			new_tableau[rows][i] = 0;
		}
		else{
			if(curr<0){
				curr = -1*curr;
				double reqd = curr - (int)curr;
				reqd = 1-reqd;
				reqd = -1*reqd;
				new_tableau[rows][i] = reqd;
			}
			else{
				new_tableau[rows][i] = -1*(tableau[eqn_index][i] - ((int)tableau[eqn_index][i]));
			}
		}
	}
	double curr = tableau[eqn_index][cols-1];
	if(is_int(curr)){
		new_tableau[rows][cols]=0;
	}
	else if(curr>0){
		new_tableau[rows][cols] = -1*(tableau[eqn_index][cols-1] - ((int)tableau[eqn_index][cols-1]));
	}
	else{
		curr*=-1;
		double reqd = curr - (int)curr;
		reqd = 1-reqd;
		reqd*=-1;
		new_tableau[rows][cols] = reqd;
	}
	return new_tableau;
}

void solve_with_dual_simplex(vvd &tableau, vd &objective_fn, bool is_maximisation, bool &is_valid_integral_found, vd& variable_vals)
{
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();
	int count=0;
	while(1){
		if(count > counter_overflow_threshold){cout<<"Excess iterations, returning!!...."<<endl; return;}
		cout<<endl<<endl<<"Iteration: "<<count<<endl;
		cout<<"The Tableau: "<<endl; print_mat(tableau);

		bool is_end = is_termination_reached_dual_simplex(tableau);
		if(is_end){cout<<"Reached the termination state since all RHS values are positive"<<endl;break;}

		vd deviations = calculate_deviations(tableau, objective_fn);
		cout<<"Deviations: "; print_vd(deviations);

		//Identifying the leaving variable
		int leaving_var_row_index=-1;
		double min_till_now = zero_threshold;
		for(int i=0;i<num_rows;i++){
			if(tableau[i][num_cols-1] < min_till_now){
				min_till_now = tableau[i][num_cols-1];
				leaving_var_row_index = i;
			}
		}
		if(leaving_var_row_index == -1){cout<<"UNBOUNDED!! {Since we could not find a valid leaving variable}"<<endl; return;}
		cout<<"Leaving variable: x"<<tableau[leaving_var_row_index][1]-1<<endl;

		//Identifying the entering variable
		double min_positive_ratio_till_now = double_MAXX;
		int entering_var_col_index=-1;
		for(int i=2;i<num_cols-1;i++){
			if(tableau[leaving_var_row_index][i] < (-1*zero_threshold)){
				double temp_ratio = deviations[i-2]/tableau[leaving_var_row_index][i];
				if(temp_ratio>0 && temp_ratio<min_positive_ratio_till_now){
					min_positive_ratio_till_now = temp_ratio;
					entering_var_col_index=i;
				}
			}
		}
		if(entering_var_col_index==-1){
			cout<<"INFEASIBLE!! {Since we could not find valid entering variable}"<<endl; 
			is_infeasible=1;
			return;
		}
		cout<<"Entering variable: x"<<entering_var_col_index-1<<endl;


		//Perform row operations on the tableau to get the tableau for the next iteration ready
		for(int row_index = 0; row_index < num_rows; row_index++){
			if(row_index == leaving_var_row_index) divide_row(tableau, tableau[row_index][entering_var_col_index], row_index);
			else{
				double ratio_temp = tableau[row_index][entering_var_col_index]/tableau[leaving_var_row_index][entering_var_col_index];
				subtract_rows(tableau, ratio_temp, row_index, leaving_var_row_index);
			}
		}

		//Modify 1st and 2nd columns of the tableau
		tableau[leaving_var_row_index][0] = objective_fn[entering_var_col_index-2];
		tableau[leaving_var_row_index][1] = entering_var_col_index;

		deviations.clear();
		count++;
	}

	int num_vars = num_cols - 3;
	for(int i=0;i<num_rows;i++)	variable_vals[ tableau[i][1]-2 ] = tableau[i][num_cols-1];

	cout<<endl<<"The Optimisation function is optimised at the point : "; 
	for(int i=0;i<num_rows;i++){
		cout<<"x"<<tableau[i][1]-1<<" = "<<tableau[i][num_cols-1]<<", ";
	}
	cout<<"Rest all xis=0"<<endl;

	double optimal_value = dot_prod(objective_fn, variable_vals);	

	if(!is_maximisation) optimal_value*=-1;
	cout<<"The value of the objective function at this point is : "<<optimal_value<<endl;

	if(are_all_vars_int(variable_vals)){
		cout<<"The solution found is a valid integral solution"<<endl;
		is_valid_integral_found = 1;
	}
	return;
}

//First use BigM and then repeated adding new constraints and using Dual Simplex method to find integral optimal solution
void solve_problem(vvd tableau, vd objective_fn, bool is_max, vi artificial_list)
{
	int num_rows = tableau.size();
	int num_cols = tableau[0].size();

	cout<<endl; int count=0;
	while(1){
		if(count > counter_overflow_threshold){cout<<"Excess iterations, returning!!...."<<endl; return;}
		cout<<endl<<endl<<"Iteration: "<<count<<endl;
		cout<<"The Tableau: "<<endl; print_mat(tableau);

		vd deviations = calculate_deviations(tableau, objective_fn);
		cout<<"Deviations: "; print_vd(deviations);
		bool is_end = is_termination_reached(deviations);
		if(is_end){cout<<"Reached the termination state"<<endl;break;}

		//Identifying the entering variable
		int entering_var_col_index = find_max_elem_index(deviations) + 2;
		if(entering_var_col_index==-1){cout<<"LAPHDA!!"<<endl; break;}

		//Identifying the leaving variable
		int leaving_var_row_index=-1;
		double min_ratio = double_MAXX;
		for(int i=0;i<num_rows;i++){
			if(tableau[i][entering_var_col_index] <= 0) continue;
			double temp = tableau[i][num_cols-1]/tableau[i][entering_var_col_index];
			if(temp < min_ratio){
				min_ratio = temp;
				leaving_var_row_index = i;
			}
		}
		if(leaving_var_row_index==-1){cout<<"Unbounded solution!!!....EXITING..."<<endl; return;}
		cout<<"leaving_var_row_index: "<<leaving_var_row_index<<", entering_var_col_index: "<<entering_var_col_index<<endl;


		//Perform row operations on the tableau to get the tableau for the next iteration ready
		for(int row_index = 0; row_index < num_rows; row_index++){
			if(row_index == leaving_var_row_index) divide_row(tableau, tableau[row_index][entering_var_col_index], row_index);
			else{
				double ratio_temp = tableau[row_index][entering_var_col_index]/tableau[leaving_var_row_index][entering_var_col_index];
				subtract_rows(tableau, ratio_temp, row_index, leaving_var_row_index);
			}
		}

		//Modify 1st and 2nd columns of the tableau
		tableau[leaving_var_row_index][0] = objective_fn[entering_var_col_index-2];
		tableau[leaving_var_row_index][1] = entering_var_col_index;

		deviations.clear();
		count++;
	}

	int num_vars = num_cols - 3;
	vd variable_vals(num_vars,0);
	for(int i=0;i<num_rows;i++)	variable_vals[ tableau[i][1]-2 ] = tableau[i][num_cols-1];

	//Identify infeasible soution
	for(int i=0;i<num_rows;i++){
		if(is_in(tableau[i][1]-1, artificial_list)){
			cout<<endl<<"The iterations have been completed and there are artificial variables in the base with values strictly greater than 0, so the problem has no solution (INFEASIBLE!!!!!)."<<endl;
			return;
		}
	}

	cout<<endl<<"The Optimisation function is optimised at the point : "; 
	for(int i=0;i<num_rows;i++){
		cout<<"x"<<tableau[i][1]-1<<" = "<<tableau[i][num_cols-1]<<", ";
	}
	cout<<"Rest all xis=0"<<endl;

	double optimal_value = dot_prod(objective_fn, variable_vals);	
	if(!is_max) optimal_value *= -1;
	cout<<"The value of the objective function at this point is : "<<optimal_value<<endl;

	if(are_all_vars_int(variable_vals)){
		cout<<"The solution found is a valid integral solution"<<endl;
		return;
	}

	//Initial solution using BIG M Method done
	//Now we will find integral solution by
	//iteratively using dual simplex method 
	// and adding new constraints

	cout<<endl<<endl<<"BIG_M done, moving on to finding an integral solution"<<endl<<endl;
	int eqn_index;
	int iteration_index =0;
	while(iteration_index<5){
		cout<<endl<<endl<<"Iteration index: "<<iteration_index<<endl;
		eqn_index = find_eqn_index_of_var_with_max_fractional(variable_vals, tableau);
		tableau = make_new_tableau(eqn_index, tableau);

		//Now apply dual simplex method to the new tableau constructed after including the new constraint
		bool is_valid_integral_found = 0;
		solve_with_dual_simplex(tableau, objective_fn, is_max, is_valid_integral_found, variable_vals);
		if(is_valid_integral_found) break;
		iteration_index++;
	}
}


int main()
{
	int n, m;
	cout<<"Number of constraints : ";
	cin>>m;
	cout<<"Number of variables : ";
	cin>>n;
	vector<vector<double> > mat(m , vector<double> (n+1, 0)); 

	for(int i=0;i<m;i++){
		cout<<"Enter Constraint "<<i<<" coefficients: "<<endl;
		for(int j=0;j<n+1;j++){
			cin>>mat[i][j];
		}
	}
	cout<<endl<<"Enter whether the equations are <= type(1), >= type(2), or = type(3): "<<endl;
	vi eqn_types(m,0);
	for(int i=0;i<m;i++){
		cout<<"Equation "<<i<<": "; cin>>eqn_types[i]; 
	}
	bool is_problem_simplex = can_simplex_be_used(eqn_types);
	vi slack_variable_list; vi surplus_variable_list; vi artificial_variable_list;
	create_mat(eqn_types, mat, slack_variable_list, surplus_variable_list, artificial_variable_list);	

	print_mat(mat);
	cout<<"slack variables: "; print_v(slack_variable_list);
	cout<<"surplus variables: "; print_v(surplus_variable_list);
	cout<<"artifiical variables: "; print_v(artificial_variable_list);

	cout<<endl<<"Enter the objective function: "<<endl;
	vd objective_fn(n,0);
	for(int i=0;i<n;i++){
		cin>>objective_fn[i];
	}
	bool is_maximisation;
	cout<<endl<<"Does the objective_fn have to be maximised(1) or minimised(0)? : ";
	cin>>is_maximisation;
	n = mat[0].size();
	create_max_obj(objective_fn, is_maximisation, n-1, slack_variable_list, surplus_variable_list, artificial_variable_list);

	cout<<"The coeficient matrix: "<<endl;
	print_mat(mat);
	cout<<"The objective function: "<<endl;
	print_vd(objective_fn);
	cout<<endl;


	vvd tableau = construct_bigM_tableau(mat, objective_fn, slack_variable_list, surplus_variable_list, artificial_variable_list);
	solve_problem(tableau, objective_fn, is_maximisation, artificial_variable_list);

}