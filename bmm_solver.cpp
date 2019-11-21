#include <utility>
#include <vector>
using namespace std;


double first_moment_beta(double a, double b){
	return a / (a + b);
}


double second_moment_beta(double a, double b){
	return a * (a + 1) / ((a+b+1)*(a+b));
}


pair<double, double> solve_beta(double first_moment, double second_moment){
	double f = first_moment;
	double s = second_moment;

	double tau = (s-f)/(f*f - s);

	pair<double, double> beta;

	beta.first = f * tau;
	beta.second = tau - a;

	return beta;
}


void positive_update(double &a, double &b, double p){
	double first_moment = first_moment_beta(a, b);
	double second_moment = second_moment_beta(a, b);

	double posterior_first_moment = (first_moment - p * first_moment_beta(a, b+1))/(1 - p);
	double posterior_second_moment = (second_moment - p * second_moment_beta(a, b+1))/(1 - p);
	
	pair<double, double> new_beta = solve(posterior_first_moment, posterior_second_moment);

	a = new_beta.first;
	b = new_beta.second;
}


void negative_update(double &a, double &b, double p){
	double first_moment = first_moment_beta(a, b);
	double second_moment = second_moment_beta(a, b);

	double posterior_first_moment = (first_moment - p * first_moment_beta(a+1, b))/(1 - p);
	double posterior_second_moment = (second_moment - p * second_moment_beta(a+1, b))/(1 - p);
	
	pair<double, double> new_beta = solve(posterior_first_moment, posterior_second_moment);

	a = new_beta.first;
	b = new_beta.second;
}

double abs(double m){
	if(m > 0.0){
		return m;
	} else{
		return -m;
	}
}

vector<double> ans;

vector<double> bns;

void update_literal(vector<int> &clause){
	float p = 1.0;

	for(int & literal : clause){
		if(literal > 0){
			p = p * (1 - first_moment_beta(ans[literal-1], bns[literal-1]));
		} else{
		        p = p * first_moment_beta(ans[abs(literal)-1], bns[abs(literal) - 1]);
		}
	}

	for(int & literal : clause){
		if(literal > 0){
			positive_update(ans[literal-1], bns[literal -1], p);
		} else{
		        negative_update(ans[abs(literal)-1], bns[abs(literal) -1], p);
		}
	}
}

