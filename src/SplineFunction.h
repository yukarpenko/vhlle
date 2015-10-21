#ifndef SPLINEFUNCTION_H
#define SPLINEFUNCTION_H
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

class SplineFunction
{
public:
    vector< pair<double, double> > vals;
    SplineFunction():vals(0) { vals.resize(0); }
    SplineFunction(vector<double> x, vector<double> y):vals(0)
    {
        for(unsigned int i=0;i<x.size();++i)
        {
            vals.push_back(make_pair(x[i], y[i]));
        }
        sort(vals.begin(), vals.end());
    }
    void add_val(double x, double val)
    {
        vals.push_back(make_pair(x, val));
        sort(vals.begin(), vals.end());
    }
    double f(double arg) const
    {
        unsigned int ind = 0;
        pair<double, double> op = make_pair(arg, 0.);
        vector< pair<double, double> >::const_iterator it = lower_bound(vals.begin(), vals.end(), op);
        ind = distance(vals.begin(), it);
        /*while (ind<vals.size() && arg>=vals[ind].first)
        {
            if (fabs(vals[ind].first-arg)<1e-10) return vals[ind].second;
            ind++;
        }*/
        if (ind==0) return vals[0].second +
            (arg - vals[0].first) *
            (vals[1].second - vals[0].second) / (vals[1].first - vals[0].first);
        if (ind==vals.size()) return vals[ind-2].second +
            (arg - vals[ind-2].first) *
            (vals[ind-1].second - vals[ind-2].second) / (vals[ind-1].first - vals[ind-2].first);
        return vals[ind-1].second +
            (arg - vals[ind-1].first) *
            (vals[ind].second - vals[ind-1].second) / (vals[ind].first - vals[ind-1].first);
    }
	double df(double arg) const
	{
		unsigned int ind = 0;
        pair<double, double> op = make_pair(arg, 0.);
        vector< pair<double, double> >::const_iterator it = lower_bound(vals.begin(), vals.end(), op);
        ind = distance(vals.begin(), it);
        if (ind==0) return (vals[1].second - vals[0].second) / (vals[1].first - vals[0].first);
        if (ind==vals.size()) return (vals[ind-1].second - vals[ind-2].second) / (vals[ind-1].first - vals[ind-2].first);
        return (vals[ind].second - vals[ind-1].second) / (vals[ind].first - vals[ind-1].first);
	}
    double fsquare(double arg)
    {
        double ret = f(arg);
        return ret*ret;
    }
    void clear()
    {
        vals.resize(2);
        vals[0].first = 0.;
        vals[0].second = 0.;
        vals[1].first = 1.;
        vals[1].second = 0.;
    }
    void clearall()
    {
        vals.resize(0);
    }
    void fill(vector<double> x, vector<double> y)
    {
        vals.resize(0);
        for(unsigned int i=0;i<x.size();++i)
        {
            vals.push_back(make_pair(x[i], y[i]));
        }
        sort(vals.begin(), vals.end());
    }
    void setConstant(double val)
    {
        vals.resize(0);
        vals.push_back(make_pair(0., val));
        vals.push_back(make_pair(1., val));
    }
    //void loadFromFile(const char *file);
};

/*SplineFunction::SplineFunction()
{
    vals.resize(0);
}*/

#endif // SPLINEFUNCTION_H
