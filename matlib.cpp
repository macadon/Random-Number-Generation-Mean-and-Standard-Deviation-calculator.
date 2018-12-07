#include "matlib.h"
#include "testing.h"
#include "geometry.h"
#include <cmath>
#include <vector>
#include <stdlib.h>
using namespace std;

static inline double Homerzero(double x, double a0) {
        return a0;
}

static inline double Homerone(double x, double a0, double a1){
        return a0+x*Homerzero(x,a1);
}

static inline double Homertwo(double x, double a0, double a1, double a2) {
        return a0+x*Homerone( x,  a1,  a2);
}

static inline double Homerthree(double x, double a0, double a1, double a2, double a3) {
        return a0+x*Homertwo( x,  a1,  a2,  a3);
}

static inline double Homerfour(double x, double a0, double a1, double a2, double a3, double a4) {
        return a0+x*Homerthree(x, a1,  a2,  a3,  a4);
}

static inline double Homerfive(double x, double a0, double a1, double a2, double a3, double a4, double a5) {
        return a0+x*Homerfour(x, a1, a2, a3, a4, a5);
        }
static inline double Homersix(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6){
        return a0+x*Homerfive(x,a1,a2,a3,a4,a5,a6);
}

static inline double Homerseven(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7) {
       return  a0+x*Homersix(x,a1,a2,a3,a4,a5,a6,a7);
}

static inline double Homereight(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8) {
        return a0+x*Homerseven(x,a1,a2,a3,a4,a5,a6,a7,a8);
}





double normpdf(double x) {

    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}
double normcdf(double x){
double k = 1.0/(1.0 + 0.2316419*x);
double k_sum = k*(0.31938153+k*(-0.356563782+k*(1.781477937+k*(-1.821255978+k*1.330274429))));
if(x>=0) {
return 1-(1/sqrt(2*PI))*exp(-x*x*0.5)*k_sum;
}
else {
return 1-normcdf(-x);
}
}
const double a0 = 2.50662823884;
const double a1 = -18.61500062529;
const double a2 = 41.39119773534;
const double a3 = -25.44106049637;
const double b1 = -8.47351093090;
const double b2 = 23.08336743743;
const double b3 = -21.06224101826;
const double b4 = 3.13082909833;
const double c0 = 0.3374754822726147;
const double c1 = 0.9761690190917186;
const double c2 = 0.1607979714918209;
const double c3 = 0.0276438810333863;
const double c4 = 0.0038405729373609;
const double c5 = 0.0003951896511919;
const double c6 = 0.0000321767881768;
const double c7 = 0.0000002888167364;
const double c8 = 0.0000003960315187;

double norminv( double x ) {
    // We use Moro's algorithm
    double y = x - 0.5;
    if (y<0.42 || y>-0.42) {
        double r = y*y;
        return y*Homerthree(r,a0,a1,a2,a3)/Homerfour(r,1.0,b1,b2,b3,b4);
    } else {
        double r;
        if (y<0.0) {
            r = x;
        } else {
            r = 1.0 - x;
        }
        double s = log( -log( r ));
        double t = Homereight(s,c0,c1,c2,c3,c4,c5,c6,c7,c8);
        if (x>0.5) {
            return t;
        } else {
            return -t;
        }
    }
}

void solvequadratic(double a, double b, double c, double &z1, double &z2) {
if(b*b<4*a*c) {
cout<<"No real roots!"<<endl;
}
else if (b*b==4*a*c) {
z1=-b/(2*a);
cout<<"Exactly one real root!"<<z1<<endl;
}
else {
cout<<"Two roots!"<<endl;
z1=(-b+sqrt(b*b-4*a*c))/2*a;
z2=(-b-sqrt(b*b-4*a*c))/2*a;
cout<<z1<<endl;
cout<<z2<<endl;
}





}
double mean(vector <double> v) {
int n=v.size();
double run;
for(int i=0;i<=n-1;i++) {
run+=v[i];
}
return run/n;

}

double sdev(vector <double> v) {
int n=v.size();
double x=mean(v);
double total;
for (int i=0;i<=n-1;i++) {
total+=(x-v[i])*(x-v[i]);


}
return sqrt(total/(n-1));



}
static void testnormcdf() {
  // test bounds
  ASSERT(normcdf(0.3)>0);
  ASSERT(normcdf(0.3)<1);
  // test extreme values
  ASSERT_APPROX_EQUAL(normcdf(-1e10), 0, 0.001);
  ASSERT_APPROX_EQUAL(normcdf(1e10), 1.0, 0.001);
  // test increasing
  ASSERT(normcdf(0.3)<normcdf(0.5));
  // test symmetry
  ASSERT_APPROX_EQUAL(normcdf(0.3),
    1 - normcdf(-0.3), 0.0001);
  ASSERT_APPROX_EQUAL(normcdf(0.0), 0.5, 0.0001);
  // test inverse
  ASSERT_APPROX_EQUAL(normcdf(norminv(0.3)),
    0.3, 0.0001);
  // test well known value
  ASSERT_APPROX_EQUAL(normcdf(1.96), 0.975, 0.001);
}

static void testnorminv() {
    ASSERT_APPROX_EQUAL(norminv(1.96), 0.675, 0.01 );
}

void testMatlib() {
    TEST( testnorminv );
    TEST( testnormcdf );
}

void sort(vector <double> &v) {
int n=v.size();
for(int j=0;j<=n-1;j++) {

int iMin=j;

for(int i=j+1;i<n;i++) {

if(v[i]<v[iMin]) {
iMin=i;
}
}

if(iMin!=j) {
std::swap(v[j],v[iMin]);
}

}
}


void swap(double a, double b) {
double temp=a;
a=b;
b=temp;
}


void randuniform(int n, vector <double> &v) {
for(int i=0;i<n;i++) {
v[i]=((double)rand())/((double)RAND_MAX+1);

}

}

double unif() {
return ((double)rand())/((double)RAND_MAX+1);

}

double BM(int N, vector <double> & t) {
for(int i=0;i<N;i++) {
double u1=unif();
double u2=unif();
double n1=sqrt(-2*log(u1))*cos(2*PI*u2);
t[i]=n1;

}
}
