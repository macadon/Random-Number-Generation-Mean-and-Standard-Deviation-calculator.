#include <iostream>
#include "matlib.cpp"
#include "geometry.cpp"
#include "testing.cpp"
using namespace std;

int main() {
        double x,y,r,a,b,c,z1,z2;
        vector<double> v(5);
        vector<double> w(100);
        vector<double> t(100);
        cout<<"normcdf(1.96) = "<<normcdf(1.96)<<"\n";
        cout<<"norminv(0.975) = "<<norminv(0.975)<<"\n";

        cout<<"Enter a radius!"<<endl;
        cin>>r;
        cout<<"The area of the circle with that radius is: "<<area(r)<<endl;
        cout<<"With circumference equal to "<<circumference(r)<<endl;
        cout<<"Normal at 1.96: "<<normcdf(1.96);
        testMatlib();
        cout<<"Enter the coefficients of the quadratic:"<<endl;
        cin>>a;
        cin>>b;
        cin>>c;
        solvequadratic(a,b,c,z1,z2);

        cout<<"Enter five coefficients: "<<endl;
        cin>>v[0];
        cin>>v[1];
        cin>>v[2];
        cin>>v[3];
        cin>>v[4];
        cout<<"The size of the vector is:  "<<v.size()<<endl;
        cout<<"Average: "<<mean(v)<<endl;
        cout<<"Standard deviation: "<<sdev(v)<<endl;
        sort(v);
        cout<<"That vector sorted is: "<<endl;
        for(int i=0;i<=5;i++) {
        cout<<v[i]<<endl;
        }
        cout<<"Enter x:"<<endl;
        cin>>x;
        cout<<"Enter y: "<<endl;
        cin>>y;
        std::swap(x,y);
        cout<<x<<endl;
        cout<<y<<endl;
        cout<<"Here is a vector with 100 uniform draws from the unit interval: "<<endl;
        randuniform(100,w);
        for(int i=0;i<100;i++) {
        cout<<w[i]<<endl;
        }
cout<<"What's the mean of these?"<<endl;
        cout<<mean(w);
        cout<<"Standard deviation is: "<<sdev(w)<<endl;
        cout<<"How about some normally distributed numbers? "<<endl;
        BM(100,t);
        for(int i=0;i<100;i++){
        cout<<t[i]<<endl;
        }
        cout<<"Their standard deviation is:"<<sdev(t)<<endl;
        cout<<"Mean is: "<<mean(t)<<endl;
        return 0;
        }



