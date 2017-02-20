#include<iostream>
#include"matplotlib-cpp-starter/matplotlibcpp.h"
using namespace std;

namespace plt = matplotlibcpp;

int main(){
    
    cout<<"matplotlib-cpp sample start"<<endl;
    
   // int n = ;
    double alpha=10;
    int c=-10;
    vector<double> x(21), y(21);
    for(int i=0; i<21; ++i) {
        x.at(i) = c;
        //if(c<5&&c>-5){
        y.at(i) = exp(-(pow(c,2)/pow(alpha,2)));
        //}
        //else y.at(i)=0;
        c++;
    }
    
    plt::plot(x, y, "--r");
    plt::show();
    
    return 0;
}
