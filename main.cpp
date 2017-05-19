#include <iostream>

#include <armadillo>
#include "params.h"
#include "functions.h"
#include <ctime>

using namespace std;

int main(int argc, char *argv[])
{
    cout << H(2) << endl;
    cout << H(-.5)<< endl;
    cout << sum_elem(3,4,3.14/10,3.14/9,0.2,.453) << endl;
    cout << sum_elem(-3,4,3.14/10,3.14/9,0.2,.453) << endl;

    clock_t begin = clock();
//    for (int i=0; i<1000;++i)
//        cout <<
//                sum_integr(3.14/10,3.14/9,0.2,.453);
//                        <<endl;
    clock_t end = clock();
    cout << "time: "<< double(end - begin) / CLOCKS_PER_SEC << endl;
    cout << "result: "<< sum_integr(3.14/10,3.14/9,0.2,.453)<<endl;

    //    cout << integr_s(3.14/10,3.14/9,0.2,0.15)<<endl;
    return 0;
}
