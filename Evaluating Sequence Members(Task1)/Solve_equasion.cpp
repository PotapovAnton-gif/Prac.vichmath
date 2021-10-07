#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;



int main() {

    long double a[50];

    double pi = M_PI;
    double e = M_E;
    double arctg = atan(1);


    a[0] = pi / (arctg + 1 - 5 * e);
    a[1] = a[0];

    for (int i = 2; i < 50; i++) {
        a[i] = (pi - a[i-1] + 5 * e * a[i-2]) / arctg;  
    }
    int stchetchik = 0;
    for (int i = 0; i < 50; i++) {
        
        cout << setprecision(16) << fixed;
        cout << a[i] << endl;

        stchetchik++;
    }
    cout << stchetchik << endl;
    
    return 0;
}




