#include <stdlib.h>
#include <time.h>
#include <iostream>
using namespace std;

int main()
{
    srand(time(NULL));
    int M, N;
    cin >> M >> N;
    cout << M << " " << N << endl;

    for (int i = 0; i < M; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            cout << rand() % 2 << " ";
        }
        cout << endl;
    }
    return 0;
}