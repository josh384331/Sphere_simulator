#include "arrow.h"
#include <chrono>

int main(){
    cout << "Program Started" << endl;
    auto start_time = chrono::high_resolution_clock::now();
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end_time - start_time;
    arrow myArrow("arrow.json");
    start_time = chrono::high_resolution_clock::now();
    myArrow.exercise_3_9();
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    printf("Arrow sim w/quat [sec]: %f\n",duration.count());
    cout << "program Ended" << endl;
    return 0;
}

