#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <fstream>
#include <chrono>

class Timer
{
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::time_point<std::chrono::system_clock> end;
    double duration;
    
public:
    Timer() 
    {
        start = std::chrono::system_clock::now();
        std::ofstream time_out;
        time_out.open("time_out.txt", std::ios_base::app);
        time_out << "Start -- ";
        time_out.close();
    }
    ~Timer() 
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        duration = elapsed.count();
        std::ofstream time_out;
        time_out.open("time_out.txt", std::ios_base::app);
        time_out << "Elapsed Time: " << duration << "s\n";
        time_out.close();
    }
};

#endif
