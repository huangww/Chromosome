#include <iostream>


int main(int argc, char *argv[])
{
    int startTime = time(NULL);

    std::cout << "good" << std::endl;

    int endTime = time(NULL);
    int elapsedTime = endTime - startTime;
    std::cout << "Runing time: " << elapsedTime 
    << "seconds" << std::endl;
    return 0;
}
