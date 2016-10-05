#ifndef INPUT_HPP_DRBPSRM7
#define INPUT_HPP_DRBPSRM7

#include <iostream>
#include <map>
#include "define.hpp"

class Input
{
public:
    Input ();
    virtual ~Input ();

    void getInput(int argc, char *argv[]);
    void file();
    void print();
    
    std::string projectName;
    std::map<std::string, double> parameter;

private:
    std::string infname;

    void parse(std::string line);
};

#endif /* end of include guard: INPUT_HPP_DRBPSRM7 */
