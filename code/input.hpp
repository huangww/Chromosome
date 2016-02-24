#ifndef INPUT_HPP_DRBPSRM7
#define INPUT_HPP_DRBPSRM7

#include <iostream>

class Input
{
public:
    Input ();
    virtual ~Input ();

    void getInput(int argc, char *argv[]);
    void file();
    
    struct Parameter {
        std::string *paraName;
        double *paraValue;
    } *parameter;

private:
    std::string inputfname;
    int nPara;      // number of parameters

    void init();
    void print();
    void parse(int index, std::string line);
    void excute();
    int getParaNumber();
};

#endif /* end of include guard: INPUT_HPP_DRBPSRM7 */
