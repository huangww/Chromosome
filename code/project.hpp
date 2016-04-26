#ifndef PROJECT_HPP_SVUCG7LM
#define PROJECT_HPP_SVUCG7LM

#include "input.hpp"

class Project
{
public:
    Project () { };
    virtual ~Project (){ };
    
    virtual void setup(Input *input) = 0;
    virtual void run() = 0;
};

#endif /* end of include guard: PROJECT_HPP_SVUCG7LM */
