#ifndef ASEP_HPP_MB6PNS13
#define ASEP_HPP_MB6PNS13

#include "project.hpp"

class Asep: public Project
{
public:
    Asep ();
    virtual ~Asep ();

    void print();
    void run();

private:
    class State *state;
};

#endif /* end of include guard: ASEP_HPP_MB6PNS13 */
