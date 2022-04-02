#include<csound/csound.hpp>


int main()
{
    Csound cs;
    cs.CompileCsd("other_test.csd");
    cs.Start();
    cs.Perform();

    return 0;
}
