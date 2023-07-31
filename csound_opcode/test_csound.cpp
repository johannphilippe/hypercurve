//#include<csound.hpp>
#include<csound/csound.hpp>
#include<iostream>
#include<string>

int main()
{
    Csound cs;
    //cs.SetOption("--verbose");
    cs.SetOption("--opcode-lib=./libcsound_hypercurve.so");
    std::string csd_path = "test.csd";
    int res = cs.Compile(csd_path.c_str());
    if(res != 0)
    {
        int c = cs.GetMessageCnt();
        if(c > 0) {
            for(int i = 0; i < c; i++)
            {
                std::cout << cs.GetFirstMessage() << std::endl;
                cs.PopFirstMessage();
            }
        }
        throw(std::runtime_error("res is not 0 : " + std::to_string(res)));
    }
    cs.Perform();
    return 0;
}
