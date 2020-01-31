#include <stdlib.h>
#include <stdio.h>
#include <solver.h>
#include "global.h"
#include "mpi.h"
#include <float.h>
#include "decomp.h"

#include "SolverZeidel.h"

static void show_usage(std::string name)
{
    std::cout << "Usage: " << name << " [<option(s)>]\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-d,--decomp PROCESSORS\tDecomposition of the computational domain\n";
}


int main(int argc, char** argv)
{
#ifdef _DEBUG
	_controlfp(~(_MCW_EM & (~_EM_INEXACT) & (~_EM_UNDERFLOW)), _MCW_EM);
#endif

    Parallel::init(&argc, &argv);

    hLog = fopen("task.log", "w"); // открываем файл для записи лога; вместо printf(...) необходимо использовать log(...)

	if (argc > 1) {
        std::vector <std::string> sources;
        std::string str;
        int procCount;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if ((arg == "-h") || (arg == "--help")) {
                show_usage(argv[0]);
                return 0;
            } else if ((arg == "-d") || (arg == "--decomp")) {
                if (i + 1 < argc) {
                    str = argv[++i];
                    try {
                        procCount = std::stoi(str);
                    }
                    catch(...) {
                        std::cerr << "--decomp option requires integer argument." << std::endl;
                        show_usage(argv[0]);
                        return 1;
                    }
                } else {
                    std::cerr << "--decomp option requires one argument." << std::endl;
                    show_usage(argv[0]);
                    return 1;
                }
            } else {
                std::cerr << "Wrong arguments." << std::endl;
                show_usage(argv[0]);
                return 1;
            }
        }

        Method* m = new Decomp(procCount);
        m->init((char*) "task.xml");
        m->run();
        m->done();
        delete m;
        return 0;
	}
	else {
        Method *method = Solver::initMethod((char *) "task.xml");
        Solver::runMethod(method);
        Solver::destroyMethod(method);
    }

	fclose(hLog);

	Parallel::done();

	return 0;
}