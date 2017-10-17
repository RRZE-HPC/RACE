#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parse.h"

my_option::my_option(const char* name_, int has_arg_, int* flag_, int val_, char*desc_): desc(desc_)
{
    gnu_opt = {name_, has_arg_, flag_, val_};
}

my_option::my_option()
{
}

parser::parser():mat_file(NULL), iter(100), cores(1), smt(1), pin(FILL), prgname("a.out"), numOptions(6)
{
    long_options = new my_option[numOptions+1];

    long_options[0] = {"matrix",  required_argument, 0,  'm', "Matrix File in MatrixMarket Format"};
    long_options[1] = {"iter",    required_argument, 0,  'i', "Iterations to be carried out" };
    long_options[2] = {"cores",   required_argument, 0,  'c', "Number of cores to be used" };
    long_options[3] = {"smt",     required_argument, 0,  't', "Number of threads per core to be used (recommended 1)" };
    long_options[4] = {"pin",     required_argument, 0,  'p', "Pinning strategy to be used; availablle options FILL or SCATTER" };
    long_options[5] = {"help",    no_argument,       0,  'h', "Prints this help informations" };
    long_options[6] = {0,         0,                 0,   0,  0 };

    gnuOptions = new option[numOptions+1];

    for(int i=0; i<numOptions; ++i)
    {
        gnuOptions[i] = long_options[i].gnu_opt;
    }

}

parser::~parser()
{
    delete[] long_options;
    delete[] gnuOptions;
}

bool parser::parse_arg(int argc, char **argv)
{
    prgname = argv[0];
    while (1) {
        int option_index = 0, c;
        c = getopt_long(argc, argv, "0:m:i:c:t:p:h",
                gnuOptions, &option_index);

        printf("c=%d\n", c);
        if (c == -1)
            break;

        switch (c) {
            case 0:
                printf("No parameters specified.\n Usage: %s -m [MATRIX] -i [ITERATIONS] -c [CORES] -t [SMT] -p [FILL/SCATTER]\n", prgname);
                return false;
                break;
            case 'm':
                mat_file = optarg;
                break;
            case 'i':
                iter = atoi(optarg);
                break;
            case 'c':
                printf("cores = %s\n", optarg);
                cores = atoi(optarg);
                break;
            case 't':
                smt = atoi(optarg);
                break;
           case 'p':
                char *pin_str = optarg;
                printf("pin strategy = %s\n", pin_str);
                if(!(strcmp(pin_str, "SCATTER")))
                {
                    pin = SCATTER;
                }
                else if(!(strcmp(pin_str, "FILL")))
                {
                    pin = FILL;
                }
                else
                {
                    printf("%s pin strategy unknown. setting to FILL\n", pin_str);
                    pin = FILL;
                }
                break;
           case 'h':
                help();
                return 0;
                break;
           default:
                break;
        }
    }

    if(mat_file == NULL)
    {
        return false;
    }
    else
    {
        return true;
    }
}

void parser::help()
{
    printf("Usage: %s [OPTION]...\n",prgname);
    printf("Valid options are:\n\n");

    printf("   %s    \t %s\n", "option", "description");
    for(int i=0; i<numOptions; ++i)
    {
        printf("-%c or --%s \t %s\n", ((char) long_options[i].gnu_opt.val), long_options[i].gnu_opt.name, long_options[i].desc);
    }
    /*    printf(" -m, --matrix=MATRIX FILE\t\tMatrix File in MatrixMarket Format\n");
    printf(" -c, --cores=CORES\t\tNumber of cores to be used\n");
    printf(" -t, --smt=THREADS PER CORE\t\tNumber of threads per core to be used (recommended 1)\n");
    printf(" -p, --pin=PIN STRATEGY\t\tPinning strategy to be used; availablle options FILL or SCATTER\n");
    printf(" -h, --help \t\t Prints this help informations\n");
    */
}
