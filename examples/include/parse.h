#ifndef RACE_PARSE_H
#define RACE_PARSE_H

#include <RACE/type.h>
#include <vector>
#include <getopt.h>

using namespace RACE;

struct my_option
{
    option gnu_opt;
    char* desc;
    my_option(const char*name, int has_arg, int *flag, int val, char* desc);
    my_option();
};

struct parser
{
        char *mat_file;
        int iter;
        std::vector<int> iter_vec;
        int cores;
        int smt;
        int nodes;
        double cache_size;
        std::vector<int> cache_size_vec;
        PinMethod pin;
        bool validate;
        double tol;
        double convTol;
        bool RCM_flag;
        char *colorType;
        char *prgname;
        int numOptions;
        my_option *long_options;
        option *gnuOptions;
        parser();
        ~parser();
        bool parse_arg(int argc, char **argv);
        int dump_arg();
        void help();
};

#endif
