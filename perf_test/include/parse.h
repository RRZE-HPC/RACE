#ifndef RACE_PARSE_H
#define RACE_PARSE_H

#include "type.h"
#include <vector>
#include <getopt.h>

using namespace RACE;

struct my_option
{
    option gnu_opt;
    char const* desc;
    my_option(const char*name, int has_arg, int *flag, int val, char const* desc);
    my_option();
};

struct parser
{
        char *mat_file;
        int iter;
        int cores;
        int smt;
        int nodes;
        PinMethod pin;
        bool validate;
        double tol;
        char const* prgname;
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
