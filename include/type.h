#ifndef RACE_TYPE_H
#define RACE_TYPE_H

#include <vector>

namespace RACE {

    enum dist{
        ONE=1,
        TWO=2,
        POWER=3
    };

    enum d2Method{
        TWO_BLOCK,
        THREE_BLOCK /*TODO: To be done*/
    };

    enum LBTarget{
        ROW,
        NNZ
    };

    enum PinMethod{
        FILL,
        SCATTER
    };

    /*Load balncing mode, based on efficiency or just minimise effRow*/
    enum LBMode{
        MIN,
        EFFICIENCY
    };

}

struct ZoneLeaf;

typedef std::vector<ZoneLeaf> tree_t;


#endif
