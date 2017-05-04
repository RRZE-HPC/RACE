#ifndef NAME_LEVEL_DATA_H
#define NAME_LEVEL_DATA_H

#include "error.h"

struct LevelData{
    int *levelRow;
    int *levelNnz;
    int totalLevel;

    LevelData();
    ~LevelData();
};

/*LevelData::LevelData():levelRow(NULL),levelNnz(NULL),totalLevel(0)
{

}

LevelData::~LevelData()
{
    if(levelRow) {
        delete[] levelRow;
    }

    if(levelNnz) {
        delete[] levelNnz;
    }
}
*/
#endif
