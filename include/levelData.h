#ifndef RACE_LEVEL_DATA_H
#define RACE_LEVEL_DATA_H

#include "error.h"

struct LevelData{
    int *levelRow;
    int *levelNnz;
    int totalLevel;
    int nrow;
    int nnz;

    LevelData();
    LevelData& operator=(const LevelData& other);
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
