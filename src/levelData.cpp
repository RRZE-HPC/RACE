#include "levelData.h"

LevelData::LevelData():levelRow(NULL),levelNnz(NULL),totalLevel(0)
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


