#include "levelData.h"

LevelData::LevelData():levelRow(NULL),levelNnz(NULL),totalLevel(0)
{

}

//assignment operator
LevelData& LevelData::operator=(const LevelData& other)
{
	if(this!=&other)
	{
		if(totalLevel != 0)
		{
			ERROR_PRINT("Cannot assign levelData")
		}
		totalLevel = other.totalLevel;
		levelRow = new int[totalLevel];
		levelNnz = new int[totalLevel];
		for(int i=0; i<totalLevel; ++i)
		{
			levelRow[i] = other.levelRow[i];
			levelNnz[i] = other.levelNnz[i];
		}
	}

	return *this;
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


