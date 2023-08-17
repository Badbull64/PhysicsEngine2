#pragma once 
#ifndef HASHGRID_H
#include <iostream>
#include <cmath>
#include <list>
#include <memory>
#include <vector>
#include <algorithm>
#include <cstdlib>

class c_HashTable {
public:
    float spacing;
    int tableSize;
    std::vector<int> cellStart;
    std::vector<int> cellEntries;
    std::vector<int> queryIds;
    int querySize;

    c_HashTable(float spacing, int maxNumObjects)
        : spacing(spacing),
        tableSize(2 * maxNumObjects),
        cellStart(tableSize + 1, 0),
        cellEntries(maxNumObjects, 0),
        queryIds(maxNumObjects, 0),
        querySize(0) {}

    int hashCoords(int xi, int yi) {
        int h = (xi * 92837111) ^ (yi * 689287499); // Modified fantasy function for 2D
        return std::abs(h) % tableSize;
    }

    int intCoord(float coord) {
        return std::floor(coord / spacing);
    }

    int hashPos(const std::vector<int>& pos, int nr) {
        return hashCoords(intCoord(pos[2 * nr]), intCoord(pos[2 * nr + 1]));
    }


    void create(const std::vector<int>& pos) {
        int numObjects = std::min(static_cast<int>(pos.size() / 2), static_cast<int>(cellEntries.size()));

        std::fill(cellStart.begin(), cellStart.end(), 0);
        std::fill(cellEntries.begin(), cellEntries.end(), 0);

        for (int i = 0; i < numObjects; i++) {
            int h = hashPos(pos, i);
            cellStart[h]++;
        }

        int start = 0;
        for (int i = 0; i < tableSize; i++) {
            start += cellStart[i];
            cellStart[i] = start;
        }
        cellStart[tableSize] = start;

        for (int i = 0; i < numObjects; i++) {
            int h = hashPos(pos, i);
            cellStart[h]--;
            cellEntries[cellStart[h]] = i;
        }
    }

    void query(const std::vector<int>& pos, int nr, int maxDist) {
        int x0 = intCoord(pos[2 * nr] - maxDist);
        int y0 = intCoord(pos[2 * nr + 1] - maxDist);

        int x1 = intCoord(pos[2 * nr] + maxDist);
        int y1 = intCoord(pos[2 * nr + 1] + maxDist);

        querySize = 0;

        for (int xi = x0; xi <= x1; xi++) {
            for (int yi = y0; yi <= y1; yi++) {
                int h = hashCoords(xi, yi);
                int start = cellStart[h];
                int end = cellStart[h + 1];

                for (int i = start; i < end; i++) {
                    queryIds[querySize] = cellEntries[i];
                    querySize++;
                }
            }
        }
        return;
    }
};

#endif // !HASHGRID_H
