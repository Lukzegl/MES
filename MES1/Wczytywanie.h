#pragma once

#include <string>
#include <vector>
#include "Struktury.h"

void read_data(
    int* simulation_time,
    int* simulation_step_time,
    int* conductivity,
    int* Alfa,
    int* Tot,
    int* InitialTemp,
    int* Density,
    int* SpecificHeat,
    int* Nodes_number,
    int* Elements_number,
    std::string filename,
    std::vector<Node>* Nodes,
    std::vector<Element>* Elements
);