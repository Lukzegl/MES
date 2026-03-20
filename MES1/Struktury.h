#pragma once
#include <vector>
#include <iostream>
#include<iomanip>
#include <string>
#include "Gaus.h"

class Node
{
private:
    double x, y;
    bool BC;

public:
    Node(double x1=0, double y1=0);
    double get_x();
    double get_y();
    bool get_BC();
    void change_BC(bool state);
    void print();
};

class Element
{
private:
    std::vector<int> node_ids;

public:
    Element(int number, int ids[]);
    int get_id(int index);
    int size();
};

class Jakobian
{
public:
	Jakobian() = default;
    Jakobian(double dxi[], double deta[], std::vector<Node>* nodes, Element* element);
    double detJ;
    double J[2][2];
    double invJ[2][2];
};


class Grid
{
private:
    std::vector<Node> nodes;
    std::vector<Element> elements;

    int simulation_time;
    int simulation_step_time;
    int conductivity;
    int Alfa;
    int Tot;
    int InitialTemp;
    int Density;
    int SpecificHeat;
    int Nodes_number;
    int Elements_number;

public:
    void Create(std::string file);
    Grid();
    void ShowData();

    int get_simulation_time();
    int get_simulation_step_time();
    int get_conductivity();
    int get_Alfa();
    int get_Tot();
    int get_InitialTemp();
    int get_Density();
    int get_SpecificHeat();
    int get_Nodes_number();
    int get_Elements_number();

    std::vector<Node>* get_Nodes_w();
    std::vector<Element>* get_Elements_w();

    void printelement(int index)
    {
        std::cout << "Element " << index << ": ";
        for (int id = 0; id < elements[index].size(); id++)
        {
            std::cout << elements[index].get_id(id) << " ";
        }
        std::cout << std::endl;

        std::cout << "Node coordinates: " << std::endl;
        for (int id = 0; id < elements[index].size(); id++)
        {
            int node_id = elements[index].get_id(id);
            std::cout << "Node " << node_id << ": (" << nodes[node_id].get_x() << ", " << nodes[node_id].get_y() << ")" << std::endl;
        }
    }
};

class ElementData
{
private:
    std::vector<Node>* nodesLocal;

    int conductivity;
    int alfa;
    int tot;
    int density;
    int specificHeat;
    int gaussPointsN;

    double H_local[4][4];
    double Hbc[4][4];
    double C_local[4][4];
    double P[4];

    double dN_dxi[4];
    double dN_deta[4];

    std::vector<int> nodeIDs;
    Jakobian jakobian;
    

public:
    ElementData(Element* elem, std::vector<Node>* nodes, int cond, int alfa_val, int tot_val, int density_val, int specificHeat_val, int gaussN = 2);
    //ElementData(Element* elem, std::vector<Node>* nodes);

    void calculate_Hbc();
	void calculate_Clocal();

	const double* getH() const;
    const double* getHbc() const;
    
    const double* getP() const;
    const double* getC() const;
    
	void printAllData() const;

    Element* element;
	const std::vector<int>& getNodeIDs() const;
};


class GlobalMatrix 
{
private:
    std::vector<std::vector<double>> H_global;
	std::vector<std::vector<double>> C_global;
	std::vector<double> P_global;
    int size;

public:
    GlobalMatrix(int n);

    void addLocalMatrix(const ElementData& elemData);
    void addLocalVector(const ElementData& elemData);
	void addLocalC(const ElementData& elemData);


    void printGlobalH();
	void printGlobalC();

	std::vector<std::vector<double>> getGlobalH() const { return H_global; }
	std::vector<std::vector<double>> getGlobalC() const { return C_global; }
	const std::vector<double>& getGlobalP() const { return P_global; }

};
