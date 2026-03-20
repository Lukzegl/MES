#include "Struktury.h"
#include "Wczytywanie.h"
#include "Gaus.h"
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
//obiektowa nuda start
Node::Node(double x1, double y1)
{
    x = x1;
    y = y1;
}

double Node::get_x()
{
    return x;
}

double Node::get_y()
{
    return y;
}

bool Node::get_BC()
{
    return BC;
}

void Node::change_BC(bool state)
{
    BC = state;
}

void Node::print()
{
    std::cout << "x: " << x << " Y: " << y << " BC: " << BC << std::endl;
}

Element::Element(int number, int ids[])
{
    for (int i = 0; i < number; i++)
    {
        node_ids.push_back(ids[i]);
    }
}

int Element::get_id(int index)
{
    return node_ids[index];
}

int Element::size()
{
    return node_ids.size();
}

void Grid::Create(std::string file)
{
    read_data(&simulation_time, &simulation_step_time, &conductivity, &Alfa, &Tot, &InitialTemp, &Density, &SpecificHeat, &Nodes_number, &Elements_number, file, &nodes, &elements);
}
//obiektowa nuda end


Grid::Grid()
{
    simulation_time = 0;
    simulation_step_time = 0;
    conductivity = 0;
    Alfa = 0;
    Tot = 0;
    InitialTemp = 0;
    Density = 0;
    SpecificHeat = 0;
    Nodes_number = 0;
    Elements_number = 0;
}

void Grid::ShowData()
{
    std::cout << "Simulation time: " << simulation_time << std::endl;
    std::cout << "Simulation step time: " << simulation_step_time << std::endl;
    std::cout << "Conductivity: " << conductivity << std::endl;
    std::cout << "Alfa: " << Alfa << std::endl;
    std::cout << "Tot: " << Tot << std::endl;
    std::cout << "Initial Temp: " << InitialTemp << std::endl;
    std::cout << "Density: " << Density << std::endl;
    std::cout << "Specific Heat: " << SpecificHeat << std::endl;
    std::cout << "Nodes number: " << Nodes_number << std::endl;
    std::cout << "Elements number: " << Elements_number << std::endl;

    if (nodes.size() > 0)
    {
        std::cout << std::endl << "NODES: " << std::endl;
        for (int i = 0; i < nodes.size(); i++)
        {
            std::cout << i << ". " << nodes[i].get_x() << " " << nodes[i].get_y() << std::endl;
        }
    }
    else
    {
        std::cout << "No nodes data" << std::endl;
    }

    if (elements.size() > 0)
    {
        std::cout << std::endl << "ELEMENTS: " << std::endl;
        for (int i = 0; i < elements.size(); i++)
        {
            for (int id = 0; id < elements[i].size(); id++)
            {
                std::cout << elements[i].get_id(id) << " ";
            }
            std::cout << std::endl;
        }
    }
    else
    {
        std::cout << "No elements data" << std::endl;
    }
}
//obiektowa nuda1 start
int Grid::get_simulation_time()
{
    return simulation_time;
}

int Grid::get_simulation_step_time()
{
    return simulation_step_time;
}

int Grid::get_conductivity()
{
    return conductivity;
}

int Grid::get_Alfa()
{
    return Alfa;
}

int Grid::get_Tot()
{
    return Tot;
}

int Grid::get_InitialTemp()
{
    return InitialTemp;
}

int Grid::get_Density()
{
    return Density;
}

int Grid::get_SpecificHeat()
{
    return SpecificHeat;
}

int Grid::get_Nodes_number()
{
    return Nodes_number;
}

int Grid::get_Elements_number()
{
    return Elements_number;
}

std::vector<Node>* Grid::get_Nodes_w()
{
    return &nodes;
}

std::vector<Element>* Grid::get_Elements_w()
{
    return &elements;
}

//obiektowa nuda end

//ElementData::ElementData(Element* elem, std::vector<Node>* nodes) : element(elem), nodesLocal(nodes)
ElementData::ElementData(Element* elem, std::vector<Node>* nodes, int cond, int alfa_val, int tot_val, int density_val, int specificHeat_val, int gaussN):element(elem), nodesLocal(nodes), conductivity(cond), alfa(alfa_val), tot(tot_val), density(density_val), specificHeat(specificHeat_val), gaussPointsN(gaussN)
{
    GausPoint gausPoints;

    //inicjalizacja
    for (int i = 0; i < 4; i++)
    {
        P[i] = 0.0;

        for (int j = 0; j < 4; j++)
        {
            H_local[i][j] = 0.0;
            Hbc[i][j] = 0.0;
			C_local[i][j] = 0.0;
        }
    }

    //id nodes
    for(int i=0;i<4;i++)
    {
        nodeIDs.push_back(element->get_id(i));
	}
    
    //SCHEMAT 2 3 4
	double* xArr;
	double* wArr;

	int n = gaussPointsN;
    switch (n)
    {
        case 2:
			xArr = GausPoint::x2;
            wArr = GausPoint::w2;
            break;
		case 3:
            xArr = GausPoint::x3;
            wArr = GausPoint::w3;
			break;
        case 4:
            xArr = GausPoint::x4;
			wArr = GausPoint::w4;
			break;
    default:
		cout << "========================ERRROR====================" << endl;
        cout << "====================USTAWIONO NA 2================" << endl;
		xArr = GausPoint::x2;
		wArr = GausPoint::w2;
        break;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double ksi = xArr[i];
            double eta = xArr[j];
			double weight = wArr[i] * wArr[j];

            dN_dxi[0] = -0.25 * (1 - eta);
            dN_dxi[1] = 0.25 * (1 - eta);
            dN_dxi[2] = 0.25 * (1 + eta);
            dN_dxi[3] = -0.25 * (1 + eta);

            dN_deta[0] = -0.25 * (1 - ksi);
            dN_deta[1] = -0.25 * (1 + ksi);
            dN_deta[2] = 0.25 * (1 + ksi);
            dN_deta[3] = 0.25 * (1 - ksi);

			jakobian = Jakobian(dN_dxi, dN_deta, nodesLocal, element);
           
            //H
            double dN_dx[4];
            double dN_dy[4];

            for (int k = 0; k < 4; k++)
            {
                dN_dx[k] = jakobian.invJ[0][0] * dN_dxi[k] + jakobian.invJ[0][1] * dN_deta[k];
                dN_dy[k] = jakobian.invJ[1][0] * dN_dxi[k] + jakobian.invJ[1][1] * dN_deta[k];
            }

            for (int a = 0; a < 4; a++)
            {
                for (int b = 0; b < 4; b++)
                {
                    H_local[a][b] += conductivity * (dN_dx[a] * dN_dx[b] + dN_dy[a] * dN_dy[b]) * jakobian.detJ * weight;
                }
	            }
            //
        }
    }

	calculate_Hbc();

	//HBC + LOCAL
    for(int i =0; i <4; i++)
    {
        for(int j=0;j<4;j++)
        {
            H_local[i][j] += Hbc[i][j];
        }
	}
    
	calculate_Clocal();
}


void ElementData::calculate_Hbc()
{
    double HBC_L[4][4] = { 0 };

    double P_l[4] = { 0 };

    
    double* xArr;
    double* wArr;
    int n = gaussPointsN;

    switch (n)
    {
    case 2:
        xArr = GausPoint::x2;
        wArr = GausPoint::w2;
        break;
    case 3:
        xArr = GausPoint::x3;
        wArr = GausPoint::w3;
        break;
    case 4:
        xArr = GausPoint::x4;
        wArr = GausPoint::w4;
        break;
    default:
        std::cout << "========================ERRROR====================" << std::endl;
        std::cout << "====================USTAWIONO NA 2================" << std::endl;
        xArr = GausPoint::x2;
        wArr = GausPoint::w2;
        n = 2;
        break;
    }


    for (int edge = 0; edge < 4; ++edge)
    {
        int idA = nodeIDs[edge];
        int idB = (edge < 3) ? nodeIDs[edge + 1] : nodeIDs[0];

        if (nodesLocal->at(idA).get_BC() == true && nodesLocal->at(idB).get_BC() == true)
        {
            
            double xA = nodesLocal->at(idA).get_x();
            double yA = nodesLocal->at(idA).get_y();
            double xB = nodesLocal->at(idB).get_x();
            double yB = nodesLocal->at(idB).get_y();

            double edgeLength = sqrt(pow(xB - xA, 2) + pow(yB - yA, 2));

           
            for (int gp = 0; gp < n; ++gp)
            {
                double s = xArr[gp];
                double weight = wArr[gp];

                double ksi = 0.0;
                double eta = 0.0;

                switch (edge)
                {
                case 0:
                    ksi = s; eta = -1.0; //DÓ£
                    break;
                case 1:
                    ksi = 1.0; eta = s; //LEWO
                    break;
                case 2:
                    ksi = -s; eta = 1.0; //GÓRA
                    break;
                case 3:
                    ksi = -1.0; eta = -s;//PRAWO
                    break;
                }

                double N[4];
                N[0] = 0.25 * (1 - ksi) * (1 - eta);
                N[1] = 0.25 * (1 + ksi) * (1 - eta);
                N[2] = 0.25 * (1 + ksi) * (1 + eta);
                N[3] = 0.25 * (1 - ksi) * (1 + eta);

                //HBC + P
                for (int a = 0; a < 4; ++a)
                {
                    for (int b = 0; b < 4; ++b)
                    {
                        HBC_L[a][b] += N[a] * N[b] * weight * (edgeLength / 2.0);
                    }
                    P_l[a] += N[a] * weight * (edgeLength / 2.0);
                }
            }
        }
    }

    
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Hbc[i][j] = HBC_L[i][j] * alfa;
        }
    }

    for (int i = 0; i < 4; i++)
    {
        P[i] = P_l[i] * tot * alfa;
    }
}

Jakobian::Jakobian(double dxi[], double deta[], std::vector<Node>* nodes, Element* element)
{
    int node0 = element->get_id(0);
    int node1 = element->get_id(1);
    int node2 = element->get_id(2);
    int node3 = element->get_id(3);

    J[0][0] = 
        dxi[0] * nodes->at(node0).get_x() +
        dxi[1] * nodes->at(node1).get_x() +
        dxi[2] * nodes->at(node2).get_x() +
        dxi[3] * nodes->at(node3).get_x();

    J[1][0] = 
        deta[0] * nodes->at(node0).get_x() +
        deta[1] * nodes->at(node1).get_x() +
        deta[2] * nodes->at(node2).get_x() +
        deta[3] * nodes->at(node3).get_x();

    J[0][1] = 
        dxi[0] * nodes->at(node0).get_y() +
        dxi[1] * nodes->at(node1).get_y() +
        dxi[2] * nodes->at(node2).get_y() +
        dxi[3] * nodes->at(node3).get_y();

    J[1][1] = 
        deta[0] * nodes->at(node0).get_y() +
        deta[1] * nodes->at(node1).get_y() +
        deta[2] * nodes->at(node2).get_y() +
        deta[3] * nodes->at(node3).get_y();

    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    
    
    invJ[0][0] = J[1][1] / detJ;
    invJ[0][1] = -J[0][1] / detJ;
    invJ[1][0] = -J[1][0] / detJ;
    invJ[1][1] = J[0][0] / detJ;
}

GlobalMatrix::GlobalMatrix(int n) : size(n) 
{
    H_global.resize(n, std::vector<double>(n, 0.0));
	C_global.resize(n, std::vector<double>(n, 0.0));
	P_global.resize(n, 0.0);
}

//H
void GlobalMatrix::addLocalMatrix(const ElementData& elemData) 
{
	const std::vector<int>& nodes = elemData.getNodeIDs();
	const double* H_local = elemData.getH();

    for(int i =0;i <4; i++) 
    {
        for(int j = 0; j < 4; j++) 
        {
            H_global[nodes[i]][nodes[j]] += H_local[i * 4 + j];
        }
	}
}

//C
void GlobalMatrix::addLocalC(const ElementData& elemData)
{
    const std::vector<int>& nodes = elemData.getNodeIDs();
    const double* C_local = elemData.getC();

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            C_global[nodes[i]][nodes[j]] += C_local[i * 4 + j];
        }
    }
}

void GlobalMatrix::addLocalVector(const ElementData& elemData)
{
    const std::vector<int>& nodes = elemData.getNodeIDs();
    const double* P_local = elemData.getP();

    for (int i = 0; i < 4; i++) 
    {
        P_global[nodes[i]] += P_local[i];
    }
}
//C
void ElementData::calculate_Clocal()
{
	GausPoint gausPoints;
	double* xArr;
	double* wArr;
	int n = gaussPointsN;

    switch (n)
    {
    case 2:
        xArr = GausPoint::x2;
        wArr = GausPoint::w2;
        break;
    case 3:
        xArr = GausPoint::x3;
        wArr = GausPoint::w3;
        break;
    case 4:
        xArr = GausPoint::x4;
        wArr = GausPoint::w4;
        break;
    default:
        cout << "========================ERRROR====================" << endl;
        cout << "====================USTAWIONO NA 2================" << endl;
        xArr = GausPoint::x2;
        wArr = GausPoint::w2;
        break;
    }

	double N[4];
	double C_temp[4][4] = { 0 };

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double ksi = xArr[i];
            double eta = xArr[j];
            double weight = wArr[i] * wArr[j];

            N[0] = 0.25 * (1 - ksi) * (1 - eta);
            N[1] = 0.25 * (1 + ksi) * (1 - eta);
            N[2] = 0.25 * (1 + ksi) * (1 + eta);
            N[3] = 0.25 * (1 - ksi) * (1 + eta);

			// TAK TRZEBA OD NOWA LICZYC JAKOBIANA BO INACZEJ NIE DZIA£A
            //Kwesia osobnych funkcjii podzia³u na segmenty
            dN_dxi[0] = -0.25 * (1 - eta);
            dN_dxi[1] = 0.25 * (1 - eta);
            dN_dxi[2] = 0.25 * (1 + eta);
            dN_dxi[3] = -0.25 * (1 + eta);

            dN_deta[0] = -0.25 * (1 - ksi);
            dN_deta[1] = -0.25 * (1 + ksi);
            dN_deta[2] = 0.25 * (1 + ksi);
            dN_deta[3] = 0.25 * (1 - ksi);

            jakobian = Jakobian(dN_dxi, dN_deta, nodesLocal, element);
            for (int a = 0; a < 4; a++)
            {
                for (int b = 0; b < 4; b++)
                {
                    C_temp[a][b] += density * specificHeat * N[a] * N[b] * jakobian.detJ * weight;
                }
            }
        }
	}

    //CL
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            C_local[i][j] = C_temp[i][j];
        }
	}

   
}


//PRINT---------------------------------------------------------------------------------------------------------

void GlobalMatrix::printGlobalH() 
{
    std::cout << "Global H Matrix:" << std::endl;
    for (int i = 0; i < size; i++) 
    {
        for (int j = 0; j < size; j++) 
        {
            std::cout << H_global[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void GlobalMatrix::printGlobalC()
{
    std::cout << "Global C Matrix:" << std::endl;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            std::cout << C_global[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void ElementData::printAllData() const
{
    std::cout << std::endl << "Element Data:" << std::endl;
    std::cout << "Node IDs: ";
    for (int i = 0; i < element->size(); i++)
    {
        std::cout << element->get_id(i) << " ";
    }
    std::cout << std::endl;
    std::cout << std::endl << "H Local Matrix:" << std::endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << H_local[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << "Hbc Matrix:" << std::endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << Hbc[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout<<std::endl << "P Vector:" << std::endl;
    for (int i = 0; i < 4; i++)
    {
        std::cout << P[i] << " ";
    }
    std::cout << std::endl;
    std::cout<<std::endl << "C Local Matrix:" << std::endl;
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << C_local[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


//END PRINT---------------------------------------------------------------------------------------------------------

//obietowa nuda

const double* getCptr(const ElementData& ed) 
{
    return &(reinterpret_cast<const ElementData&>(ed).getH()[0]);
}


const double* ElementData::getH() const 
{
    return &H_local[0][0];
}

const double* ElementData::getHbc() const 
{
    return &Hbc[0][0];
}

const double* ElementData::getP() const 
{
    return &P[0];
}

const double* ElementData::getC() const 
{
    return &C_local[0][0];
}

const std::vector<int>& ElementData::getNodeIDs() const 
{
    return nodeIDs;
}
