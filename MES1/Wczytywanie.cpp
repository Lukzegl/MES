#include "wczytywanie.h"
#include "Struktury.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;

void read_data(int* simulation_time, int* simulation_step_time, int* conductivity, int* Alfa, int* Tot, int* InitialTemp, int* Density, int* SpecificHeat, int* Nodes_number, int* Elements_number, string filename, vector<Node>* Nodes, vector<Element>* Elements)
{
	ifstream file(filename);

	if (!file.is_open())
	{
		cout << "BRAK PLIKU!!!!" << endl;
		return;
	}
	//1 linijka
	string line;
	getline(file, line);
	int pos = line.find(" ");
	*simulation_time = stoi(line.substr(pos + 1));

	//2 linijka
	getline(file, line);
	pos = line.find(" ");
	*simulation_step_time = stoi(line.substr(pos + 1));

	//3 linijka
	getline(file, line);
	pos = line.find(" ");
	*conductivity = stoi(line.substr(pos + 1));

	//4 linijka
	getline(file, line);
	pos = line.find(" ");
	*Alfa = stoi(line.substr(pos + 1));

	//5 linijka
	getline(file, line);
	pos = line.find(" ");
	*Tot = stoi(line.substr(pos + 1));

	//6 linijka
	getline(file, line);
	pos = line.find(" ");
	*InitialTemp = stoi(line.substr(pos + 1));

	//7 linijka
	getline(file, line);
	pos = line.find(" ");
	*Density = stoi(line.substr(pos + 1));

	//8 linijka
	getline(file, line);
	pos = line.find(" ");
	*SpecificHeat = stoi(line.substr(pos + 1));

	//9 linijka
	getline(file, line);
	pos = line.find(" ");
	line = line.substr(pos + 1);
	pos = line.find(" ");
	*Nodes_number = stoi(line.substr(pos + 1));

	//10 linijka
	getline(file, line);
	pos = line.find(" ");
	line = line.substr(pos + 1);
	pos = line.find(" ");
	*Elements_number = stoi(line.substr(pos + 1));

	getline(file, line); //*Nodes

	for (int i = 0; i < *Nodes_number; i++)
	{
		getline(file, line);

		pos = line.find(",");
		line = line.substr(pos + 1);

		pos = line.find(",");

		double x = stod(line.substr(0, pos));
		double y = stod(line.substr(pos + 1));

		Node node(x, y);
		Nodes->push_back(node);
	}

	getline(file, line); //*Elements

	for (int i = 0; i < *Elements_number; i++)
	{
		getline(file, line);
		pos = line.find(",");
		line = line.substr(pos + 1);
		pos = line.find(",");

		int a = stoi(line.substr(0, pos));

		line = line.substr(pos + 1);
		pos = line.find(",");

		int b = stoi(line.substr(0, pos));

		line = line.substr(pos + 1);
		pos = line.find(",");

		int c = stoi(line.substr(0, pos));
		int d = stoi(line.substr(pos + 1));

		a = a - 1;
		b = b - 1;
		c = c - 1;
		d = d - 1;
		int ids[4] = { a, b, c, d };

		Element element(4, ids);
		Elements->push_back(element);
	}

	//BC
	getline(file, line);

	getline(file, line);
	//cout<<endl << "AAAA " << line<<endl;
	string nodeid = "";
	int ID_i;

	

	for (int i = 0; i < line.length()+1; i++)
	{
		if (line[i] == ',')
		{
			ID_i = stoi(nodeid);
			Nodes->at(ID_i-1).change_BC(true);
			//cout << ID_i - 1 << " ";
			nodeid = "";
		}
		else
		{
			nodeid = nodeid + line[i];
		}
	}
	ID_i = stoi(nodeid);
	Nodes->at(ID_i - 1).change_BC(true);

	/*
	cout << "TEST_____BCCC" << endl;

	for (int i = 0; i < *Nodes_number; i++)
	{
		cout << i << " ";
		Nodes->at(i).print();
	}
	*/
	file.close();
}