#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include "Struktury.h"
#include "Wczytywanie.h"
#include "Gaus.h"
using namespace std;

int main()
{
    Grid grid;
    grid.Create("dane1.txt");

    grid.ShowData();
    cout << endl;

    int gaussN = 4; //2,3,4
    
    std::string output_filename = "wyniki_temperatury.txt";
    

    GlobalMatrix globalH(grid.get_Nodes_number());

    for (int i = 0; i < grid.get_Elements_number(); i++)
    {
        //
		// MIEJSCE NA STAWIANIE ELEMENTU
		// POD DZIWNE RZECZY Z KONDUKTYWNOSCIA I INNYMI PARAMETRAMI

        
        //
        //

        ElementData elemData(&grid.get_Elements_w()->at(i),
            grid.get_Nodes_w(),
            grid.get_conductivity(),
            grid.get_Alfa(),
            grid.get_Tot(),
            grid.get_Density(),
            grid.get_SpecificHeat(),
            gaussN);

        globalH.addLocalMatrix(elemData);
        globalH.addLocalVector(elemData);
        globalH.addLocalC(elemData);

        std::cout << "element " << i << ":" << std::endl;
		elemData.printAllData();

    }
    
    std::cout << "\nGlobalna macierz H:" << std::endl;
    globalH.printGlobalH();

    std::cout << "\nGlobalny wektor P:" << std::endl;
    const std::vector<double>& P_global = globalH.getGlobalP();
    for (int i = 0; i < P_global.size(); i++) {
        std::cout << P_global[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "\nGlobalna macierz C:" << std::endl;
    globalH.printGlobalC();
   
    int simulation_time = grid.get_simulation_time();
    int time_step = grid.get_simulation_step_time();
    int initial_temp = grid.get_InitialTemp();
    int n_nodes = grid.get_Nodes_number();

    std::vector<double> t_current(n_nodes, initial_temp);  
    std::vector<double> t_previous(n_nodes, initial_temp); 
    double delta_tau = time_step;

    //Zapis////////////////////////////////////////////////////////
    std::ofstream results_file(output_filename.c_str());
    results_file << std::fixed << std::setprecision(4);

    results_file << "Symulacja temperatury - Metoda Elementow Skonconych" << std::endl;
    results_file << "=================================================" << std::endl;
    results_file << "Czas symulacji: " << simulation_time << " s" << std::endl;
    results_file << "Krok czasowy: " << time_step << " s" << std::endl;
    results_file << "Temperatura poczatkowa: " << initial_temp << " K" << std::endl;
    results_file << "Liczba wezlow: " << n_nodes << std::endl;
    results_file << "Liczba elementow: " << grid.get_Elements_number() << std::endl;
    results_file << "=================================================" << std::endl << std::endl;

    
    results_file << "Czas t = 0 s" << std::endl;
    results_file << "-------------------" << std::endl;
    for (int i = 0; i < n_nodes; i++) results_file << "  Wezel " << i << ": T = " << t_previous[i] << " K" << std::endl;
    results_file << std::endl;

    ////////////////////////////////////////////////////////////////
    for (int t = time_step; t <= simulation_time; t += time_step)
    {
        results_file << "Czas t = " << t << " s" << std::endl;
        results_file << "-------------------" << std::endl;

        
        std::vector< std::vector<double> > A = globalH.getGlobalH();
        std::vector< std::vector<double> > C = globalH.getGlobalC();
        std::vector<double> b(n_nodes, 0.0);

        const std::vector<double>& P = globalH.getGlobalP();

        //Ciekawe rzeczy
		//Wzór: (H + C/delta_tau) * t_current = (C/delta_tau) * t_previous + P_global
        for (int i = 0; i < n_nodes; i++) 
        {
            for (int j = 0; j < n_nodes; j++) 
            {
                A[i][j] += C[i][j] / delta_tau;
            }
            
            b[i] = 0.0;

            for (int j = 0; j < n_nodes; j++) 
            {
                b[i] += C[i][j] / delta_tau * t_previous[j];
            }

            b[i] += P[i];
        }

        t_current = solveLU(A, b);

		double max_temp = t_current[0];
		double min_temp = t_current[0];

        for(int i =0; i<n_nodes; i++)
        {
            if(t_current[i] > max_temp)
                max_temp = t_current[i];
            if(t_current[i] < min_temp)
                min_temp = t_current[i];
		}

		//cout << "Czas t = " << t << " s: Max Temp = " << max_temp << " K, Min Temp = " << min_temp << " K" << endl;
		cout  << max_temp << " " << min_temp << endl;

        for (int i = 0; i < n_nodes; i++) 
        {
            results_file << "  Wezel " << i << ": T = " << t_current[i] << " K" << std::endl;
        }
        results_file << std::endl;

        
        t_previous = t_current;
    }

    
    results_file.close();

    std::cout << "\nSAVED: wyniki_temperatury.txt" << std::endl;

    return 0;
}