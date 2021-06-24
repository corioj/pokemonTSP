//
//  poke.cpp
//  p4-poke
//
//  Created by John Corio on 11/23/19.
//  Copyright © 2019 John Corio. All rights reserved.
//

// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0

#include <stdio.h>
#include <iostream>
#include <iomanip> // for std::precision and std::fixed
#include <algorithm> // for std::fill
#include "xcode_redirect.hpp"
#include <getopt.h>
#include "coordinate.hpp"
#include <queue>
#include <vector>
#include <limits> // for infinity
#include <cmath> // for sqrt
#include "string.h"

// or do i just do linear iteration thru the vector, pq is O(E*logV) --> E = V^2 since graph is dense

double euclidean(const Coordinate &c1, const Coordinate &c2)
{
    double diff1 = (c1.x - c2.x);
    double diff2 = (c1.y - c2.y);
    double result = sqrt(diff1 * diff1 + diff2 * diff2);
    return result;
}

class MST {
public:
    // MST constructor, initialize vector to # of vertices
    MST(size_t size_in)
    {
        coordinates.reserve(size_in);
        minspanDist = 0;
    }
    
    void insert_coordinate(Coordinate &coord)
    {
        coordinates.emplace_back(coord);
    }
    
    bool check_terrains(const Coordinate &c1, const Coordinate &c2)
    {
        // land = land, water = land, coast = coast
        if (c1.terr == c2.terr)
        {
            return true;
        }
        // either one is coast
        else if (c1.terr == Terrain::Coast || c2.terr == Terrain::Coast)
        {
            return true;
        }
        // neither are coast and they are not equal, must be land -> water
        return false;
        
    }
    
    // iterate over all coordinates, edit the edge in graph vector at current coordinate index
    // if undiscovered, i.e land->land, water->water, coast->land/water or vice versa, go into distance calcs and overwrite. keep track of minimum index (HAVE TO FIGURE OUT LOGISTICS OF THAT, i.e. COMPARE IT TO WHAT SMALLEST DISTANCE IN FIRST PLACE, INFINITE?)
    void prims_algorithm()
    {
        // our table of doubles for distances, could store double within Coordinate and make it easier
        std::vector<double> graph;
        graph.resize(coordinates.size());
        //std::fill(graph.begin(), graph.begin() + coordinates.size(), std::numeric_limits<double>::infinity());
        
        for (size_t i = 0; i < graph.size(); i++)
        {
            graph[i] = std::numeric_limits<double>::infinity();
        }
        
        // extra variables to make switching between points easier
        int index_of_min = 0;
        double min_dist = std::numeric_limits<double>::infinity();
        int previndex = 0;
        // mark first distance as 0 and first coord as discovered
        graph[0] = 0;
        coordinates[0].disc = true;

        // start at coordinate recently marked true, will be updated at end by prevIndex, HAVE TO UPDATE LOOP LATER
        // run (# of V) - 1 times, marked 0 as discovered already
        for (size_t i = 1; i < coordinates.size(); i++)
        {
            // loop over rest of coordinates
            for (size_t j = 0; j < coordinates.size(); j++)
            {
                // if undiscovered and not examining the same vertex, go further
                if (!coordinates[j].disc)
                {
                    // checking if terrain matches up
                    if (check_terrains(coordinates[previndex], coordinates[j]))
                    {
                        // calc euclidean distance, if smaller than smallest found so far distance to coordinates[j] then we overwrite graph[j].first with distance btw coordinates[i] & [j] and overwrite graph[j].second with previndex
                        double distance = euclidean(coordinates[previndex], coordinates[j]);
                        if (distance < graph[j])
                        {
                            graph[j] = distance;
                            coordinates[j].prev = previndex;
                        }
                    }
                    // if distance between j and any point previously discovered is smaller than min_dist recorded in the undiscovered nodes so far, overwrite BOTH min_dist and index_of_min
                    if (!coordinates[j].disc && graph[j] < min_dist)
                    {
                        min_dist = graph[j];
                        index_of_min = static_cast<int>(j);
                    }
                }
            }
            if (min_dist == std::numeric_limits<double>::infinity())
            {
                std::cerr << "Cannot construct MST";
                exit(1);
            }
            
            coordinates[index_of_min].disc = true;
            previndex = index_of_min;
            minspanDist += min_dist;
            // once this loop is done we have the next path to take and its distance, add to totaldistance
            min_dist = std::numeric_limits<double>::infinity();
        }
    }
    
    void promising_prims_algorithm()
    {
        // our table of doubles for distances, could store double within Coordinate and make it easier
        std::vector<double> graph;
        graph.resize(coordinates.size());
        //std::fill(graph.begin(), graph.begin() + coordinates.size(), std::numeric_limits<double>::infinity());
        
        for (size_t i = 0; i < graph.size(); i++)
        {
            graph[i] = std::numeric_limits<double>::infinity();
        }
        
        // extra variables to make switching between points easier
        int index_of_min = 0;
        double min_dist = std::numeric_limits<double>::infinity();
        int previndex = 0;
        // mark first distance as 0 and first coord as discovered
        graph[0] = 0;
        coordinates[0].disc = true;
        
        // start at coordinate recently marked true, will be updated at end by prevIndex, HAVE TO UPDATE LOOP LATER
        // run (# of V) - 1 times, marked 0 as discovered already
        for (size_t i = 1; i < coordinates.size(); i++)
        {
            // loop over rest of coordinates
            for (size_t j = 0; j < coordinates.size(); j++)
            {
                // if undiscovered and not examining the same vertex, go further
                if (!coordinates[j].disc)
                {
                    double distance = euclidean(coordinates[previndex], coordinates[j]);
                    if (distance < graph[j])
                    {
                        graph[j] = distance;
                        coordinates[j].prev = previndex;
                    }
                }
                // if distance between j and any point previously discovered is smaller than min_dist recorded in the undiscovered nodes so far, overwrite BOTH min_dist and index_of_min
                if (!coordinates[j].disc && graph[j] < min_dist)
                {
                    min_dist = graph[j];
                    index_of_min = static_cast<int>(j);
                }
            }
            coordinates[index_of_min].disc = true;
            previndex = index_of_min;
            minspanDist += min_dist;
            // once this loop is done we have the next path to take and its distance, add to totaldistance
            min_dist = std::numeric_limits<double>::infinity();
        }
    }
    
    double get_distance()
    {
        return minspanDist;
    }
    
    void print_spanning_tree()
    {
        // STUFF IS WRONG HERE
        std::cout << get_distance() << "\n";
        for (size_t i = 1; i < coordinates.size(); i++)
        {
            if (coordinates[i].prev < static_cast<int>(i))
            {
                std::cout << coordinates[i].prev << " " << i << "\n";
            }
            else
            {
                std::cout << i << " " << coordinates[i].prev << "\n";
            }
        }
    }
    
    double min_PartC_distance(Coordinate &tourCoord, Coordinate &endCoord)
    {
        // compare each coordinate in the part C MST to a coordinate passed in from the OPTTSP path
        // will always be either the starting coordinate or the last coordinate in the path so far
        double minBeg = std::numeric_limits<double>::infinity();
        double minEnd = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < coordinates.size(); i++)
        {
            // do our distance calcs
            double eucBeg = euclidean(tourCoord, coordinates[i]);
            if (eucBeg < minBeg)
            {
                minBeg = eucBeg;
            }
            double eucEnd = euclidean(endCoord, coordinates[i]);
            if (eucEnd < minEnd)
            {
                minEnd = eucEnd;
            }
        }
        // return our minimum distance to the starting point
        return minEnd + minBeg;
    }
    
    Coordinate& get_coordinate(size_t i)
    {
        return coordinates[i];
    }
    
private:
    // vector of all our coordinates
    std::vector<Coordinate> coordinates;
    // spanning distance is useful to have lying around if we need MST for Part C (or Part B)
    double minspanDist;
    // the example in P4 Piazza has each Coordinate containing a distance, think about this
};

class FASTTSP {
public:
    FASTTSP(size_t x)
    {
        coordinates.reserve(x);
        graph.reserve(x);
        mintourDist = 0;
    }
    
    FASTTSP(std::vector<Coordinate> &cool)
    {
        coordinates.reserve(cool.size());
        // THIS COULD CAUSE ISSUES SOON
        coordinates.insert(coordinates.begin(), cool.begin(), cool.end());
        graph.reserve(cool.size());
        mintourDist = 0;
    }
    
    void insert_coordinate(Coordinate &coord)
    {
        coordinates.emplace_back(coord);
    }
    
    // POINTS MAYBE, TIME TABLE NOT LOOKING TOO HOT THO THANKS 281 FOR RUINING LIFE FOR ENJOYING BREAK :))
    void arbitrary_insertion()
    {
        graph.push_back(0);
        coordinates[0].disc = true;
        
        // NEAREST NEIGHBOR METHOD
        
        double nearestneighbor = std::numeric_limits<double>::infinity();
        size_t indexofmin = 0;
        for (size_t i = 1; i < coordinates.size(); i++)
        {
            double x = partialEuclidean(coordinates[0], coordinates[i]);
            if (x < nearestneighbor)
            {
                nearestneighbor = x;
                indexofmin = i;
            }
        }
        // add nearest neighbor
        graph.push_back(indexofmin);
        coordinates[indexofmin].disc = true;
        mintourDist += finishEuclidean(nearestneighbor);
         
        
        /*
        // 0-1-2 METHOD
        graph.push_back(1);
        graph.push_back(2);
        mintourDist += finishEuclidean(partialEuclidean(coordinates[0], coordinates[1]));
        mintourDist += finishEuclidean(partialEuclidean(coordinates[1], coordinates[2]));
         */

        // iterate thru rest of our coordinates
        for (size_t i = 1/* THIS WAS 3 IN 0-1-2 METHOD*/; i < coordinates.size(); i++)
        {
            // if not discovered, mark
            if (!coordinates[i].disc)
            {
                coordinates[i].disc = true;
                size_t insertHere = 1;
                double currentMin = std::numeric_limits<double>::infinity();
                // find minimum insertion cost
                for (size_t insertSpot = 1; insertSpot < graph.size(); insertSpot++)
                {
                    // if EDGE(Vcoords[graph[x]], Vcoords[i]) + EDGE(Vcoords[graph[x]], Vcoords[i]) - EDGE(Vcoords[graph[x]], Vcoords[graph[x+1]]) is smallest insertion, insert it
                    double edge1 = finishEuclidean(partialEuclidean(coordinates[graph[insertSpot - 1]], coordinates[i]));
                    double edge2 = finishEuclidean(partialEuclidean(coordinates[graph[insertSpot]], coordinates[i]));
                    double edgeSub = finishEuclidean(partialEuclidean(coordinates[graph[insertSpot - 1]], coordinates[graph[insertSpot]]));
                    if (edge1 + edge2 - edgeSub < currentMin)
                    {
                        insertHere = insertSpot;
                        currentMin = (edge1 + edge2 - edgeSub);
                    }
                }
                // COMPARE TO LAST EDGE IN THE TOUR BETWEEN END POINT AND 0 MAYBE???
                
                // insert at position of minimum intrusion
                //if (insertHere != graph.size() - 1)
                //{
                mintourDist -= finishEuclidean(partialEuclidean(coordinates[graph[insertHere - 1]], coordinates[graph[insertHere]]));
                mintourDist += finishEuclidean(partialEuclidean(coordinates[graph[insertHere - 1]], coordinates[i]));
                mintourDist += finishEuclidean(partialEuclidean(coordinates[graph[insertHere]], coordinates[i]));
                //}
                /*else
                {
                    mintourDist += finishEuclidean(partialEuclidean(coordinates[graph[insertHere]], coordinates[i]));
                }*/
                graph.insert(graph.begin() + insertHere, i);
            }
        }
        // final edge
        mintourDist += finishEuclidean(partialEuclidean(coordinates[graph[0]], coordinates[graph[graph.size() - 1]]));
    }

    void print_TSP()
    {
        std::cout << mintourDist << "\n";
        for (size_t g = 0; g < graph.size(); g++)
        {
            std::cout << graph[g] << " ";
        }
    }
    
    double partialEuclidean(Coordinate &c1, Coordinate &c2)
    {
        double diff1 = (c1.x - c2.x);
        double diff2 = (c1.y - c2.y);
        return (diff1 * diff1) + (diff2 * diff2);
    }
    
    double finishEuclidean(double tobeSQRT)
    {
        return sqrt(tobeSQRT);
    }
    
    double get_mintourDist()
    {
        return mintourDist;
    }
    
    std::vector<size_t> copy_partC()
    {
        return graph;
    }
    
private:
    // vector of all our coordinates
    std::vector<Coordinate> coordinates;
    // double vector for other class functions
    std::vector<size_t> graph;
    // handy to have this lying around
    double mintourDist;
};

class OPTTSP {
public:
    OPTTSP(size_t x)
    {
        coordinates.reserve(x);
        solution.reserve(x);
        finalsol.reserve(x);
        distanceMatrix.resize(x, std::vector<double>(x, 0));
        upperbound = std::numeric_limits<double>::infinity();
        currdistance = 0;
        minspanDist = 0;
    }
    
    void insert_coordinate(Coordinate &x)
    {
        coordinates.emplace_back(x);
    }
    
    void filldistancematrix()
    {
        for (size_t row = 0; row < coordinates.size(); row++)
        {
            for (size_t col = 0; col < coordinates.size(); col++)
            {
                if (row != col && distanceMatrix[row][col] == 0 && distanceMatrix[col][row] == 0)
                {
                    distanceMatrix[row][col] = euclidean(coordinates[row], coordinates[col]);
                    distanceMatrix[col][row] = distanceMatrix[row][col];
                }
            }
        }
    }
    
    void genUpperbound()
    {
        // THIS MAY CAUSE ISSUES LATER ON, COULD FUCK UP OUR ORIGINAL COORDINATE ORDERING SINCE WE PASS IN OUR MEMBER VARIABLE VECTOR
        FASTTSP cool(coordinates);
        cool.arbitrary_insertion();
        upperbound = cool.get_mintourDist();
        solution = cool.copy_partC();
    }
    
    void genPerms(size_t permLength)
    {
        // SAMPLE E OPTTSP
        // 328.77
        // 0 9 7 8 4 3 10 5 2 1 6
        
        if (coordinates.size() == permLength)
        {
            // edit upperbound
            if (currdistance + distanceMatrix[0][solution[permLength - 1]] < upperbound)
            {
                upperbound = currdistance + distanceMatrix[0][solution[permLength - 1]];
                finalsol = solution;
            }
            return;
        }
        if (!promising(permLength))
        {
            return;
        }
        for (size_t i = permLength; i < solution.size(); ++i)
        {
            // swap coordinates within solution
            std::swap(solution[permLength], solution[i]);
            
            // EUCLIDEAN CALC UPDATE
            //currdistance += euclidean(coordinates[solution[permLength]], coordinates[solution[permLength-1]]);
            
            
            // DISTANCE MATRIX UPDATE
            currdistance += distanceMatrix[solution[permLength]][solution[permLength-1]];
            
            
            
            // go into new recursive layer
            genPerms(permLength + 1);
            
            // EUCLIDEAN CALC UPDATE
            //currdistance -= euclidean(coordinates[solution[permLength]], coordinates[solution[permLength-1]]);
            
            
            // DISTANCE MATRIX UPDATE
            currdistance -= distanceMatrix[solution[permLength]][solution[permLength-1]];
            
            
            
            // unswap coordinates within solution
            std::swap(solution[permLength], solution[i]);
        }
     
    }
     
    
    // MAKE ALL CALLS TO COORDINATE VECTOR BE REFERENCED THRU THE SOLUTION VECTOR
    bool promising(size_t permLength)
    {
        //std::cout << "entered promising\n";
        if (solution.size() - permLength < 5)
        {
            //std::cout << "entered promising size dif " << permLength << "\n";
            return true;
        }
        // NEED TO CLEAN UP MY TREE HERE, COULD MARK COORDINATES ADDED BY GENPERMS AS DISCOVERED AND THEN DO A MEMBER VARIABLE MST THAT ONLY PERFORMS A PRIMS ALGORITHM ON THOSE INSTEAD OF CREATING AND PUSHING AN MST FOR EVERY SINGLE PERMUTATION, THEN I MAYBE CAN USE DISTANCE MATRIX FOR THIS OR SOMETHING ELSE
        
        MST partCTree(coordinates.size() - permLength);
        for (size_t i = permLength; i < solution.size(); i++)
        {
            partCTree.insert_coordinate(coordinates[solution[i]]);
        }
        // run Part C Prim's
        partCTree.promising_prims_algorithm();
        // find closest edge weights to starting and ending point of permLength vector
        return currdistance + partCTree.get_distance() + partCTree.min_PartC_distance(coordinates[0], coordinates[solution[permLength - 1]]) < upperbound;
         
        //std::cout << "entered promising calcs " << permLength <<"\n";
        //return currdistance + promising_prims_algorithm_insideC(permLength) + findmindist_to_MST(permLength) < upperbound;
    }
    
    void print_OPTTSP()
    {
        // add distance between final edge and beginning
        std::cout << upperbound << "\n";
        for (size_t i = 0; i < finalsol.size(); i++)
        {
            std::cout << finalsol[i] << " ";
        }
    }
    
    // could pass in the permlength variable and then iterate over the solution vector using it instead of iterating over coordinates, use distance matrix to ease all this, make separate functions for the min partC distance, JUST HAVE TO BE CAREFUL TO RESET THE DISCOVERY AND PRIM VALUES EACH TIME
    double promising_prims_algorithm_insideC(size_t permLength)
    {
        // reset the distance
        minspanDist = 0;
        //std::cout << "entered promising prims\n";
        // extra variables to make switching between points easier
        int index_of_min = 0;
        double min_dist = std::numeric_limits<double>::infinity();
        int previndex = 0;
        // mark first coord in solution as discovered
        coordinates[solution[permLength]].disc = true;
        
        // start at coordinate recently marked true, will be updated at end by prevIndex
        // run (# of V) - 1 times, marked solution[permlength] as discovered already
        for (size_t i = permLength + 1; i < solution.size(); i++)
        {
            // loop over rest of coordinates
            for (size_t j = permLength; j < solution.size(); j++)
            {
                // if distance between j and any point previously discovered is smaller than min_dist recorded in the undiscovered nodes so far, overwrite BOTH min_dist and index_of_min
                if (!coordinates[j].disc && distanceMatrix[solution[previndex]][solution[j]] < min_dist)
                {
                    min_dist = distanceMatrix[solution[previndex]][solution[j]];
                    index_of_min = static_cast<int>(j);
                }
            }
            coordinates[solution[index_of_min]].disc = true;
            previndex = index_of_min;
            minspanDist += min_dist;
            // once this loop is done we have the next path to take and its distance, add to totaldistance
            min_dist = std::numeric_limits<double>::infinity();
        }
        undo_prims(permLength);
        return minspanDist;
    }
    
    double findmindist_to_MST(size_t permLength)
    {
        //std::cout << "entered mindistMST\n";
        double minBeg = std::numeric_limits<double>::infinity();
        double minEnd = std::numeric_limits<double>::infinity();
        for (size_t i = permLength; i < solution.size(); i++)
        {
            // do our distance calcs
            if (distanceMatrix[0][i] < minBeg)
            {
                minBeg = distanceMatrix[0][solution[i]];
            }
            if (distanceMatrix[permLength-1][i] < minEnd)
            {
                minEnd = distanceMatrix[permLength-1][solution[i]];
            }
        }
        // return our minimum distance to the starting point
        return minBeg + minEnd;
    }
    
    void undo_prims(size_t permlength)
    {
        //std::cout << "entered undo prims\n";
        for (size_t idx = permlength; idx < solution.size(); idx++)
        {
            coordinates[solution[idx]].disc = false;
        }
    }
    
private:
    // distance matrix
    std::vector<std::vector<double>> distanceMatrix;
    // coordinates
    std::vector<Coordinate> coordinates;
    // running solution vector
    std::vector<size_t> solution;
    // final solution vector
    std::vector<size_t> finalsol;
    // running distance count
    double currdistance;
    // upper bound
    double upperbound;
    // MST dist
    double minspanDist;
};

// PRESSING QUESTIONS:

// How to create an MST of just undiscovered variables?

// HOW TO JUGGLE THE VECTOR OF SIZE_T FOR THE PATHS AND COORDINATE VECTOR??? TAKE SOME HINTS FROM PART B

// GO FIX ALL OF THE CALLS TO COORDINATES IN THIS, DO A SIMILAR THING TO PART B WHERE YOU USE SOLUTION AS INDICES INTO THE COORDINATES VECTOR

std::string optarg_manip(char* optarg)
{
    std::string mode(optarg);
    return mode;
}


int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    xcode_redirect(argc, argv);
    std::cout << std::setprecision(2);
    std::cout << std::fixed;
    
    // get opt handling
    int gotopt;
    int opt_index = 0;
    option long_opts[] = {
        {"help", no_argument, nullptr, 'h'},
        {"mode", required_argument, nullptr, 'm'},
        { nullptr, 0, nullptr, '\0'}
    }; // option, REVISE LATER

    // yeah
    while((gotopt = getopt_long(argc, argv, "m:h", long_opts, &opt_index)) != -1) {
        switch (gotopt) {
        // HELP OPTION
        case 'h': {
            std::cerr << "hehehehe im not giving u any help whatsoever!\n";
            exit(0);
        }
            
        case 'm': {
            // do some stuff with either MST, FASTTSP, or OPTTSP
            //std::string mode = optarg_manip(optarg);
            // check whether it's MST, FASTTSP, or OPTTSP, OPTIMIZE LATER
            
            if (optarg[0] == '\0')
            {
                 // do something
            }
            else
            {
                std::string mode = optarg_manip(optarg);
                if (mode == "MST")
                {
                    size_t numCoords = 0;
                    std::cin >> numCoords;
                    MST minimumSpanningTree(numCoords);
                    int x, y;
                    while (std::cin >> x >> y)
                    {
                        Coordinate piano = Coordinate(x, y);
                        minimumSpanningTree.insert_coordinate(piano);
                    }
                    minimumSpanningTree.prims_algorithm();
                    minimumSpanningTree.print_spanning_tree();
                    
                }
                else if (mode == "FASTTSP")
                {
                    size_t numCoords = 0;
                    std::cin >> numCoords;
                    FASTTSP fastTSPsol(numCoords);
                    int x, y;
                    while (std::cin >> x >> y)
                    {
                        Coordinate piano = Coordinate(x, y);
                        fastTSPsol.insert_coordinate(piano);
                    }
                    //fastTSPsol.nearestNeighbor();
                    fastTSPsol.arbitrary_insertion();
                    fastTSPsol.print_TSP();
                }
                else if (mode == "OPTTSP")
                {
                    size_t numCoords = 0;
                    std::cin >> numCoords;
                    OPTTSP optTSPsol(numCoords);
                    int x, y;
                    while (std::cin >> x >> y)
                    {
                        Coordinate piano = Coordinate(x, y);
                        optTSPsol.insert_coordinate(piano);
                    }
                    // no integrated part A
                    optTSPsol.genUpperbound();
                    optTSPsol.filldistancematrix();
                    optTSPsol.genPerms(1);
                    optTSPsol.print_OPTTSP();
                }
                else
                {
                    std::cerr << "ur string is wrong\n";
                    exit(1);
                }
            }
            break;
        }
                
                
        default: {
            std::cerr << "UHHHHH NICE ERROR LOSER LOL\n";
            exit(1);
            break;
        }
                
        }
    }
}
