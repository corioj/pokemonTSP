//
//  coordinate.hpp
//  p4-poke
//
//  Created by John Corio on 11/24/19.
//  Copyright © 2019 John Corio. All rights reserved.
//

// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0

#ifndef coordinate_hpp
#define coordinate_hpp

#include <stdio.h>

enum class Terrain : char {
    Land,
    Water,
    Coast
};

struct Coordinate {
    // debate putting a distance here
    int x;
    int y;
    int prev = -1;
    bool disc = false;
    Terrain terr;
    
    Coordinate(int x_in, int y_in)
    : x(x_in), y(y_in)
    {
        
        // JUST
        // REWRITE
        // THIS
        // AT
        // SOME
        // POINT
        
        // x or y positive, water
        if (x > 0 || y > 0)
        {
            terr = Terrain::Land;
        }
        // y is 0, possible coast
        else if (y == 0)
        {
            // x <= 0, must be coast
            if (x > 0)
            {
                terr = Terrain::Land;
            }
            // x > 0, land
            else
            {
                terr = Terrain::Coast;
            }
        }
        // x is 0, possible coast
        else if (x == 0)
        {
            // y > 0, land
            if (y > 0)
            {
                terr = Terrain::Land;
            }
            // y <= 0, coast
            else
            {
                terr = Terrain::Coast;
            }
        }
        // both negative, water
        else
        {
            terr = Terrain::Water;
        }
    }
    
};

#endif /* coordinate_hpp */
