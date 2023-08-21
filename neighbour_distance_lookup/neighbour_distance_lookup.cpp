
#include <limits>
#include <queue>
#include <tuple>
#include <compare>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "data_grid.h"
#include "distance_lookup.h"

constexpr int test_size = 1024;
constexpr int lookup_size = 16;
//template<int Tsize,>
using data_grid_type = data_grid<test_size,float> ;

int main(int argc, char* argv[])
{
    //create gid
    auto ground_truth_grid = std::make_unique<data_grid_type>(std::numeric_limits<float>::max());
    auto walked_grid = std::make_unique<data_grid_type>(std::numeric_limits<float>::max());
    auto estimated_from_manhattan_grid = std::make_unique<data_grid_type>(std::numeric_limits<float>::max());

    auto relative_error_grid = std::make_unique<data_grid_type>(std::numeric_limits<float>::max());

    tile_coordinate center_tile = tile_coordinate{ test_size / 2, test_size / 2 };
    
    //fill the ground trooth grid with the true distance to all the tiles
    for (short iy = 0; iy < data_grid_type::width; ++iy)
    {
        for (short ix = 0; ix < data_grid_type::width; ++ix)
        {
            float dist_to_tile =  sqrtf(static_cast<float>(((ix - center_tile.x) * (ix - center_tile.x)) + ((iy - center_tile.y) * (iy - center_tile.y))));
    
            (*ground_truth_grid)[tile_coordinate{ ix,iy}] = dist_to_tile;
        }
    }

    //fill the manhattan grid with the estimated distances
    for (short iy = 0; iy < data_grid_type::width; ++iy)
    {
        for (short ix = 0; ix < data_grid_type::width; ++ix)
        {
            (*estimated_from_manhattan_grid)[tile_coordinate{ ix,iy }] = distance_lookup<lookup_size>::approximate_vector_size(center_tile - tile_coordinate{ ix,iy });
        }
    }
    
    //set center tile to o
    (*walked_grid)[center_tile] = 0.0f;
    
    //struct to hold items in open queue
    struct open_tile_entry: public std::tuple<float,tile_coordinate>
    {
        //i want to set the tuple to these values passed into the constructor
        constexpr open_tile_entry(float priority,tile_coordinate coorndinate):std::tuple<float, tile_coordinate>(priority,coorndinate){} 
    
        auto operator<=>(const open_tile_entry& other) const
        {
           return std::get<0>(*this) <=> std::get<0>(other);
        }
        
    };
    
    if (false)
    {

        //create open queue
        std::priority_queue<open_tile_entry> open_queue = std::priority_queue<open_tile_entry>();

        //add the first tile to the open list
        open_queue.emplace(open_tile_entry{ 0,center_tile });

        //start flood fill in all directions 
        while (open_queue.size())
        {
            //get the smalles value on the open queue
            open_tile_entry active_cell = open_queue.top();

            //remove smalles value from open queue
            open_queue.pop();

            //get active cell cord
            auto& active_cell_cord = std::get<1>(active_cell);

            //get distance to active cell
            auto active_cell_distance = std::get<0>(active_cell);

            //get manhatten dist to cell
            auto active_manhattan_distance = abs(active_cell_cord.x - center_tile.x) + abs(active_cell_cord.y - center_tile.y);

            //check if the move is in the y axis 
            constexpr bool is_y_move[static_cast<int>(NEIGHBOUR_DIRECTIONS::CARDINAL_COUNT)] = { false,true,false,true };

            //loop through neighbouring cells 
            for (int i = 0; i < static_cast<int>(NEIGHBOUR_DIRECTIONS::CARDINAL_COUNT); ++i)
            {
                auto neighbour_cord = tile_coordinate{
                     static_cast<short>(std::get<1>(active_cell).x + neighbour_address_offset<short, data_grid_type::width>::x_offset_for_direction[i]),
                     static_cast<short>(std::get<1>(active_cell).y + neighbour_address_offset<short, data_grid_type::width>::y_offset_for_direction[i]) };

                //skip if out of bounds
                if (!data_grid_type::is_valid(neighbour_cord))
                {
                    continue;
                }

                auto neighbour_manhattan_distance = abs(neighbour_cord.x - center_tile.x) + abs(neighbour_cord.y - center_tile.y);

                //skip if it is close to origin or has already been writen to
                if (neighbour_manhattan_distance > active_manhattan_distance)
                {
                    auto distance_to_neighbour = distance_lookup<lookup_size>::get_neighbour_distance(center_tile, active_cell_cord, active_cell_distance, is_y_move[i]);

                    if ((*walked_grid)[neighbour_cord] > distance_to_neighbour)
                    {

                        //write out to target cell
                        (*walked_grid)[neighbour_cord] = distance_to_neighbour;

                        //add neighbour to open list 
                        open_queue.push(open_tile_entry(distance_to_neighbour, neighbour_cord));
                    }
                }
            }
        }
    }

    float max_error = 0;
    tile_coordinate coordiate_of_max_error = {};

    float max_relative_error = 0;
    tile_coordinate coordiate_of_max_relative_error = {};

    data_grid_type* ground_truth_to_compare = ground_truth_grid.get();
    data_grid_type* estimation_to_compare = estimated_from_manhattan_grid.get();

    //compare the 2 arrays to calculate the error
    for (short iy = 0; iy < data_grid_type::width; ++iy)
    {
        for (short ix = 0; ix < data_grid_type::width; ++ix)
        {
            tile_coordinate cord_to_compare = tile_coordinate{ ix,iy };

            float error_at_cord = (*estimation_to_compare)[cord_to_compare] - (*ground_truth_to_compare)[cord_to_compare];
            
            float max_tile_neighbour_error = 0;

            for (int idir = 0; idir < static_cast<int>(NEIGHBOUR_DIRECTIONS::CARDINAL_COUNT); ++idir)
            {
                short neighbour_x = ix + neighbour_address_offset<short, data_grid_type::width>::x_offset_for_direction[idir];
                short neighbour_y = iy + neighbour_address_offset<short, data_grid_type::width>::y_offset_for_direction[idir];
            
                tile_coordinate neighbour_cord_to_compare = tile_coordinate{ neighbour_x , neighbour_y };
            
                if (!data_grid_type::is_valid(neighbour_cord_to_compare))
                {
                    continue;
                }

                float error_at_neighbour = (*estimation_to_compare)[neighbour_cord_to_compare] - (*ground_truth_to_compare)[neighbour_cord_to_compare];
            
                float relative_error = error_at_cord - error_at_neighbour;
            
                max_tile_neighbour_error = fmax(abs(relative_error), max_tile_neighbour_error);

            }

            (*relative_error_grid)[cord_to_compare] = max_tile_neighbour_error;

            if (abs(error_at_cord) > max_error)
            {
                max_error = abs(error_at_cord);
                coordiate_of_max_error = cord_to_compare;
            }  

            if (abs(max_tile_neighbour_error) > max_relative_error)
            {
                max_relative_error = abs(max_tile_neighbour_error);
                coordiate_of_max_relative_error = cord_to_compare;
            }
        }
    }

    std::cout << std::endl;
    std::cout << "printing out error results" << std::endl;

    std::stringstream string_out = {};

    for (short iy = 0; iy < data_grid_type::width; ++iy)
    {
        for (short ix = 0; ix < data_grid_type::width; ++ix)
        {
            tile_coordinate cord_to_compare = tile_coordinate{ ix,iy };

            float error_at_cord = (*estimation_to_compare)[cord_to_compare] - (*ground_truth_to_compare)[cord_to_compare];

            int rounded_error = static_cast<int>(abs(error_at_cord * 100.0f));

            int hundreds = rounded_error / 100;
            int tens = (rounded_error - hundreds) / 10;
            int singles = rounded_error % 10;

            string_out << std::setprecision(2) << hundreds << "." << tens << singles << ",";
        }

        string_out << std::endl;

    }

    std::stringstream neighbour_error_string_out = {};

    for (short iy = 0; iy < data_grid_type::width; ++iy)
    {
        for (short ix = 0; ix < data_grid_type::width; ++ix)
        {
            tile_coordinate cord_to_compare = tile_coordinate{ ix,iy };

            float error_at_cord = (*relative_error_grid)[cord_to_compare] ;

            int rounded_error = static_cast<int>(abs(error_at_cord * 100.0f));

            int hundreds = rounded_error / 100;
            int tens = (rounded_error - hundreds) / 10;
            int singles = rounded_error % 10;

            neighbour_error_string_out << std::setprecision(2) << hundreds << "." << tens << singles << ",";
        }

        neighbour_error_string_out << std::endl;

    }

    std::cout << string_out.str();

    std::cout << std::endl;
    std::cout << "printing out relative error results" << std::endl;
    std::cout << neighbour_error_string_out.str();
    

    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "results from neighbour dist estimation" << std::endl;
    std::cout << "biggest error was : " << max_error << std::endl;
    std::cout << "biggest error was in tile: (" << coordiate_of_max_error.x << "," << coordiate_of_max_error.y << ")" << std::endl;

    std::cout << "biggest relative error was : " << max_relative_error << std::endl;
    std::cout << "biggest relative error was in tile: (" << coordiate_of_max_relative_error.x << "," << coordiate_of_max_relative_error.y << ")" << std::endl;

    return 0;
}
