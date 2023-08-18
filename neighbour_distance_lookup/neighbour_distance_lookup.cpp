
#include <limits>
#include <queue>
#include <tuple>
#include <compare>

#include "data_grid.h"

//template<int Tsize,>
using data_grid_type = data_grid<1024,float> ;

int main(int argc, char* argv[])
{
    //create gid
    data_grid_type walked_grid = data_grid_type(std::numeric_limits<float>::max());
    data_grid_type ground_truth_grid = data_grid_type(std::numeric_limits<float>::max());

    tile_coordinate center_tile = tile_coordinate{ 512,512 };

    //fill the ground trooth grid with the true distance to all the tiles
    for (short iy = 0; iy < data_grid_type::width; ++iy)
    {
        for (short ix = 0; ix < data_grid_type::width; ++ix)
        {
            float dist_to_tile = ((ix - center_tile.x) * (ix - center_tile.x)) + ((iy - center_tile.y) * (iy - center_tile.y));

            ground_truth_grid[tile_coordinate{ ix,iy}] = dist_to_tile;
        }
    }

    //set center tile to o
    walked_grid[center_tile] = 0.0f;

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

    
    std::priority_queue<open_tile_entry> open_queue;

    //add the first tile to the open list
    open_queue.emplace(open_tile_entry{ 0,center_tile });
    
    //start flood fill in all directions 
    while (open_queue.size())
    {
        //get the smalles value on the open queue
        auto active_cell = open_queue.top();

        //get active cell cord
        auto& active_cell_cord = std::get<1>(active_cell);

        //get manhatten dist to cell
        auto active_manhattan_distance = abs(active_cell_cord.x) + abs(active_cell_cord.y);

        //loop through neighbouring cells 
        for (int i = 0; i < static_cast<int>(NEIGHBOUR_DIRECTIONS::CARDINAL_COUNT); ++i)
        {
           auto neighbour_cord = tile_coordinate{
                static_cast<short>(std::get<1>(active_cell).x + neighbour_address_offset<short, data_grid_type::width>::x_offset_for_direction[i]),
                static_cast<short>(std::get<1>(active_cell).y + neighbour_address_offset<short, data_grid_type::width>::x_offset_for_direction[i])};

           //skip if out of bounds
           if (data_grid_type::is_valid(neighbour_cord))
           {
               continue;
           }

           auto neighbour_manhattan_distance = abs(neighbour_cord.x) + abs(neighbour_cord.y);

           //skip if it is close to origin or has already been writen to
           if (walked_grid[neighbour_cord] != std::numeric_limits<float>::max() && neighbour_manhattan_distance > active_manhattan_distance)
           {
              auto difference = 
           }
        }
    }
    return 0;
}
