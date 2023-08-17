
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
    data_grid_type grid = data_grid_type(std::numeric_limits<float>::max());

    //set center tile to o
    grid[tile_coordinate{512,512}] = 0.0f;

    //struct to hold items in open queue
    struct open_tile_entry: std::tuple<float,tile_coordinate>
    {
        constexpr open_tile_entry(float priority,tile_coordinate coorndinate)
        {
            
            priority,coorndinate
        } { }
        auto operator<=>(const open_tile_entry&) const = default;
        
    };
    
    std::priority_queue<open_tile_entry> open_queue;
    
    //start flood fill in all directions 
    
    return 0;
}
