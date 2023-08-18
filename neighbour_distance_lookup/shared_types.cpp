#include "shared_types.h"

tile_coordinate_difference tile_coordinate::operator-(const tile_coordinate& cord_to_subtract) const
{
    tile_coordinate_difference out = {
        static_cast<int>(x) - static_cast<int>(cord_to_subtract.x),
        static_cast<int>(y) - static_cast<int>(cord_to_subtract.y)};

    return  out;
}
