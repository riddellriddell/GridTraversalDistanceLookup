#pragma once

//not able to store tile coordinate differences
struct tile_coordinate
{
    short x;
    short y;

    struct tile_coordinate_difference operator - (const tile_coordinate& cord_to_subtract) const;
    
};

struct tile_coordinate_difference
{
    int x;
    int y;
};

enum class NEIGHBOUR_DIRECTIONS
{
    LEFT,
    UP,
    RIGHT,
    DOWN,
    UP_LEFT,
    UP_RIGHT,
    DOWN_RIGHT,
    DOWN_LEFT,
    COUNT,

    CARDINAL_COUNT = UP_LEFT,
    DIAGONAL_COUNT = COUNT - UP_LEFT
};