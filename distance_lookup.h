#pragma once
#include <array>
#include <cmath>

#include "shared_types.h"

// CPP program to find MSB
// number for a given POSITIVE n.
#ifdef __GNUC__
#define clz(x) __builtin_clz(x)
#define ctz(x) __builtin_ctz(x)
#else
#include <immintrin.h>
#define clz(x) _tzcnt_u32(x)
#define ctz(x) _tzcnt_u64(x) 
#endif


template<uint32_t TLookupSize>
class distance_lookup
{
    static constexpr size_t total_lookup_array_size = TLookupSize * TLookupSize * 2;
    
    static std::array<float,total_lookup_array_size> dist_to_neighbour;
public:
    static float get_neighbour_distance(
        const tile_coordinate& from,
        const tile_coordinate& to,
        float current_dist,
        bool is_perpendicular_neighbour);

private:

};

template <uint32_t TLookupSize>
float distance_lookup<TLookupSize>::get_neighbour_distance(
    const tile_coordinate& from,
    const tile_coordinate& to,
    float current_dist,
    bool is_short_axis_neighbour)
{
    //get difference
    auto cord_diff = from - to;

    //remove sign values and convert biggest to 
    uint32_t x = static_cast<uint32_t>(abs(cord_diff.x));
    uint32_t y = static_cast<uint32_t>(abs(cord_diff.x));

    //get the bit shift amount
    int left_most_bit = clz( x | y);

    constexpr int lookup_max_val = clz(TLookupSize);

    int bit_shift_amount = std::max(left_most_bit - lookup_max_val, 0);

    //shift the coordinate difference into the same space as the lookup table
    x = x >> bit_shift_amount;
    y = y >> bit_shift_amount;

    //this code is to pack both the distance to the short axis neighbour
    //and long axis neighbour into the same table
    //this is a trade off of extra instructions vs cache coherence
    int lookup_index = {};

    //start of long axis lookup
    constexpr int short_axis_lookup_start = TLookupSize * TLookupSize;
        
    lookup_index = (x + y * TLookupSize) + (short_axis_lookup_start * is_short_axis_neighbour);


    auto distance_difference_to_neighbour = dist_to_neighbour[lookup_index];

    return current_dist + distance_difference_to_neighbour;

}

