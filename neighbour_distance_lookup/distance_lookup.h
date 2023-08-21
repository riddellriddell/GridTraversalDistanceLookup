#pragma once
#include <array>
#include <cmath>

#include "shared_types.h"

#include "../ExternalLibraries/gcem-master/include/gcem.hpp"

// CPP program to find MSB
// number for a given POSITIVE n.
#ifdef __GNUC__
#define clz(x) __builtin_clz(x)
#define ctz(x) __builtin_ctz(x)
#else
#include <immintrin.h>
#define clz(x) _lzcnt_u32(x)
#define ctz(x) _lzcnt_u64(x) 
#endif

template<uint32_t TLookupSize>
static consteval std::array<float, TLookupSize* TLookupSize * 2> calculate_dist_to_neighbour_lookup()
{
    std::array<float, TLookupSize * TLookupSize * 2>  out_array = {};

    //calculate x axis distance changes
    for (int iy = 0; iy < TLookupSize; ++iy)
    {
        for (int ix = 0; ix < TLookupSize; ++ix)
        {
            float dist_to_current_cell = gcem::sqrt(static_cast<float>((ix * ix) + (iy * iy)));
            float dist_to_neigbour_cell = gcem::sqrt(static_cast<float>(((ix + 1) * (ix + 1)) + ((iy) * (iy))));
            float dist_dif = dist_to_neigbour_cell - dist_to_current_cell;

            out_array._Elems[(iy * TLookupSize) + ix] = dist_dif;
        }
    }

    //start of long axis lookup
    constexpr int y_axis_lookup_start = TLookupSize * TLookupSize;

    //calculate y axis distance changes
    for (int iy = 0; iy < TLookupSize; ++iy)
    {
        for (int ix = 0; ix < TLookupSize; ++ix)
        {
            float dist_to_current_cell = gcem::sqrt(static_cast<float>((ix * ix) + (iy * iy)));
            float dist_to_neigbour_cell = gcem::sqrt(static_cast<float>((ix * ix) + ((iy + 1) * (iy + 1))));
            float dist_dif = dist_to_neigbour_cell - dist_to_current_cell;

            out_array._Elems[(iy * TLookupSize) + ix + y_axis_lookup_start] = dist_dif;
        }
    }

    return out_array;
}


template<uint32_t TLookupSize>
static consteval std::array<float, TLookupSize* TLookupSize> calculate_convert_manhattan_lookup()
{
    std::array<float, TLookupSize* TLookupSize >  out_array = {};

    //calculate x axis distance changes
    for (int iy = 0; iy < TLookupSize; ++iy)
    {
        for (int ix = 0; ix < TLookupSize; ++ix)
        {
            float manhattan_dist =static_cast<float>(ix + iy);

            float dist_to_current_cell = gcem::sqrt(static_cast<float>((ix * ix) + (iy * iy)));
           
            if (manhattan_dist == 0)
            {
                out_array._Elems[(iy * TLookupSize) + ix] = 0;
            }
            else
            {
                out_array._Elems[(iy * TLookupSize) + ix] = dist_to_current_cell / manhattan_dist;
            }
        }
    }

    return out_array;
}

#pragma region compile_time_checks_for_calculate_dist_to_neighbour_lookup

//check that moving across one on the x is always 1
//due to float point problems added check with epsilon 
static_assert( gcem::abs(calculate_dist_to_neighbour_lookup<16>()[1] - 1) < std::numeric_limits<float>::epsilon());
static_assert( calculate_dist_to_neighbour_lookup<16>()[15] == ( gcem::sqrt(static_cast<float>(16 * 16)) - gcem::sqrt(static_cast<float>(15 * 15))));

//y move from (1,0) to (1,1) should give the distance of sqrt( 1 * 1 + 1 * 1 )
static_assert(calculate_dist_to_neighbour_lookup<16>()[(16 * 16) + 1] == gcem::sqrt(2.0f) - gcem::sqrt(1.0f));
static_assert(gcem::abs(calculate_dist_to_neighbour_lookup<16>()[(16 * 16) + 1] + calculate_dist_to_neighbour_lookup<16>()[1] - gcem::sqrt(2.0f)) < std::numeric_limits<float>::epsilon());

#pragma endregion

static consteval int calculate_leading_zeros(uint32_t value)
{
    int leading_zeros = 0;

    int bits_in_type = sizeof(uint32_t) * 8;

    for (int i = 0; i < sizeof(uint32_t) * 8; ++i)
    {
        if (value >> ((bits_in_type - 1) - i))
        {
            return i;
        }
    }

    return bits_in_type;
}

#pragma region compile_time_checks_for_calculate_left_most_bit
static_assert((sizeof(uint32_t) * 8) == 32);
static_assert(calculate_leading_zeros(0) == 32);
static_assert(calculate_leading_zeros(1) == 31);
static_assert(calculate_leading_zeros(std::numeric_limits<uint32_t>::max()) == 0);
static_assert(calculate_leading_zeros(std::numeric_limits<uint32_t>::max()>> 1) == 1);
#pragma endregion

template<uint32_t TLookupSize>
class distance_lookup
{
    static constexpr size_t total_lookup_array_size = TLookupSize * TLookupSize * 2;
    
    static constexpr std::array<float,total_lookup_array_size> dist_to_neighbour = calculate_dist_to_neighbour_lookup<TLookupSize>();

    static constexpr size_t convert_manhattan_to_crow_array_size = TLookupSize * TLookupSize;

    static constexpr std::array<float, convert_manhattan_to_crow_array_size> convert_manhattan_to_crow = calculate_convert_manhattan_lookup<TLookupSize>();
public:
    static float get_neighbour_distance(
        const tile_coordinate& from,
        const tile_coordinate& to,
        float current_dist,
        bool is_perpendicular_neighbour);

    static float approximate_vector_size(tile_coordinate_difference vector);


private:

};

template <uint32_t TLookupSize>
float distance_lookup<TLookupSize>::get_neighbour_distance(
    const tile_coordinate& from,
    const tile_coordinate& to,
    float current_dist,
    bool is_y_axis_neighbour)
{
    //get difference
    auto cord_diff = from - to;

    //remove sign values and convert biggest to 
    uint32_t x = static_cast<uint32_t>(abs(cord_diff.x));
    uint32_t y = static_cast<uint32_t>(abs(cord_diff.y));

    //get the bit shift amount
    int leading_zeros = clz( x | y);

    constexpr int lookup_leading_zeros = calculate_leading_zeros(TLookupSize - 1);

    int bit_shift_amount = std::max(lookup_leading_zeros - leading_zeros, 0);

    //shift the coordinate difference into the same space as the lookup table
    x = x >> bit_shift_amount;
    y = y >> bit_shift_amount;

    //this code is to pack both the distance to the short axis neighbour
    //and long axis neighbour into the same table
    //this is a trade off of extra instructions vs cache coherence
    int lookup_index = {};

    //start of long axis lookup
    constexpr int y_axis_lookup_start = TLookupSize * TLookupSize;
        
    lookup_index = (x + y * TLookupSize) + (y_axis_lookup_start * is_y_axis_neighbour);

    auto distance_difference_to_neighbour = dist_to_neighbour[lookup_index];

    return current_dist + distance_difference_to_neighbour;

}

template<uint32_t TLookupSize>
inline float distance_lookup<TLookupSize>::approximate_vector_size(tile_coordinate_difference vector)
{
    //get largest axis
     //remove sign values and convert biggest to
    uint32_t x = static_cast<uint32_t>(abs(vector.x));
    uint32_t y = static_cast<uint32_t>(abs(vector.y));
    
    //get the bit shift amount
    int leading_zeros = clz(x | y);

    constexpr int lookup_leading_zeros = calculate_leading_zeros(TLookupSize - 1);

    int bit_shift_amount = std::max(lookup_leading_zeros - leading_zeros, 0);

    //shift the coordinate difference into the same space as the lookup table
    uint32_t x_lookup = x >> bit_shift_amount;
    uint32_t y_lookup = y >> bit_shift_amount;

    //lookup the manhattan distance multiplier that converts from manhatan to true 
    float dist_multiplier = convert_manhattan_to_crow[(y_lookup * TLookupSize) + x_lookup];
    
    return (x + y) * dist_multiplier;

}