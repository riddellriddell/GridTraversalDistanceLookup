#pragma once
#include <array>
#include <cmath>
#include <tuple>

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

#pragma region compile_time_checks_for_calculate_dist_to_neighbour_lookup

//check that moving across one on the x is always 1
//due to float point problems added check with epsilon 
static_assert(gcem::abs(calculate_dist_to_neighbour_lookup<16>()[1] - 1) < std::numeric_limits<float>::epsilon());
static_assert(calculate_dist_to_neighbour_lookup<16>()[15] == (gcem::sqrt(static_cast<float>(16 * 16)) - gcem::sqrt(static_cast<float>(15 * 15))));

//y move from (1,0) to (1,1) should give the distance of sqrt( 1 * 1 + 1 * 1 )
static_assert(calculate_dist_to_neighbour_lookup<16>()[(16 * 16) + 1] == gcem::sqrt(2.0f) - gcem::sqrt(1.0f));
static_assert(gcem::abs(calculate_dist_to_neighbour_lookup<16>()[(16 * 16) + 1] + calculate_dist_to_neighbour_lookup<16>()[1] - gcem::sqrt(2.0f)) < std::numeric_limits<float>::epsilon());

#pragma endregion


template<uint32_t TLookupSize>
static consteval std::array<float, TLookupSize* TLookupSize> calculate_convert_manhattan_lookup()
{
    std::array<float, TLookupSize* TLookupSize >  out_array = {};

    //calculate x axis distance changes
    for (int iy = 0; iy < TLookupSize; ++iy)
    {
        for (int ix = 0; ix < TLookupSize; ++ix)
        {
            bool use_average_center = (iy > (TLookupSize / 2)) || (ix > (TLookupSize / 2));

            float x_center = static_cast<float>(ix) + (0.5f * use_average_center);
            float y_center = static_cast<float>(iy) + (0.5f * use_average_center);


            float manhattan_dist =static_cast<float>(x_center + y_center);

            float dist_to_current_cell = gcem::sqrt(static_cast<float>((x_center * x_center) + (y_center * y_center)));
           
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


template<uint32_t TLookupSize>
static consteval  std::array<std::tuple<float, float>, TLookupSize* TLookupSize >  calculate_convert_x_y_lookup()
{
    std::array<std::tuple<float,float>, TLookupSize* TLookupSize >  out_array = {};

    //calculate x axis distance changes
    for (int iy = 0; iy < TLookupSize; ++iy)
    {
        for (int ix = 0; ix < TLookupSize; ++ix)
        {
            if ((ix + iy) == 0)
            {
                out_array._Elems[(iy * TLookupSize) + ix] = { 1.0f,1.0f };

                continue;
            }

            bool use_average_center = (iy > (TLookupSize / 2)) || (ix > (TLookupSize / 2));

            float x_center = static_cast<float>(ix) + (0.5f * use_average_center);
            float y_center = static_cast<float>(iy) + (0.5f * use_average_center);

            float dist_sqr = (x_center * x_center) + (y_center * y_center);

            float dist = gcem::sqrt(dist_sqr);

            float x_percent_effect_on_dist = (x_center * x_center) / ((x_center * x_center) + (y_center * y_center));
            float y_percent_effect_on_dist = (y_center * y_center) / ((x_center * x_center) + (y_center * y_center));

            float dist_caused_by_x = dist * x_percent_effect_on_dist;
            float dist_caused_by_y = dist * y_percent_effect_on_dist;

            if (x_center <= std::numeric_limits<float>::epsilon())
            {
                x_center = 1.0f;
            }

            if (y_center <= std::numeric_limits<float>::epsilon())
            {
                y_center = 1.0f;
            }

            float x_mul = dist_caused_by_x / x_center;
            float y_mul = dist_caused_by_y / y_center;
            
            out_array._Elems[(iy * TLookupSize) + ix] = std::tuple<float,float>{ x_mul ,y_mul };

        }
    }

    return out_array;
}

template<uint32_t TLookupSize>
std::array<std::tuple<float, float>, TLookupSize* TLookupSize >  calculate_convert_x_y_lookup_dynamic()
{
    std::array<std::tuple<float, float>, TLookupSize* TLookupSize >  out_array = {};

    //calculate x axis distance changes
    for (int iy = 0; iy < TLookupSize; ++iy)
    {
        for (int ix = 0; ix < TLookupSize; ++ix)
        {
            if ((ix + iy) == 0)
            {
                out_array._Elems[(iy * TLookupSize) + ix] = { 1.0f,1.0f };

                continue;
            }

            //bool use_average_center = (iy > (TLookupSize / 2)) || (ix > (TLookupSize / 2));
           // bool use_average_center = false;


            float x_center = static_cast<float>(ix);//+ (0.5f * use_average_center);
            float y_center = static_cast<float>(iy);//+ (0.5f * use_average_center);

            float dist_sqr = (x_center * x_center) + (y_center * y_center);

            float dist = gcem::sqrt(dist_sqr);

            float x_percent_effect_on_dist = (x_center * x_center) / ((x_center * x_center) + (y_center * y_center));
            float y_percent_effect_on_dist = (y_center * y_center) / ((x_center * x_center) + (y_center * y_center));

            float dist_caused_by_x = dist * x_percent_effect_on_dist;
            float dist_caused_by_y = dist * y_percent_effect_on_dist;

            if (x_center <= std::numeric_limits<float>::epsilon())
            {
                x_center = 1.0f;
            }

            if (y_center <= std::numeric_limits<float>::epsilon())
            {
                y_center = 1.0f;
            }

            float x_mul = dist_caused_by_x / x_center;
            float y_mul = dist_caused_by_y / y_center;

            out_array._Elems[(iy * TLookupSize) + ix] = std::tuple<float, float>{ x_mul ,y_mul };

        }
    }

    return out_array;
}

template<uint32_t TLookupSize>
static consteval std::array<float, TLookupSize* TLookupSize * 2> calculate_neighbour_manhattan_difference_lookup()
{
    std::array<float, TLookupSize* TLookupSize * 2>  out_array = {};

    //calculate x axis distance changes
    for (int iy = 0; iy < TLookupSize; ++iy)
    {
        for (int ix = 0; ix < TLookupSize; ++ix)
        {
            float manhattan_mul_to_current_cell = 1;
            if ((gcem::abs(ix) + gcem::abs(iy)) != 0)
            {
                float dist_to_current = static_cast<float>(gcem::sqrt(static_cast<float>((ix * ix) + (iy * iy))));
                float man_to_current = static_cast<float>(gcem::abs(ix) + gcem::abs(iy));
                
                manhattan_mul_to_current_cell = dist_to_current / man_to_current;
            }
           
            float dist_to_neighbour = gcem::sqrt(static_cast<float>(((ix + 1) * (ix + 1)) + ((iy) * (iy))));
            float man_to_neighbour = static_cast<float>(gcem::abs(ix + 1) + gcem::abs(iy));
            float  manhattan_mul_to_neigbour_cell = dist_to_neighbour / man_to_neighbour;
           
            float mul_dif = manhattan_mul_to_neigbour_cell - manhattan_mul_to_current_cell;

            out_array._Elems[(iy * TLookupSize) + ix] = mul_dif;
        }
    }

    //start of long axis lookup
    constexpr int y_axis_lookup_start = TLookupSize * TLookupSize;

    //calculate y axis distance changes
    for (int iy = 0; iy < TLookupSize; ++iy)
    {
        for (int ix = 0; ix < TLookupSize; ++ix)
        {
            float manhattan_mul_to_current_cell = 1;
            if ((gcem::abs(ix) + gcem::abs(iy)) != 0)
            {
                float dist_to_current = static_cast<float>(
                    gcem::sqrt(
                        static_cast<float>(
                            (ix * ix) + (iy * iy)
                            )
                    )
                    );
                float man_to_current = static_cast<float>(gcem::abs(ix) + gcem::abs(iy));

                manhattan_mul_to_current_cell = dist_to_current / man_to_current;
            }

            float dist_to_neighbour = gcem::sqrt(static_cast<float>(((ix) * (ix)) + ((iy + 1) * (iy + 1))));
            float man_to_neighbour = static_cast<float>(gcem::abs(ix) + gcem::abs(iy + 1));
            float  manhattan_mul_to_neigbour_cell = dist_to_neighbour / man_to_neighbour;

            float mul_dif = manhattan_mul_to_neigbour_cell - manhattan_mul_to_current_cell;

            out_array._Elems[(iy * TLookupSize) + ix + y_axis_lookup_start] = mul_dif;
        }
    }

    return out_array;
}

#pragma region compile_time_checks_for_calculate_neighbour_manhattan_difference_lookup

//check that moving across one on the x is always 1
//due to float point problems added check with epsilon 
static_assert( gcem::abs(calculate_neighbour_manhattan_difference_lookup<16>()[1]) < std::numeric_limits<float>::epsilon());


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
    
    static constexpr std::array<float, total_lookup_array_size> dist_to_neighbour = {};//calculate_dist_to_neighbour_lookup<TLookupSize>();

    static constexpr size_t convert_manhattan_to_dist_array_size = TLookupSize * TLookupSize;

    static constexpr std::array<float, convert_manhattan_to_dist_array_size> convert_manhattan_to_dist = calculate_convert_manhattan_lookup<TLookupSize>();

    static constexpr std::array<float, total_lookup_array_size> mul_change_to_neighbour = {};//calculate_neighbour_manhattan_difference_lookup<TLookupSize>();
public:
    static constexpr std::array<std::tuple<float, float>, convert_manhattan_to_dist_array_size > convert_x_y_to_dist =  calculate_convert_x_y_lookup<TLookupSize>();

    static float get_neighbour_distance(
        const tile_coordinate& from,
        const tile_coordinate& to,
        float current_dist,
        bool is_perpendicular_neighbour);

    static float approximate_vector_size(tile_coordinate_difference vector);

    static float approximate_vector_size_v2(tile_coordinate_difference vector);

    static float calculate_new_manhattan_multiplier(
        const tile_coordinate& from,
        const tile_coordinate& to,
        float current_manhattan_multiplier,
        bool is_perpendicular_neighbour);
  

private:

};

static_assert(std::get<0>(distance_lookup<16>::convert_x_y_to_dist[0]) == 1.0f);
static_assert(std::get<0>(distance_lookup<16>::convert_x_y_to_dist[1]) == 1.0f);


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
    float dist_multiplier = convert_manhattan_to_dist[(y_lookup * TLookupSize) + x_lookup];
    
    return (x + y) * dist_multiplier;

}

template<uint32_t TLookupSize>
inline float distance_lookup<TLookupSize>::approximate_vector_size_v2(tile_coordinate_difference vector)
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
    auto dist_multiplier = convert_x_y_to_dist[(y_lookup * TLookupSize) + x_lookup];

    float x_mul = std::get<0>(dist_multiplier);
    float y_mul = std::get<1>(dist_multiplier);

    return (x * x_mul) + (y * y_mul);

}

template<uint32_t TLookupSize>
inline float distance_lookup<TLookupSize>::calculate_new_manhattan_multiplier(
    const tile_coordinate& from, 
    const tile_coordinate& to, 
    float current_manhattan_multiplier, 
    bool is_y_axis_neighbour)
{
    //get difference
    auto cord_diff = from - to;

    //remove sign values and convert biggest to 
    uint32_t x = static_cast<uint32_t>(abs(cord_diff.x));
    uint32_t y = static_cast<uint32_t>(abs(cord_diff.y));

    //get the bit shift amount
    int leading_zeros = clz(x | y);

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

    auto distance_difference_to_neighbour = mul_change_to_neighbour[lookup_index];

    //scale it down based on how much the input coordinate had to be scaled 
    distance_difference_to_neighbour = distance_difference_to_neighbour / (1 << bit_shift_amount);

    return current_manhattan_multiplier + distance_difference_to_neighbour;
}
