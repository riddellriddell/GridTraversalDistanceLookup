#pragma once
#include <array>

#include "shared_types.h"

//helper template that give you the amount to add or subtract from an address to get the index of a neighbouring tile
//that is stored in a 1D array that maps to a 2D array
//everything is calculated at compile time and should be a cheap lookup
template<typename ToffsetType, int Twidth>
struct neighbour_address_offset
{	
	static constexpr ToffsetType left = static_cast<ToffsetType>(-1);								//left
	static constexpr ToffsetType up = -static_cast<ToffsetType>(Twidth);		//up
	static constexpr ToffsetType right = static_cast<ToffsetType>(1);								//right
	static constexpr ToffsetType down = static_cast<ToffsetType>(Twidth);		//down
	static constexpr ToffsetType up_left = static_cast<ToffsetType>(up + left);						//up left
	static constexpr ToffsetType up_right = static_cast<ToffsetType>(up + right);       			//up right
	static constexpr ToffsetType down_left = static_cast<ToffsetType>(down + left);					//down left
	static constexpr ToffsetType down_right = static_cast<ToffsetType>(down + right);				//down right
			
	constexpr static std::array<ToffsetType, static_cast<size_t>(NEIGHBOUR_DIRECTIONS::COUNT)> offset_for_direction =
		{{left, up, right, down, up_left, up_right, down_left, down_right}};

	//this helper array calculates the offset to get the tile at the opposite side of a grid
	//for example if you have the a tile on the left side of the grid and want to get at the same row
	//on the right side of the grid this will provide the offset
	//use this to get the sub sector index of a tile in a neighbouring sector given the index of an edge tile 

	static constexpr ToffsetType wrap_left =		static_cast<ToffsetType>(Twidth -1);						//left
	static constexpr ToffsetType wrap_up =			static_cast<ToffsetType>((Twidth * Twidth) - Twidth);		//up
	static constexpr ToffsetType wrap_right =		static_cast<ToffsetType>(-(Twidth-1));						//right
	static constexpr ToffsetType wrap_down =		static_cast<ToffsetType>(-((Twidth * Twidth) - Twidth));	//down
	static constexpr ToffsetType wrap_up_left =		static_cast<ToffsetType>(wrap_up + wrap_left);										//up left
	static constexpr ToffsetType wrap_up_right =	static_cast<ToffsetType>(wrap_up + wrap_right);       								//up right
	static constexpr ToffsetType wrap_down_left =	static_cast<ToffsetType>(wrap_down + wrap_left);									//down left
	static constexpr ToffsetType wrap_down_right =	static_cast<ToffsetType>(wrap_down + wrap_right);									//down right
			
	constexpr static std::array<ToffsetType, static_cast<size_t>(NEIGHBOUR_DIRECTIONS::COUNT)> wrap_around_offset_for_direction =
		{{wrap_left, wrap_up, wrap_right, wrap_down, wrap_up_left, wrap_up_right, wrap_down_left, wrap_down_right}};



	static constexpr ToffsetType x_left = static_cast<ToffsetType>(-1);											//left
	static constexpr ToffsetType x_up = static_cast<ToffsetType>(0);											//up
	static constexpr ToffsetType x_right = static_cast<ToffsetType>(1);											//right
	static constexpr ToffsetType x_down = static_cast<ToffsetType>(0);											//down
	static constexpr ToffsetType x_up_left = static_cast<ToffsetType>(-1);										//up left
	static constexpr ToffsetType x_up_right = static_cast<ToffsetType>(1);       								//up right
	static constexpr ToffsetType x_down_left = static_cast<ToffsetType>(-1);									//down left
	static constexpr ToffsetType x_down_right = static_cast<ToffsetType>(1);									//down right

	constexpr static std::array<ToffsetType, static_cast<size_t>(NEIGHBOUR_DIRECTIONS::COUNT)> x_offset_for_direction =
	{ {x_left, x_up, x_right, x_down, x_up_left, x_up_right, x_down_left, x_down_right} };

	static constexpr ToffsetType y_left = static_cast<ToffsetType>(0);											//left
	static constexpr ToffsetType y_up = static_cast<ToffsetType>(-1);											//up
	static constexpr ToffsetType y_right = static_cast<ToffsetType>(0);											//right
	static constexpr ToffsetType y_down = static_cast<ToffsetType>(1);											//down
	static constexpr ToffsetType y_up_left = static_cast<ToffsetType>(-1);										//up left
	static constexpr ToffsetType y_up_right = static_cast<ToffsetType>(-1);       								//up right
	static constexpr ToffsetType y_down_left = static_cast<ToffsetType>(1);										//down left
	static constexpr ToffsetType y_down_right = static_cast<ToffsetType>(1);									//down right

	constexpr static std::array<ToffsetType, static_cast<size_t>(NEIGHBOUR_DIRECTIONS::COUNT)> y_offset_for_direction =
	{ {y_left, y_up, y_right, y_down, y_up_left, y_up_right, y_down_left, y_down_right} };
};

template<int TSize, typename TGridDataType>
struct data_grid
{
    static constexpr int width = TSize;
    static constexpr int total_tile_count = width  * width;
    
    std::array<TGridDataType,total_tile_count> data = {};

	data_grid(TGridDataType init_value);

	static int convert_xy_to_index(const tile_coordinate& cord);

	static tile_coordinate convert_index_to_xy(int index);

	static bool is_valid(const tile_coordinate& cord);

	//TGridDataType& operator [] (tile_coordinate& cord);

	TGridDataType& operator [] (const tile_coordinate& cord);

	TGridDataType& operator [] (int index);
	
};

template <int TSize, typename TGridDataType>
data_grid<TSize, TGridDataType>::data_grid(TGridDataType init_value)
{
	for(size_t i = 0; i < data.size(); ++i)
	{
		data[i] = init_value;
	}
}

template <int TSize, typename TGridDataType>
int data_grid<TSize, TGridDataType>::convert_xy_to_index(const tile_coordinate& cord)
{
	return (cord.y * width) + cord.x;
}

template <int TSize, typename TGridDataType>
tile_coordinate data_grid<TSize, TGridDataType>::convert_index_to_xy(int index)
{
	tile_coordinate out_cord = {};

	out_cord.y = index / width;
	out_cord.x = index - (out_cord.y * width);

	return out_cord;	
}

template<int TSize, typename TGridDataType>
inline bool data_grid<TSize, TGridDataType>::is_valid(const tile_coordinate& cord)
{
	if (cord.x >= width || cord.x < 0 || cord.y >= width || cord.y < 0)
	{
		return false;
	}

	return true;
}

template <int TSize, typename TGridDataType>
TGridDataType& data_grid<TSize, TGridDataType>::operator[](const tile_coordinate& cord)
{
	int index = convert_xy_to_index(cord);
	return (*this)[index];
}


//template <int TSize, typename TGridDataType>
//TGridDataType& data_grid<TSize, TGridDataType>::operator[](tile_coordinate&& cord)
//{
//	int index = convert_xy_to_index(cord);
//	return (*this)[index];
//}

template <int TSize, typename TGridDataType>
TGridDataType& data_grid<TSize, TGridDataType>::operator[](int index)
{
	return data[index];
}
