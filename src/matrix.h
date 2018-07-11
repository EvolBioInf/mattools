/*
 * Copyright (C) 2017  Fabian Klötzl
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cassert>
#include <err.h>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>
// #include <optional>

class matrix
{
  public:
	using size_type = size_t;

  protected:
	template <typename T>
	auto make_index_map(const std::vector<T> &container)
		-> std::unordered_map<T, size_type>
	{
		auto ret = std::unordered_map<T, size_type>{};
		ret.reserve(container.size());
		for (size_type i = 0; i < container.size(); ++i) {
			ret[container[i]] = i;
		}
		return ret;
	}
	/// The size of the matrix
	size_type size = 0;
	/// A list of names (of sequences)
	std::vector<std::string> names = {};
	/// The matrix itself
	std::vector<double> values = {};

	std::unordered_map<std::string, size_t> name_map{};

	bool m_has_coverages = false;
	std::vector<double> coverages = {};

  public:
	matrix() = default;
	/** @brief Create a new matrix from a set of names and values.
	 *
	 * Take arguments by value, because compilers can optimize this anyway.
	 *
	 * @param _names - The new set of names
	 * @param _values - The new values
	 * @returns the new matrix
	 */
	matrix(std::vector<std::string> _names, std::vector<double> _values)
		: size{_names.size()}, names{std::move(_names)}, values{std::move(
															 _values)},
		  name_map{make_index_map(names)}, m_has_coverages{false}, coverages{}

	{
		assert(size == names.size());
		assert(size * size == values.size());
	}

	matrix(std::vector<std::string> _names, std::vector<double> _values,
		   std::vector<double> _coverages)
		: size{_names.size()}, names{std::move(_names)},
		  values{std::move(_values)}, name_map{make_index_map(names)},
		  m_has_coverages{true}, coverages{std::move(_coverages)}

	{
		assert(size == names.size());
		assert(size * size == values.size());
		assert(size * size == coverages.size());
	}

	double &entry(const std::string &ni, const std::string &nj)
	{
		return entry(name_map.at(ni), name_map.at(nj));
	}

	const double &entry(const std::string &ni, const std::string &nj) const
	{
		return entry(name_map.at(ni), name_map.at(nj));
	}

	/** @brief Access an entry by coordinates.
	 *
	 * @param i - the row index
	 * @param j - the column index
	 * @returns a mutable reference to the entry
	 */
	double &entry(size_type i, size_type j)
	{
		return values[i * size + j];
	}

	/** @brief Access an entry by coordinates.
	 *
	 * @param i - the row index
	 * @param j - the column index
	 * @returns an immutable reference to the entry
	 */
	const double &entry(size_type i, size_type j) const
	{
		return values[i * size + j];
	}

	/** @brief Get an iterator to the beginning of a row.
	 *
	 * @param i - the row index
	 * @returns An iterator.
	 */
	auto row(size_type i) noexcept
	{
		return values.begin() + i * size;
	}

	/** @brief Get an iterator to the beginning of a row.
	 *
	 * @param i - the row index
	 * @returns A read-only iterator.
	 */
	auto row(size_type i) const noexcept
	{
		return values.cbegin() + i * size;
	}

	/** @brief Get an iterator to the end of a row.
	 *
	 * @param i - the row index
	 * @returns An iterator.
	 */
	auto row_end(size_type i) noexcept
	{
		return values.begin() + (i + 1) * size;
	}

	/** @brief Get an iterator to the beginning of a row.
	 *
	 * @param i - the row index
	 * @returns A read-only iterator.
	 */
	auto row_end(size_type i) const noexcept
	{
		return values.cbegin() + (i + 1) * size;
	}

	/** @brief Get the size of the matrix
	 *
	 * @returns The size.
	 */
	auto get_size() const noexcept
	{
		return size;
	}

	/** @brief Get the list of names
	 *
	 * @returns Returns a read-only reference of the names.
	 */
	auto get_names() const noexcept -> const std::vector<std::string> &
	{
		return names;
	}

	auto get_values() const noexcept -> const std::vector<double> &
	{
		return values;
	}

	auto has_coverages() const noexcept
	{
		return m_has_coverages;
	}

	auto get_coverages() const -> const std::vector<double> &
	{
		if (!has_coverages()) throw "no coverages.";
		return coverages;
	}

	double &cov_entry(size_type i, size_type j)
	{
		if (!has_coverages()) throw "no coverages.";
		return coverages[i * size + j];
	}

	const double &cov_entry(size_type i, size_type j) const
	{
		if (!has_coverages()) throw "no coverages.";
		return coverages[i * size + j];
	}

	// defined in matrix.cxx
	std::string to_string() const;
};

/** @brief Sample a distance matrix. Only the names given by the input list are
 * included in the new submatrix. The implied order of names from the index list
 * is preserved.
 *
 * @param self - The matrix to sample.
 * @param first - An iterator to a list of names.
 * @param last - An iterator past the list of names.
 * @returns a new matrix with only the names specified in the list.
 */
template <typename ForwardIt>
matrix sample2(const matrix &self, const ForwardIt first, const ForwardIt last)
{
	auto new_size = distance(first, last);
	auto new_names = std::vector<std::string>(first, last);

	auto new_values = std::vector<double>(new_size * new_size);
	auto ret = matrix{new_names, new_values};

	// Copy and rearrange the values
	for (auto it = first; it != last; it++) {
		auto name1 = *it;
		for (auto jk = first; jk != last; jk++) {
			auto name2 = *jk;
			ret.entry(name1, name2) = self.entry(name1, name2);
		}
	}

	return ret;
}

// defined in matrix.cxx
std::vector<matrix> parse_all(const char *const *);
std::vector<matrix> parse_all(std::vector<std::string> file_names);
std::string format(const matrix &, char, const char *, bool);

class square_iterator_helper
{
  public:
	using size_type = matrix::size_type;
	size_type row = 0;
	size_type col = 0;
	size_type size = 0;

	void next()
	{
		col++;
		if (col >= size) {
			col = 0;
			row++;
		}
	}

	void prev()
	{
		if (col > 0) {
			col--;
		} else {
			row--;
			col = size;
		}
	}

	static square_iterator_helper begin(size_type size)
	{
		return square_iterator_helper(0, 0, size);
	}

	static square_iterator_helper end(size_type size)
	{
		return square_iterator_helper(size, 0, size);
	}

  public:
	square_iterator_helper();
	square_iterator_helper(size_type _row, size_type _col, size_type _size)
		: row(_row), col(_col), size(_size){};
};

class ltriangle_iterator_helper
{
  public:
	using size_type = matrix::size_type;
	size_type row = 0;
	size_type col = 0;
	size_type size = 0;

	void next()
	{
		col++;
		if (col >= row) {
			col = 0;
			row++;
		}
	}

	void prev()
	{
		if (col > 0) {
			col--;
		} else {
			row--;
			col = row - 1;
		}
	}

	static ltriangle_iterator_helper begin(size_type size)
	{
		return ltriangle_iterator_helper(1, 0, size);
	}

	static ltriangle_iterator_helper end(size_type size)
	{
		return ltriangle_iterator_helper(size, 0, size);
	}

  public:
	ltriangle_iterator_helper();
	ltriangle_iterator_helper(size_type _row, size_type _col, size_type _size)
		: row(_row), col(_col), size(_size){};
};

template <class Matrix, class Helper>
class matrix_iterator
{
  private:
	using my_double = typename std::conditional<std::is_const<Matrix>::value,
												const double, double>::type;

  public:
	using my_type = matrix_iterator<Matrix, Helper>;
	using matrix_type = Matrix;
	using size_type = matrix::size_type;
	using value_type = double;
	using reference = my_double &;
	using pointer = my_double *;

	using iterator_category = std::bidirectional_iterator_tag;
	using difference_type = ssize_t;

  private:
	matrix_type *data = nullptr;
	Helper helper = {};

	matrix_iterator(matrix_type &_data, Helper start)
		: data(&_data), helper(std::move(start))
	{
	}

	reference value() const
	{
		return data->entry(helper.row, helper.col);
	}

	void next()
	{
		helper.next();
	}

	void prev()
	{
		helper.prev();
	}

  public:
	matrix_iterator(){};

	reference operator*() const
	{
		return value();
	}

	auto &operator++()
	{
		next();
		return *this;
	}

	auto operator++(int)
	{
		auto ret = *this;
		next();
		return ret;
	}

	auto &operator--()
	{
		prev();
		return *this;
	}

	auto operator--(int)
	{
		auto ret = *this;
		prev();
		return ret;
	}

	bool operator==(my_type other) const
	{
		return data == other.data && helper.row == other.helper.row &&
			   helper.col == other.helper.col;
	}

	bool operator!=(my_type other) const
	{
		return !(*this == other);
	}

	bool operator<(my_type other) const
	{
		if (helper.row < other.helper.row) {
			return true;
		} else if (helper.row == other.helper.row) {
			return helper.col < other.helper.col;
		} else {
			return false;
		}
	}

	static auto begin(matrix_type &self)
	{
		return my_type(self, Helper::begin(self.get_size()));
	}

	static auto end(matrix_type &self)
	{
		return my_type(self, Helper::end(self.get_size()));
	}
};

template <typename T>
auto begin_square(T &self)
{
	return matrix_iterator<T, square_iterator_helper>::begin(self);
}

template <typename T>
auto end_square(T &self)
{
	return matrix_iterator<T, square_iterator_helper>::end(self);
}

template <typename T>
auto begin_lower_triangle(T &self)
{
	return matrix_iterator<T, ltriangle_iterator_helper>::begin(self);
}

template <typename T>
auto end_lower_triangle(T &self)
{
	return matrix_iterator<T, ltriangle_iterator_helper>::end(self);
}

template <typename T>
class square_agg
{
	T &mat;

  public:
	square_agg(T &_mat) : mat(_mat)
	{
	}

	auto begin()
	{
		return matrix_iterator<T, square_iterator_helper>::begin(mat);
	}

	auto end()
	{
		return matrix_iterator<T, square_iterator_helper>::end(mat);
	}
};

template <typename T>
auto square(T &self)
{
	return square_agg<T>(self);
}

template <typename T>
class lower_triangle_agg
{
	T &mat;

  public:
	lower_triangle_agg(T &_mat) : mat(_mat)
	{
	}

	auto begin()
	{
		return matrix_iterator<T, ltriangle_iterator_helper>::begin(mat);
	}

	auto end()
	{
		return matrix_iterator<T, ltriangle_iterator_helper>::end(mat);
	}
};

template <typename T>
auto lower_triangle(T &self)
{
	return lower_triangle_agg<T>(self);
}

inline std::vector<std::string>
common_names(std::vector<std::string> self_names,
			 std::vector<std::string> other_names)
{
	auto ret = std::vector<std::string>{};

	std::sort(self_names.begin(), self_names.end());
	std::sort(other_names.begin(), other_names.end());
	std::set_intersection(self_names.begin(), self_names.end(),
						  other_names.begin(), other_names.end(),
						  std::back_inserter(ret));

	return ret;
}
