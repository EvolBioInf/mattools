/*
 * Copyright (C) 2018  Fabian Klötzl
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
#include <cmath>
#include <cstdio>
#include <err.h>
#include <getopt.h>
#include <iostream>
#include <random>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "matrix.h"

template <class T> class lower_tr_iterator
{
  public:
	using my_type = lower_tr_iterator<T>;
	using matrix_type = T;
	using size_type = matrix::size_type;
	using value_type = double;
	using reference =
		typename std::conditional<std::is_const<matrix_type>::value,
								  const double &, double &>::type;
	using pointer = typename std::conditional<std::is_const<matrix_type>::value,
											  const double *, double *>::type;

	using iterator_category = std::bidirectional_iterator_tag;
	using difference_type = ssize_t;

  private:
	matrix_type *data = nullptr;
	size_type row = 0;
	size_type col = 0;

	lower_tr_iterator(matrix_type &_data, size_type _row, size_type _col)
		: data(&_data), row(_row), col(_col)
	{
	}

	reference value() const
	{
		return data->entry(row, col);
	}

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

  public:
	lower_tr_iterator(){};

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
		return data == other.data && row == other.row && col == other.col;
	}

	bool operator!=(my_type other) const
	{
		return !(*this == other);
	}

	bool operator<(my_type other) const
	{
		if (row < other.row) {
			return true;
		} else if (row == other.row) {
			return col < other.col;
		} else {
			return false;
		}
	}

	static auto begin(matrix_type &self)
	{
		return my_type(self, 1, 0);
	}

	static auto end(matrix_type &self)
	{
		auto size = self.get_size();
		return my_type(self, size, 0);
	}
};

template <typename T> auto begin_lower_triangle(T &self)
{
	return lower_tr_iterator<T>::begin(self);
}

template <typename T> auto end_lower_triangle(T &self)
{
	return lower_tr_iterator<T>::end(self);
}

/**
 * @brief A vector maps an index to a value. Invert this relation.
 *
 * @returns An unordered map with inverse relationship.
 */
template <typename T>
auto make_index_map(const std::vector<T> &container)
	-> std::unordered_map<T, size_t>
{
	auto ret = std::unordered_map<T, size_t>{};
	ret.reserve(container.size());
	for (size_t i = 0; i < container.size(); ++i) {
		ret[container[i]] = i;
	}
	return ret;
}

template <typename T>
auto make_random_index_map(const std::vector<T> &container)
	-> std::unordered_map<T, size_t>
{
	auto indices = std::vector<size_t>();
	std::iota(begin(indices), end(indices), 0);
	std::random_device rd{};
	std::shuffle(begin(indices), end(indices), rd);

	auto ret = std::unordered_map<T, size_t>{};
	ret.reserve(container.size());
	for (size_t i = 0; i < container.size(); ++i) {
		ret[container[i]] = i;
	}
	return ret;
}

double lower_triangle_avg(const matrix &self)
{
	double ret = 0;

	auto size = self.get_size();
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < i; j++) {
			ret += self.entry(i, j);
		}
	}

	double count = size * (size - 1) / 2.0;
	return ret / count;
}

double lower_triangle_avg_it(const matrix &self)
{
	auto ret = std::accumulate(begin_lower_triangle(self),
							   end_lower_triangle(self), 0.0);

	auto size = self.get_size();
	double count = size * (size - 1) / 2.0;
	return ret / count;
}

double lower_triangle_stddvt(const matrix &self, double avg)
{
	auto sum =
		std::accumulate(begin_lower_triangle(self), end_lower_triangle(self),
						0.0, [=](double sum, double value) {
							auto t = value - avg;
							return sum + t * t;
						});

	auto size = self.get_size();
	double count = size * (size - 1) / 2.0 - 1;
	return std::sqrt(sum / count);
}

matrix normalize(const matrix &self)
{
	auto ret = self;
	auto avg = lower_triangle_avg(self);
	auto sd = lower_triangle_stddvt(self, avg);
	auto begin = begin_lower_triangle(ret);
	auto end = end_lower_triangle(ret);

	// for (auto it = begin; it != end; it++) {
	// 	*it = (*it - avg) / sd;
	// }

	for (size_t i = 0; i < ret.get_size(); i++) {
		for (size_t j = 0; j < ret.get_size(); j++) {
			ret.entry(i, j) = (ret.entry(i, j) - avg) / sd;
		}
	}

	return ret;
}

double Z(const matrix &self, const matrix &other)
{
	const auto &self_names = self.get_names();
	auto self_map = make_index_map(self_names);
	auto other_map = make_index_map(other.get_names());

	// compute the intersection of names
	auto common_names = std::vector<std::string>{};
	for (const auto &name : self_names) {
		if (other_map.find(name) != other_map.end()) {
			common_names.push_back(name);
		}
	}

	double dist = 0;

	for (size_t i = 0; i < common_names.size() - 1; i++) {
		const auto &name1 = common_names[i];
		for (size_t j = i + 1; j < common_names.size(); j++) {
			const auto &name2 = common_names[j];
			auto d1 = self.entry(name1, name2);
			auto d2 = other.entry(name1, name2);
			dist += d1 * d2;
		}
	}

	return dist;
}

double mantel(const matrix &self, const matrix &other, bool donormalize)
{

	auto self_names = self.get_names();
	auto other_names = other.get_names();
	std::sort(begin(self_names), end(self_names));
	std::sort(begin(other_names), end(other_names));
	auto common_names = std::vector<std::string>{};
	std::set_intersection(begin(self_names), end(self_names),
						  begin(other_names), end(other_names),
						  std::back_inserter(common_names));

	auto size = common_names.size();

	auto new_self = sample2(self, common_names.begin(), common_names.end());
	auto new_other = sample2(other, common_names.begin(), common_names.end());

	if (donormalize) {
		new_self = normalize(new_self);
		new_other = normalize(new_other);
	}

	auto orig = Z(new_self, new_other);
	std::cerr << "orig: " << orig << std::endl;

	auto montecarlo = std::vector<double>();
	auto indices = std::vector<size_t>(size);
	std::iota(begin(indices), end(indices), 0);

	std::random_device rd{};
	std::mt19937 g(rd());

	for (size_t runs = 0; runs < 100000; runs++) {
		std::shuffle(begin(indices), end(indices), g);
		// randomize
		double dist = 0;

		for (size_t i = 0; i < size - 1; i++) {
			for (size_t j = i + 1; j < size; j++) {
				auto d1 = new_self.entry(i, j);
				auto d2 = new_other.entry(indices[i], indices[j]);
				dist += d1 * d2;
			}
		}

		std::cerr << "Z: " << dist << std::endl;
		montecarlo.push_back(dist);
	}

	std::sort(begin(montecarlo), end(montecarlo));
	auto it = std::lower_bound(begin(montecarlo), end(montecarlo), orig);

	return (end(montecarlo) - it) / (double)montecarlo.size();
}

static void mat_mantel_usage(int status);

/**
 * @brief The main function of `mat compare`.
 *
 * @param argc - It's argc, stupid. (Advanced by one from global argc.)
 * @param argv - It's argv, stupid. (Advanced by one from global argv.)
 * @returns 0 iff successful.
 */
int mat_mantel(int argc, char **argv)
{
	static struct option long_options[] = {
		{"help", no_argument, 0, 0},   // print help
		{"full", no_argument, 0, 'f'}, // full matrix
		{"normalize", no_argument, 0, 'n'},
		{0, 0, 0, 0} //
	};

	bool full_matrix = false;
	bool donormalize = false;

	while (true) {
		int long_index;
		int c = getopt_long(argc, argv, "fn", long_options, &long_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: // long option
			{
				auto option_string = std::string(long_options[long_index].name);
				if (option_string == "help") {
					mat_mantel_usage(EXIT_SUCCESS);
				}
				break;
			}
			case 'f': full_matrix = true; break;
			case 'n': donormalize = true; break;
			default: /* intentional fall-through */
			case '?': mat_mantel_usage(EXIT_FAILURE);
		}
	}

	argc -= optind, argv += optind; // hack

	auto matrices = parse_all(argv);

	if (matrices.size() < 2) {
		errx(1, "At least two matrices must be provided.");
	}

	if (!full_matrix) {
		std::cout << mantel(matrices[0], matrices[1], donormalize) << std::endl;
	} else {
		// compute a full distance matrix
		auto size = matrices.size();
		auto names = std::vector<std::string>();
		auto values = std::vector<double>(size * size);
		names.reserve(size);

		// come up with a new name for each matrix
		for (size_t i = 1; i <= size; i++) {
			names.push_back(std::string("M") + std::to_string(i));
		}

		auto cmpmat = matrix{names, values};
		for (size_t i = 0; i < size; i++) {
			for (size_t j = 0; j < size; j++) {
				if (i == j) continue;
				cmpmat.entry(i, j) = cmpmat.entry(j, i) =
					mantel(matrices[i], matrices[j], donormalize);
			}
		}

		std::cout << cmpmat.to_string();
	}

	return 0;
}

static void mat_mantel_usage(int status)
{
	static const char str[] = {
		"usage: mat mantel [OPTIONS] [FILE...]\n" // this comment is a hack
		"Compare matrices using the mantel test.\n\n"
		"Available options:\n"
		" -f, --full          output a full distance matrix\n"
		"     --help          print this help\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}
