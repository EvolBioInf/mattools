/*
 * Copyright (C) 2017-2019  Fabian Klötzl
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
#include <experimental/array>
#include <getopt.h>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>
#include "matrix.h"

// not exposed, yet
double p1_norm(const matrix &self, const matrix &other)
{
	auto new_names = common_names(self.get_names(), other.get_names());

	auto new_self = sample2(self, new_names.begin(), new_names.end());
	auto new_other = sample2(other, new_names.begin(), new_names.end());

	double dist = 0;
	auto other_it = begin_lower_triangle(new_other);

	for (auto entry : lower_triangle(new_self)) {
		dist += std::fabs(entry - *other_it++);
	}

	return dist;
}

/**
 * @brief Treat two distance matrices as vectors and compute their euklidian
 * distance. To avoid errors from different arrangements, the set of common
 * names is computed first and then the corresponding sub matrices are used.
 *
 * @param self - One matrix
 * @param other - The other matrix, duh.
 * @returns the euklidian distance.
 */
double p2_norm(const matrix &self, const matrix &other)
{
	auto new_names = common_names(self.get_names(), other.get_names());

	auto new_self = sample2(self, new_names.begin(), new_names.end());
	auto new_other = sample2(other, new_names.begin(), new_names.end());

	double dist = 0;
	auto other_it = begin_lower_triangle(new_other);

	for (auto entry : lower_triangle(new_self)) {
		auto d = entry - *other_it++;
		dist += d * d;
	}

	auto n = new_names.size() * (new_names.size() - 1) / 2;

	return std::sqrt(dist / n);
}

double rel(const matrix &self, const matrix &other)
{
	auto new_names = common_names(self.get_names(), other.get_names());

	auto new_self = sample2(self, new_names.begin(), new_names.end());
	auto new_other = sample2(other, new_names.begin(), new_names.end());

	double dist = 0;
	auto other_it = begin_lower_triangle(new_other);

	for (auto entry : lower_triangle(new_self)) {
		auto d = 2 * (entry - *other_it);
		auto f = entry + *other_it;
		dist += std::abs(d / f);
		other_it++;
	}

	auto n = new_names.size() * (new_names.size() - 1) / 2;

	return dist / n;
}

using function_type = double(double, double);

template <function_type numerator_fn, function_type denominator_fn>
double delta(const matrix &self, const matrix &other)
{
	auto new_names = common_names(self.get_names(), other.get_names());

	auto new_self = sample2(self, new_names.begin(), new_names.end());
	auto new_other = sample2(other, new_names.begin(), new_names.end());

	double dist = 0;
	auto other_it = begin_lower_triangle(new_other);

	for (auto entry : lower_triangle(new_self)) {
		auto numerator = numerator_fn(entry, *other_it);
		auto denominator = denominator_fn(entry, *other_it);

		dist += numerator / denominator;
		other_it++;
	}

	return dist;
}

double just_Dij(double Dij, double)
{
	return Dij;
}

double just_Dij_squared(double Dij, double)
{
	return Dij * Dij;
}

double difference_squared(double Dij, double dij)
{
	return (Dij - dij) * (Dij - dij);
}

double average_squared(double Dij, double dij)
{
	auto temp = (Dij + dij) / 2.0;
	return temp * temp;
}

double just_average(double Dij, double dij)
{
	auto temp = (Dij + dij) / 2.0;
	return temp;
}

double just_one(double, double)
{
	return 1.0;
}

double difference_abs(double Dij, double dij)
{
	return std::fabs(Dij - dij);
}

double hausdorff(const matrix &self, const matrix &other)
{
	auto new_names = common_names(self.get_names(), other.get_names());

	auto new_self = sample2(self, new_names.begin(), new_names.end());
	auto new_other = sample2(other, new_names.begin(), new_names.end());

	auto my_max = [](double a, double b) { return std::max(a, b); };
	auto other_it = begin_lower_triangle(new_other);

	double dist = std::inner_product(begin_lower_triangle(new_self),
									 end_lower_triangle(new_self), other_it,
									 0.0, my_max, difference_abs);

	return dist;
}

const auto delta1 = delta<difference_squared, just_Dij_squared>;
const auto delta2 = delta<difference_squared, average_squared>;
const auto delta3 = delta<difference_squared, just_Dij>;
const auto delta4 = delta<difference_squared, just_average>;
const auto delta5 = delta<difference_abs, just_average>;
const auto delta6 = delta<difference_squared, just_one>;

static void mat_compare_usage(int status);

/**
 * @brief The main function of `mat compare`.
 *
 * @param argc - It's argc, stupid. (Advanced by one from global argc.)
 * @param argv - It's argv, stupid. (Advanced by one from global argv.)
 * @returns 0 iff successful.
 */
int mat_compare(int argc, char **argv)
{
	int fn_index;
	auto functions = std::experimental::make_array( //
		delta1, delta2, delta3, delta4, delta5, rel, delta6, hausdorff);

	static struct option long_options[] = {
		{"delta1", no_argument, &fn_index, 0},
		{"delta2", no_argument, &fn_index, 1},
		{"delta3", no_argument, &fn_index, 2},
		{"delta4", no_argument, &fn_index, 3},
		{"delta5", no_argument, &fn_index, 4},
		{"delta6", no_argument, &fn_index, 6},
		{"hausdorff", no_argument, &fn_index, 7},
		{"help", no_argument, 0, 0},
		{"rel", no_argument, &fn_index, 5},
		{0, 0, 0, 0} //
	};

	while (true) {
		int long_index;
		int c = getopt_long(argc, argv, "", long_options, &long_index);

		if (c == -1) {
			break;
		}

		if (c == 0) { // long option
			auto option_string = std::string(long_options[long_index].name);
			if (option_string == "help") {
				mat_compare_usage(EXIT_SUCCESS);
			}
			// fn_index is set by getopt_long
		} else {
			mat_compare_usage(EXIT_FAILURE);
		}
	}

	argc -= optind, argv += optind; // hack

	if (argc < 2) mat_compare_usage(EXIT_FAILURE);

	auto first_file_name = std::string(argv[0]);
	auto second_file_name = std::string(argv[1]);

	auto first_matrices = parse(first_file_name);
	auto second_matrices = parse(second_file_name);

	// check first and second matrices
	auto count = std::min(first_matrices.size(), second_matrices.size());
	for (size_t i = 0; i < count; i++) {
		auto fn = functions[fn_index];
		double dist = fn(first_matrices[i], second_matrices[i]);
		std::cout << dist << std::endl;
	}

	return 0;
}

static void mat_compare_usage(int status)
{
	static const char str[] = {
		"usage: mat compare [OPTIONS] FILE1 FILE2\n" // this comment is a hack
		"Measure the distance of distances matrices from two files.\n\n"
		"Available options:\n"
		"  --delta1        Compute directed Fitch-Margoliash distance\n"
		"  --delta2        Compute undirected Fitch-Margoliash distance\n"
		"  --delta3        \n"
		"  --delta4        \n"
		"  --delta5        \n"
		"  --hausdorff     Find the biggest absolute difference\n"
		"  --help          Print this help\n"
		"  --rel           Compute the average relative dissimilarity\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}
