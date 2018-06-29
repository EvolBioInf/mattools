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
#include <cmath>
#include <cstdio>
#include <err.h>
#include <getopt.h>
#include <iostream>
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
	static struct option long_options[] = {
		{"help", no_argument, 0, 0},   // print help
		{"full", no_argument, 0, 'f'}, // full matrix
		{0, 0, 0, 0}				   //
	};

	bool full_matrix = false;

	while (true) {
		int long_index;
		int c = getopt_long(argc, argv, "f", long_options, &long_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: // long option
			{
				auto option_string = std::string(long_options[long_index].name);
				if (option_string == "help") {
					mat_compare_usage(EXIT_SUCCESS);
				}
				break;
			}
			case 'f': full_matrix = true; break;
			default: /* intentional fall-through */
			case '?': mat_compare_usage(EXIT_FAILURE);
		}
	}

	argc -= optind, argv += optind; // hack

	auto matrices = parse_all(argv);

	if (matrices.size() < 2) {
		errx(1, "At least two matrices must be provided.");
	}

	if (!full_matrix) {
		std::cout << p2_norm(matrices[0], matrices[1]) << std::endl;
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
					p2_norm(matrices[i], matrices[j]);
			}
		}

		std::cout << cmpmat.to_string();
	}

	return 0;
}

static void mat_compare_usage(int status)
{
	static const char str[] = {
		"usage: mat compare [OPTIONS] [FILE...]\n" // this comment is a hack
		"Compute euclidean distance of two distances matrices.\n\n"
		"Available options:\n"
		" -f, --full          output a full distance matrix\n"
		"     --help          print this help\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}
