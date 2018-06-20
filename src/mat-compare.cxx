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

// not exposed, yet
double p1_norm(const matrix &self, const matrix &other)
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
			auto d1 = self.entry(self_map[name1], self_map[name2]);
			auto d2 = other.entry(other_map[name1], other_map[name2]);
			dist += std::fabs(d1 - d2);
		}
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
	auto this_map = make_index_map(self.get_names());
	auto other_map = make_index_map(other.get_names());

	auto common_names = std::vector<std::string>{};
	for (const auto &name : self.get_names()) {
		if (other_map.find(name) != other_map.end()) {
			common_names.push_back(name);
		}
	}

	double dist = 0;

	for (size_t i = 0; i < common_names.size() - 1; i++) {
		auto name1 = common_names[i];
		for (size_t j = i + 1; j < common_names.size(); j++) {
			auto name2 = common_names[j];
			auto d1 = self.entry(this_map[name1], this_map[name2]);
			auto d2 = other.entry(other_map[name1], other_map[name2]);
			auto d = d1 - d2;
			dist += d * d;
		}
	}

	auto n = common_names.size() * (common_names.size() - 1) / 2;

	return std::sqrt(dist / n);
}

double rel(const matrix &self, const matrix &other)
{
	auto this_map = make_index_map(self.get_names());
	auto other_map = make_index_map(other.get_names());

	auto common_names = std::vector<std::string>{};
	for (const auto &name : self.get_names()) {
		if (other_map.find(name) != other_map.end()) {
			common_names.push_back(name);
		}
	}

	double dist = 0;

	for (size_t i = 0; i < common_names.size() - 1; i++) {
		auto name1 = common_names[i];
		for (size_t j = i + 1; j < common_names.size(); j++) {
			auto name2 = common_names[j];
			auto d1 = self.entry(this_map[name1], this_map[name2]);
			auto d2 = other.entry(other_map[name1], other_map[name2]);
			auto d = 2 * (d1 - d2);
			auto f = d1 + d2;
			// std::cerr << std::abs(d / f) << std::endl;
			dist += std::abs(d / f);
		}
	}

	auto n = common_names.size() * (common_names.size() - 1) / 2;

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
		std::cout << rel(matrices[0], matrices[1]) << std::endl;
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
