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

matrix diff(const matrix &self, const matrix &other)
{
	auto this_map = make_index_map(self.get_names());
	auto other_map = make_index_map(other.get_names());

	auto common_names = std::vector<std::string>{};
	for (const auto &name : self.get_names()) {
		if (other_map.find(name) != other_map.end()) {
			common_names.push_back(name);
		}
	}

	auto size = common_names.size();
	auto ret = matrix(common_names, std::vector<double>(size * size));

	for (size_t i = 0; i < common_names.size() - 1; i++) {
		auto name1 = common_names[i];
		for (size_t j = i + 1; j < common_names.size(); j++) {
			auto name2 = common_names[j];
			auto d1 = self.entry(this_map[name1], this_map[name2]);
			auto d2 = other.entry(other_map[name1], other_map[name2]);
			ret.entry(i, j) = ret.entry(j, i) = d1 - d2;
		}
	}

	return ret;
}

static void mat_diff_usage(int status);

/**
 * @brief The main function of `mat compare`.
 *
 * @param argc - It's argc, stupid. (Advanced by one from global argc.)
 * @param argv - It's argv, stupid. (Advanced by one from global argv.)
 * @returns 0 iff successful.
 */
int mat_diff(int argc, char **argv)
{
	static struct option long_options[] = {
		{"help", no_argument, 0, 0}, // print help
		{0, 0, 0, 0}				 //
	};

	bool full_matrix = false;

	while (true) {
		int long_index;
		int c = getopt_long(argc, argv, "", long_options, &long_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: // long option
			{
				auto option_string = std::string(long_options[long_index].name);
				if (option_string == "help") {
					mat_diff_usage(EXIT_SUCCESS);
				}
				break;
			}
			default: /* intentional fall-through */
			case '?': mat_diff_usage(EXIT_FAILURE);
		}
	}

	argc -= optind, argv += optind; // hack

	auto matrices = parse_all(argv);

	if (matrices.size() < 2) {
		errx(1, "At least two matrices must be provided.");
	}

	std::cout << diff(matrices[0], matrices[1]).to_string();

	return 0;
}

static void mat_diff_usage(int status)
{
	static const char str[] = {
		"usage: mat diff [OPTIONS] [FILE...]\n" // this comment is a hack
		"Compute euclidean distance of two distances matrices.\n\n"
		"Available options:\n"
		"     --help          print this help\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}
