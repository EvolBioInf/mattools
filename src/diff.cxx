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

matrix diff(const matrix &self, const matrix &other)
{
	auto new_names = common_names(self.get_names(), other.get_names());
	auto new_self = sample2(self, new_names.begin(), new_names.end());
	auto new_other = sample2(other, new_names.begin(), new_names.end());

	auto size = new_names.size();
	auto ret = matrix(new_names, std::vector<double>(size * size));

	auto self_it = begin_square(new_self);
	auto other_it = begin_square(new_other);

	for (auto it = begin_square(ret); it != end_square(ret); it++) {
		*it = *self_it++ - *other_it++;
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
