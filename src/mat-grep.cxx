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
#include <cstdio>
#include <cstdlib>
#include <err.h>
#include <getopt.h>
#include <iostream>
#include <regex>
#include <string>
#include "matrix.h"

/** @brief Grep names and remove names that don't match the given pattern.
 *
 * @param self - The matrix to subsample.
 * @param rpattern - The regex pattern to search for.
 * @param invert - Iff true, invert the pattern.
 * @returns the submatrix.
 */
matrix grep(const matrix &self, const std::regex &rpattern, bool invert)
{
	auto size = self.get_size();
	const auto &old_names = self.get_names();

	auto indices = std::vector<matrix::size_type>(size);
	std::iota(begin(indices), end(indices), size_t(0));

	// remove all names that do not match the pattern
	auto split = std::remove_if(
		begin(indices), end(indices), [&](matrix::size_type index) {
			return !std::regex_search(old_names[index], rpattern) ^ invert;
		});

	// get the submatrix, so only the matching names and values remain
	return sample(self, begin(indices), split);
}

static void mat_grep_usage(int status);

/**
 * @brief The main function of `mat grep`.
 *
 * @param argc - It's argc, stupid. (Advanced by one from global argc.)
 * @param argv - It's argv, stupid. (Advanced by one from global argv.)
 * @returns 0 iff successful.
 */
int mat_grep(int argc, char **argv)
{
	if (argc < 2) {
		mat_grep_usage(EXIT_FAILURE);
	}

	static struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"invert-match", no_argument, 0, 'v'},
		{0, 0, 0, 0}};

	auto invert = false;

	while (true) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "hv", long_options, &option_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 'h': mat_grep_usage(EXIT_SUCCESS);
			case 'v': {
				invert = true;
				break;
			}
			case '?': // intentional fall-through
			default: mat_grep_usage(EXIT_FAILURE);
		}
	}

	argc -= optind, argv += optind;

	if (argc == 0) {
		errx(EXIT_FAILURE, "missing pattern");
	}

	auto pattern = std::string{argv[0]};
	argv++;
	auto rpattern = std::regex(pattern);

	auto matrices = parse_all(argv);

	for (const auto &mat : matrices) {
		std::cout << grep(mat, rpattern, invert).to_string();
	}

	return 0;
}

static void mat_grep_usage(int status)
{
	static const char str[] = {
		"usage: mat grep [OPTIONS] PATTERN [FILE...]\n"
		"Print submatrix for names matching the PATTERN.\n"
		"The PATTERN can be a regular expression using ECMAScript format.\n\n"
		"Available options:\n"
		"  -h, --help           print this help\n"
		"  -v, --invert-match   select non-matching names\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}
