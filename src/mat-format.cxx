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
#include <cassert>
#include <cmath>
#include <err.h>
#include <errno.h>
#include <getopt.h>
#include <iostream>
#include <numeric>
#include <regex>
#include <string>
#include <vector>
#include "matrix.h"

static void mat_format_usage(int);
static char unescape(const char *);

static double PRECISION = 0.05;
static const /*expr*/ auto close_enough = [](double a, double b) {
	return a * (1.0 - PRECISION) <= b && b <= a * (1.0 + PRECISION);
};

/**
 * @brief Rearrange matrix by sorted names
 *
 * @param self - the matrix to rearrange
 * @returns sorted matrix
 */
static matrix sort(const matrix &self)
{
	auto names = self.get_names();
	std::sort(begin(names), end(names));

	return sample2(self, begin(names), end(names));
}

/**
 * @brief Unescape a sequence to its corresponding character. Why is there no
 * standard function for this?
 *
 * @param str - The c-style escape sequence.
 * @returns the corresponding character.
 */
static char unescape(const char *str)
{
	char ret = '?';
	if (str[0] != '\\') return str[0];

	assert(str[0] == '\\');
	switch (str[1]) {
		case '\'': ret = '\''; break;
		case '"': ret = '"'; break;
		case '\\': ret = '\\'; break;
		case 'a': ret = '\a'; break;
		case 'b': ret = '\b'; break;
		case 'f': ret = '\f'; break;
		case 'n': ret = '\n'; break;
		case 'r': ret = '\r'; break;
		case 't': ret = '\t'; break;
		case 'v': ret = '\v'; break;
		default: ret = '?'; break;
	}

	return ret;
}

/**
 * @brief Validate that the matrix is a proper distance matrix, and fix issues
 * where possible.
 *
 * @param self - The matrix to validate.
 * @returns the fixed matrix.
 */
static matrix fix(const matrix &original)
{
	auto self = matrix(original);
	auto size = self.get_size();

	// check positivity
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			if (self.entry(i, j) < 0) {
				warnx("Fixed entry (%zu,%zu); was negative: %lf, now 0.", i, j,
					  self.entry(i, j));
				self.entry(i, j) = 0.0;
			}
		}
	}

	// check main diagonal
	for (size_t i = 0; i < size; i++) {
		if (self.entry(i, i) != 0) {
			warnx("Fixed entry (%zu,%zu); was %lf, now is 0.", i, i,
				  self.entry(i, i));
			self.entry(i, i) = 0.0;
		}
	}

	// check symmetry
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < i; j++) {
			if (!close_enough(self.entry(i, j), self.entry(j, i))) {
				warnx("Fixed asymmetric cells (%zu,%zu) and (%zu,%zu); entries "
					  "are now averaged.",
					  i, j, j, i);

				auto avg = (self.entry(i, j) + self.entry(j, i)) / 2.0;
				self.entry(i, j) = self.entry(j, i) = avg;
			}
		}
	}

	return self;
}

/**
 * @brief Validate that the matrix is a proper distance matrix. Errors are
 * non-recoverable.
 *
 * @param self - The matrix to validate.
 * @param truncate_names - True iff names are truncated.
 * @returns the fixed matrix.
 */
static matrix validate(const matrix &original, bool truncate_names)
{
	auto self = matrix{original};
	auto size = self.get_size();

	// maybe only check first ten chars?
	auto equal = [=](const std::string &a, const std::string &b) {
		return truncate_names ? a.compare(0, 10, b, 0, 10) == 0 : a == b;
	};

	// check name uniqueness
	auto names_copy = self.get_names();
	std::sort(names_copy.begin(), names_copy.end());
	for (size_t i = 0; i < size - 1; i++) {
		if (equal(names_copy[i], names_copy[i + 1])) {
			// I don't know what to do — panic!
			auto str = truncate_names
						   ? "The truncated name %.10s appears twice."
						   : "The name %s appears twice.";
			errx(1, str, names_copy[i].c_str());
		}
	}

	// check nan and zero beyond main diagonal
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < i; j++) {
			if (close_enough(self.entry(i, j), 0)) {
				errx(1, "Zero entry beyond the main diagonal (%zu,%zu).", i, j);
			}
			if (std::isnan(self.entry(i, j))) {
				errx(1, "Not a Number (%zu,%zu)", i, j);
			}
		}
	}

	// check triangle inequality
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < i; j++) {
			for (size_t k = 0; k < j; k++) {
				// ensure d(i,j) ≤ d(i,k) + d(k,j)
				if (self.entry(i, j) > self.entry(i, k) + self.entry(j, k) &&
					!close_enough(self.entry(i, j),
								  self.entry(i, k) + self.entry(j, k))) {
					// panic
					errx(1,
						 "Violation of triangle inequality for (%zu,%zu) "
						 "and (%zu,%zu)+(%zu,%zu)",
						 i, j, i, k, k, j);
				}
			}
		}
	}

	return self;
}

/**
 * @brief The main function of `mat format`.
 *
 * @param argc - It's argc, stupid. (Advanced by one from global argc.)
 * @param argv - It's argv, stupid. (Advanced by one from global argv.)
 * @returns 0 iff successful.
 */
int mat_format(int argc, char **argv)
{
	static struct option long_options[] = {
		{"help", no_argument, 0, 'h'}, // print help
		{"truncate-names", no_argument, 0, 0},
		{"sort", no_argument, 0, 's'}, // sort by name
		{"validate", no_argument, 0, 'v'},
		{"fix", no_argument, 0, 'f'},
		{"precision", required_argument, 0, 0},
		{"separator", required_argument, 0, 0},
		{"format", required_argument, 0, 0},
		// {"apply-jc", no_argument, 0, 0},
		// {"unapply-jc", no_argument, 0, 0},
		{0, 0, 0, 0} //
	};

	auto fix_flag = false;
	auto format_flag = false;
	auto format_specifier = "%9.3e";
	auto separator = ' ';
	auto sort_flag = false;
	auto truncate_names = false;
	auto validate_flag = false;

	while (true) {
		int long_index;
		int c = getopt_long(argc, argv, "fhsv", long_options, &long_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: // long option
			{
				auto option_str = std::string{long_options[long_index].name};
				if (option_str == "separator") {
					separator = unescape(optarg);
					format_flag = true;
					break;
				}

				if (option_str == "format") {
					/* GCC's libstdc++ has a bug here. Beware, when rewriting
					the regex.
					https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77469 */
					auto rvalid = std::regex("^%([#0 +\\-]*)"
											 "(\\-?[1-9][0-9]*)?(\\.[0-9]*)?"
											 "([eE]|[lL]?[fF])$");

					format_specifier = optarg;
					if (!std::regex_match(format_specifier, rvalid)) {
						errno = EINVAL;
						err(errno, "invalid format specifier: %s", optarg);
					}

					format_flag = true;
					break;
				}

				if (option_str == "precision") {
					auto prec = std::stod(optarg);
					PRECISION = prec;
					break;
				}

				if (option_str == "truncate-names") {
					truncate_names = true;
					format_flag = true;
					break;
				}
				break;
			}
			case 'f': fix_flag = true; break;
			case 'h': mat_format_usage(EXIT_SUCCESS);
			case 'v':
				validate_flag = true;
				fix_flag = true;
				break;
			case 's': sort_flag = true; break;
			case '?': // intentional fall-through
			default: mat_format_usage(EXIT_FAILURE);
		}
	}

	argc -= optind, argv += optind;

	auto matrices = parse_all(argv);

	for (auto &m : matrices) {
		if (fix_flag) {
			m = fix(m);
		}

		if (validate_flag) {
			m = validate(m, truncate_names);
		}

		if (sort_flag) {
			m = sort(m);
		}

		auto str = format_flag
					   ? format(m, separator, format_specifier, truncate_names)
					   : m.to_string();
		std::cout << str;
	}

	return 0;
}

/**
 * @brief Print usage and exit.
 *
 * @param status - The exit code
 * @return It doesn't.
 */
static void mat_format_usage(int status)
{
	static const char str[] = {
		"usage: mat format [OPTIONS] [FILE...]\n" // this comment is a hack
		"Format the distance matrix.\n\n"
		"Available options:\n"
		"  -f, --fix             fix small errors\n"
		"      --format <str>    use <str> as the format string; default: "
		"%%1.4e\n"
		"      --precision <flt> precision to use in comparisons; default: "
		"0.05\n"
		"      --separator <c>   set the cell separator to <c>; default: ' ' "
		"aka. space\n"
		"  -s, --sort            sort by name\n"
		"      --truncate-names  truncate names to ten characters\n"
		"  -v, --validate        validate for correctness (implies -f)\n"
		"  -h, --help            print this help\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}
