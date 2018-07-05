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
#include "matrix.h"
#include <algorithm>
#include <err.h>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>
// #include <boost/config/warning_disable.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

/**
 * @brief Print the given matrix into a string. Allows for modified formatting.
 *
 * The default for format_specifier is chosen, so that four significant digits
 * are displayed and NaNs are right-aligned.
 *
 * @param self - The matrix to be printed.
 * @param separator - The character printed in between two cells.
 * @param format_specifier - A printf-style format specifier
 * @returns the formatted string
 */
std::string format(const matrix &self, char separator = ' ',
				   const char *format_specifier = "%9.3e",
				   bool truncate_names = false)
{
	std::string ret{};
	auto size = self.get_size();
	const auto &names = self.get_names();
	auto name_format = truncate_names ? "%-10.10s" : "%-10s";

	char buf[100];
	buf[0] = '\0';
	// a rough estimate of the resulting string size
	ret.reserve(size * 10 + size * size * 5 + size);

	ret += std::to_string(size);
	ret += "\n";

	for (size_t i = 0; i < size; i++) {
		snprintf(buf, 100, name_format, names[i].c_str());
		ret += buf;
		for (size_t j = 0; j < size; j++) {
			ret += separator;
			snprintf(buf, 100, format_specifier, self.entry(i, j));
			ret += buf;
		}
		ret += "\n";
	}

	return ret;
}

/** @brief Convert a matrix into a human-readable string (phylip format).
 *
 * @returs a string representing the matrix.
 */
std::string matrix::to_string() const
{
	return format(*this);
}

/** @brief Read a single "line" from a distance matrix.
 *
 * @param file_name - The current file, for error messages.
 * @param max_values - Read at most this many values.
 *
 * @returns a pair, with the name in the first component, followed by the
 * values.
 */
template <typename ForwardIt>
auto parse_line_spirit(const std::string &file_name, ForwardIt &first,
					   ForwardIt &last, size_t max_values)
{
	using namespace boost::spirit::x3;

	auto name = std::string{};
	auto values = std::vector<double>{};
	values.reserve(max_values);

	auto set_name = [&](const auto &ctx) { name = _attr(ctx); };
	auto push_back = [&](const auto &ctx) {
		values.push_back(_attr(ctx));
		_pass(ctx) = values.size() <= max_values;
	};

	const auto name_rule = lexeme[+graph][set_name] >> omit[*blank];
	const auto values_rule = double_[push_back] % *blank;
	const auto line_rule = name_rule >> -values_rule >> omit[eol];

	bool r = parse(first, last, line_rule);

	if (!r) {
		errx(1, "%s: parse error", file_name.c_str());
		std::terminate();
	}

	return std::make_pair(name, values);
}

/** @brief Parse a distance matrix from full or lower triangle format.
 *
 * The *phylip distance matrix* is a poorly defined file format. It exists in
 * multiple variations. These include full, lower triangle, and upper triangle
 * version. As the latter is barely used, it is ignored here. Furthermore, we
 * have to be very tolerant w.r.t. whitespace.
 *
 * As the format is poorly defined and engineered, the following situation is
 * ambiguous.
 *
 * > 1
 * > S
 * > 0
 *
 * It could mean either one distance matrix, with a sole zero on the main
 * diagonal (new lines are allowed within values). Or it could be two matrices,
 * one in lower triangular format with the diagonal values stripped and a second
 * empty matrix.
 *
 * http://evolution.genetics.washington.edu/phylip/doc/distance.html
 * https://www.mothur.org/wiki/Phylip-formatted_distance_matrix
 * http://condor.depaul.edu/waguirre/phylip_making_trees.pdf
 *
 * @param file_name - The file name, for error messages
 * @param input - The input stream to read from
 * @returns The matrix
 */
template <typename InputIt>
matrix parse_tolerant_internal(const std::string &file_name, InputIt &first,
							   InputIt &last)
{
	using namespace boost::spirit::x3;

	size_t size;
	const auto size_rule = long_;
	bool r = parse(first, last, size_rule >> eol, size);

	if (size == 0) {
		errx(1, "%s: matrix of size 0", file_name.c_str());
	}

	// prevent overflow
	static const auto MAX_SIZE = (size_t(1) << ((sizeof(size_t) >> 1) * 4)) - 1;
	if (size > MAX_SIZE) {
		errx(EINVAL, "%s: given matrix size is too big", file_name.c_str());
	}

	auto names = std::vector<std::string>{};
	auto values = std::vector<double>(size * size);
	names.reserve(size);

	/* The first line is special. We can use it to determine whether the input
	 * is in lower-triangle or full format. */

	auto lower_triangle = false;  // true iff given format is lower triangular
	auto diagonal_values = false; // true iff diagonal values are included
	auto first_line = parse_line_spirit(file_name, first, last, size);
	if (first_line.second.size() < size) {
		// assume lower triangle
		lower_triangle = true;
		if (first_line.second.size() == 1) {
			diagonal_values = true;
		}
	}

	names.push_back(first_line.first);
	std::copy(first_line.second.begin(), first_line.second.end(),
			  values.begin());

	// parse the rest of the matrix
	for (size_t i = 1; i < size; i++) {
		// read a different number of values, depending on the format.
		auto line_length = lower_triangle ? i + size_t(diagonal_values) : size;
		auto line = parse_line_spirit(file_name, first, last, line_length);

		names.push_back(line.first);

		std::copy(line.second.begin(), line.second.end(),
				  values.begin() + i * size);
	}

	auto ret = matrix{names, values};

	if (lower_triangle) {
		// fix upper triangle
		for (size_t i = 0; i < size; i++) {
			for (size_t j = i + 1; j < size; j++) {
				ret.entry(i, j) = ret.entry(j, i);
			}
		}
	}

	return ret;
}

/** @brief Parse the first matrix from a file and write it to a structure.
 *
 * @param file_name - The file to read.
 * @param out - An output iterator we write the matrices to.
 * @returns the position one past the last written matrix.
 */
template <typename OutputIt>
OutputIt parse_tolerant(const std::string &file_name, OutputIt out)
{
	std::ifstream file;
	std::istream *input = &std::cin;

	if (file_name != "-") {
		file = std::ifstream{file_name};
		input = &file;
	}

	input->unsetf(std::ios::skipws);
	if (!input || !*input) {
		err(errno, "%s", file_name.c_str());
	}

	boost::spirit::istream_iterator first(*input), last;

	if (input->good() && !input->eof()) {
		*out++ = parse_tolerant_internal(file_name, first, last);
	}

	return out;
}

/** @brief Parse all given file names into many matrices.
 *
 * @param argv - argv
 * @returns a list of matrices
 */
std::vector<matrix> parse_all(const char *const *argv)
{
	auto file_names = std::vector<std::string>{};
	while (*argv != nullptr) {
		file_names.push_back(*argv++);
	}

	return parse_all(file_names);
}

/** @brief Parse all given file names into many matrices.
 *
 * @param argv - argv
 * @returns a list of matrices
 */
std::vector<matrix> parse_all(std::vector<std::string> file_names)
{
	if (file_names.empty()) {
		// warn on reading from stdin
		if (isatty(STDIN_FILENO) != 0) {
			// Tell user we are expecting input …
			warnx("Reading from stdin…");
		}
		file_names.push_back("-");
	}

	auto matrices = std::vector<matrix>();
	auto inserter = std::back_inserter(matrices);
	matrices.reserve(file_names.size());

	for (const auto &file_name : file_names) {
		parse_tolerant(file_name, inserter);
	}

	return matrices;
}
