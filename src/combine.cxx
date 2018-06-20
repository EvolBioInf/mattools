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

matrix parse_tolerant_internal(const std::string &file_name,
							   std::istream *input);

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

matrix combine(const matrix &self, const matrix &other)
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

			if (self.has_coverages() && other.has_coverages()) {
				auto c1 = self.cov_entry(this_map[name1], this_map[name2]);
				auto c2 = other.cov_entry(other_map[name1], other_map[name2]);

				auto val = c1 > c2 ? d1 : d2;
				// auto val = (d1 * c1 + d2 * c2 ) / (c1 + c2);

				ret.entry(i, j) = ret.entry(j, i) = val;
			} else {
				ret.entry(i, j) = ret.entry(j, i) = std::max(d1, d2);
			}
		}
	}

	return ret;
}

static void mat_combine_usage(int status);

auto parse_line_without_name(const std::string &file_name, std::istream *input,
							 size_t max_values)
{
	auto values = std::vector<double>{};

	while (max_values-- > 0) {
		auto str = std::string();
		auto value = 0.0;

		*input >> str;
		try {
			value = std::stod(str); // also parses nan
		} catch (std::exception &e) {
			/* If error is recoverable, reset. This happens when we read past
			 * the end of a line and tried to interpret the next name as a value
			 * instead.
			 */
			if (input->fail()) {
				input->clear();
			}

			input->putback(' ');
			// not a double, probably a name. push back
			for (auto it = str.crbegin(); it != str.crend(); it++) {
				input->putback(*it);
			}

			input->putback(' ');
			break; // out of the loop
		}

		values.push_back(value);
	}

	return values;
}

template <typename OutputIt>
OutputIt parse_tolerant_with_coverage(const std::string &file_name,
									  OutputIt out)
{
	std::ifstream file;
	std::istream *input = &std::cin;

	if (file_name != "-") {
		file = std::ifstream{file_name};
		input = &file;
	}

	if (!input || !*input) {
		err(errno, "%s", file_name.c_str());
	}

	const auto skip_blank_lines = [](std::istream *input) {
		while (input->peek() == '\n') {
			input->get();
		}
	};

	if (skip_blank_lines(input), (input->good() && !input->eof())) {
		auto distance_matrix = parse_tolerant_internal(file_name, input);

		auto size = distance_matrix.get_size();
		auto coverages = std::vector<double>(size * size);
		skip_blank_lines(input);

		auto dummy = std::string();
		*input >> dummy; // "Coverages:"

		for (size_t i = 0; i < size; i++) {
			auto values = parse_line_without_name(file_name, input, size);

			std::copy(values.begin(), values.end(),
					  coverages.begin() + (i * size));
		}

		*out++ = matrix(distance_matrix.get_names(),
						distance_matrix.get_values(), coverages);
	}

	return out;
}

/**
 * @brief The main function of `mat compare`.
 *
 * @param argc - It's argc, stupid. (Advanced by one from global argc.)
 * @param argv - It's argv, stupid. (Advanced by one from global argv.)
 * @returns 0 iff successful.
 */
int mat_combine(int argc, char **argv)
{
	static struct option long_options[] = {
		{"help", no_argument, 0, 0}, // print help
		{0, 0, 0, 0}				 //
	};

	bool full_matrix = false;

	while (true) {
		int long_index;
		int c = getopt_long(argc, argv, "f", long_options, &long_index);

		if (c == -1) {
			break;
		}

		if (c == 0 && std::string(long_options[long_index].name) == "help") {
			mat_combine_usage(EXIT_SUCCESS);
		} else {
			mat_combine_usage(EXIT_FAILURE);
		}
	}

	argc -= optind, argv += optind; // hack

	// auto matrices = parse_all(argv);
	auto file_names = std::vector<std::string>(argv, argv + argc);

	// check for STDIN

	auto matrices = std::vector<matrix>();
	auto inserter = std::back_inserter(matrices);
	matrices.reserve(file_names.size());

	for (const auto &file_name : file_names) {
		parse_tolerant_with_coverage(file_name, inserter);
	}

	if (matrices.size() < 2) {
		errx(1, "At least two matrices must be provided.");
	}

	std::cout << combine(matrices[0], matrices[1]).to_string() << std::endl;

	return 0;
}

static void mat_combine_usage(int status)
{
	static const char str[] = {
		"usage: mat combine [OPTIONS] [FILE...]\n" // this comment is a hack
		"Combine two distances matrices.\n\n"
		"Available options:\n"
		"     --help          print this help\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}
