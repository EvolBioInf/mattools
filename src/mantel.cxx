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
#include <iterator>
#include <random>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "matrix.h"

template<class T = std::mt19937, std::size_t N = T::state_size>
auto ProperlySeededRandomEngine () -> typename std::enable_if<!!N, T>::type {
    typename T::result_type random_data[N];
    std::random_device source;
    std::generate(std::begin(random_data), std::end(random_data), std::ref(source));
    std::seed_seq seeds(std::begin(random_data), std::end(random_data));
    T seededEngine (seeds);
    return seededEngine;
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
	auto begin = begin_square(ret);
	auto end = end_square(ret);

	std::transform(begin, end, begin,
				   [=](double value) { return (value - avg) / sd; });

	return ret;
}

double Z(const matrix &self, const matrix &other)
{
	auto names = common_names(self.get_names(), other.get_names());

	double dist = 0;

	for (size_t i = 0; i < names.size() - 1; i++) {
		const auto &name1 = names[i];
		for (size_t j = i + 1; j < names.size(); j++) {
			const auto &name2 = names[j];
			auto d1 = self.entry(name1, name2);
			auto d2 = other.entry(name1, name2);
			dist += d1 * d2;
		}
	}

	return dist;
}

double rmsd(const matrix &self, const matrix &other)
{
	const auto& names = self.get_names();

	double dist = 0;

	for (size_t i = 0; i < names.size() - 1; i++) {
		const auto &name1 = names[i];
		for (size_t j = i + 1; j < names.size(); j++) {
			const auto &name2 = names[j];
			auto d1 = self.entry(name1, name2);
			auto d2 = other.entry(name1, name2);
			auto t = d1 - d2;
			dist += t * t;
		}
	}

	auto count = names.size() * (names.size() - 1) / 2.0;
	return std::sqrt(dist / count);
}

double mantel(const matrix &self, const matrix &other, bool donormalize)
{
	using std::begin;
	using std::end;
	auto new_names = common_names(self.get_names(), other.get_names());

	auto size = new_names.size();
	auto new_self = sample2(self, new_names.begin(), new_names.end());
	auto new_other = sample2(other, new_names.begin(), new_names.end());

	if (donormalize) {
		new_self = normalize(new_self);
		new_other = normalize(new_other);
	}

	auto orig = rmsd(new_self, new_other);
	std::cout << "orig: " << orig << std::endl;

	auto montecarlo = std::vector<double>();
	auto indices = std::vector<size_t>(size);
	std::iota(begin(indices), end(indices), 0);

	auto g = ProperlySeededRandomEngine();

	auto count = size * (size - 1) / 2.0;

	for (size_t runs = 0; runs < 100000; runs++) {
		std::shuffle(begin(indices), end(indices), g);
		// randomize
		double dist = 0;

		for (size_t i = 0; i < size - 1; i++) {
			for (size_t j = i + 1; j < size; j++) {
				auto d1 = new_self.entry(i, j);
				auto d2 = new_other.entry(indices[i], indices[j]);
				auto t = d1 - d2;
				dist += t * t;
			}
		}

		std::cerr << "Z: " << std::sqrt(dist / count) << std::endl;
		montecarlo.push_back(std::sqrt(dist / count));
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
