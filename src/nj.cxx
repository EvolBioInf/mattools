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

/*
 * Note, this source is adapted from afra.
 * Original Source: https://github.com/EvolBioInf/afra
 * Copyright (C) 2015 - 2016 Fabian Klötzl, GPLv3+
 */

#include <cassert>
#include <getopt.h>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>
#include "matrix.h"

static size_t sample_size = 0;
using support_fn_type = double(const matrix &, const std::vector<uint8_t> &);
static support_fn_type *support_fn = nullptr;
static unsigned long seed = 0;

class tree_node
{
  public:
	tree_node *left_branch = nullptr, *right_branch = nullptr;
	double left_dist = 0.0, right_dist = 0.0;
	double left_support = 0.0, right_support = 0.0;
	ssize_t index = 0;

	tree_node() = default;
	explicit tree_node(ssize_t _index) noexcept : index(_index)
	{
	}

	tree_node(tree_node *lb, tree_node *rb, double ld, double rd) noexcept
		: left_branch{lb}, right_branch{rb}, left_dist{ld},
		  right_dist{rd}, index{-1}
	{
	}

	template <typename Func>
	void traverse(const Func &process)
	{
		if (left_branch) {
			left_branch->traverse(process);
		}
		process(this);
		if (right_branch) {
			right_branch->traverse(process);
		}
	}

	template <typename Func1, typename Func2, typename Func3>
	void traverse(Func1 *pre = nullptr, Func2 *process = nullptr,
				  Func3 *post = nullptr)
	{
		if (pre) {
			(*pre)(this);
		}
		if (left_branch) {
			left_branch->traverse(pre, process, post);
		}
		if (process) {
			(*process)(this);
		}
		if (right_branch) {
			right_branch->traverse(pre, process, post);
		}
		if (post) {
			(*post)(this);
		}
	}
};

class tree_root : public tree_node
{
  public:
	tree_node *extra_branch = nullptr;
	double extra_dist = 0.0;
	double extra_support = 0.0;

	tree_root() = default;
	tree_root(tree_node *lb, tree_node *rb, tree_node *eb, double ld, double rd,
			  double ed) noexcept
		: tree_node(lb, rb, ld, rd), extra_branch{eb}, extra_dist{ed}
	{
	}
};

class tree
{
  public:
	size_t size{0};
	std::vector<tree_node> pool{};
	tree_root root{};

	tree() = default;
	explicit tree(size_t s) : size{s}, pool(2 * size)
	{
	}
};

void colorize(tree_node *self, std::vector<uint8_t> &buffer, uint8_t color);
double support_full(const matrix &distance, const std::vector<uint8_t> &buffer);

tree nj(const matrix &m)
{
	auto ret = tree{m.get_size()};

	auto matrix_size = m.get_size();

	auto node_pool = ret.pool.data();
	auto empty_node_ptr = node_pool + matrix_size;
	auto unjoined_nodes = std::vector<tree_node *>{};
	unjoined_nodes.reserve(matrix_size);

	for (size_t i = 0; i < matrix_size; i++) {
		node_pool[i] = tree_node{static_cast<ssize_t>(i)}; // leaf
		unjoined_nodes[i] = &node_pool[i];
	}

	auto r = std::vector<double>(matrix_size);
	auto local_copy = matrix{m};

	auto n = matrix_size;
	while (n > 3) {
		for (size_t i = 0; i < n; i++) {
			auto row = local_copy.row(i);
			// row_end() != row() + n
			auto rr = std::accumulate(row, row + n, 0.0);
			r[i] = rr / (n - 2);
		}

		size_t min_i = 0, min_j = 1;
		double min_value = local_copy.entry(0, 1) - r[0] - r[1];

		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				if (i == j) continue;

				double value = local_copy.entry(i, j) - r[i] - r[j];
				if (value < min_value) {
					min_i = i;
					min_j = j;
					min_value = value;
				}
			}
		}

		// force i < j
		if (min_j < min_i) {
			std::swap(min_i, min_j);
		}

		auto branch = tree_node{
			unjoined_nodes[min_i], unjoined_nodes[min_j],
			(local_copy.entry(min_i, min_j) + r[min_i] - r[min_j]) / 2.0,
			(local_copy.entry(min_i, min_j) - r[min_i] + r[min_j]) / 2.0};

		*empty_node_ptr++ = branch;
		unjoined_nodes[min_i] = empty_node_ptr - 1;
		unjoined_nodes[min_j] = unjoined_nodes[n - 1];

		double row_k[n];
		double M_ij = local_copy.entry(min_i, min_j);

		for (size_t m = 0; m < n; m++) {
			if (m == min_i || m == min_j) continue;

			row_k[m] = (local_copy.entry(min_i, m) +
						local_copy.entry(min_j, m) - M_ij) /
					   2.0;
			// if( row_k[m] < 0) row_k[m] = 0;
		}

		// row_k[min_i] and row_k[min_j] are undefined!
		row_k[min_i] = 0.0;
		row_k[min_j] = row_k[n - 1];

#define M(i, j) local_copy.entry(i, j)

		memmove(&M(min_i, 0), row_k, n * sizeof(double));
		memmove(&M(min_j, 0), &M((n - 1), 0), n * sizeof(double));

		// zero main diagonal
		M(min_i, min_i) = M(min_j, min_j) = 0.0;

		// restore symmetry
		for (size_t i = 0; i < n; i++) {
			M(i, min_i) = M(min_i, i);
		}

		for (size_t i = 0; i < n; i++) {
			M(i, min_j) = M(min_j, i);
		}

		n--;
	}

	// join three remaining nodes
	auto root = tree_root{unjoined_nodes[0],
						  unjoined_nodes[1],
						  unjoined_nodes[2],
						  (M(0, 1) + M(0, 2) - M(1, 2)) / 2.0,
						  (M(0, 1) + M(1, 2) - M(0, 2)) / 2.0,
						  (M(0, 2) + M(1, 2) - M(0, 1)) / 2.0};

	//*empty_node_ptr++ = root;
	ret.root = root;

	return ret;
}

std::string to_newick(const tree &t, const matrix &m)
{
	auto ret = std::string{};
	auto root = &t.root;

	auto pre = [&ret](const tree_node *self) {
		if (self->left_branch) {
			ret += "(";
		}
	};
	auto process = [&ret, &m](const tree_node *self) {
		if (self->left_branch) {
			if (self->left_branch->left_branch) {
				ret += std::to_string((int)(self->left_support * 100));
			}
			ret += ":" + std::to_string(self->left_dist) + ",";
		} else {
			ret += m.get_names()[self->index];
		}
	};
	auto post = [&ret](const tree_node *self) {
		if (!self->right_branch) return;
		if (self->right_branch->right_branch) {
			ret += std::to_string((int)(self->right_support * 100));
		}

		char buf[20];
		snprintf(buf, sizeof(buf), "%1.4e", self->right_dist);
		ret += std::string(":") + buf + ")";
	};

	ret += "(";
	root->left_branch->traverse(&pre, &process, &post);
	process(root);

	root->right_branch->traverse(&pre, &process, &post);
	if (root->right_branch) {
		if (root->right_branch->right_branch) {
			ret += std::to_string((int)(root->right_support * 100));
		}

		char buf[20];
		snprintf(buf, sizeof(buf), "%1.4e", root->right_dist);
		ret += std::string(":") + buf + ",";
	}

	root->extra_branch->traverse(&pre, &process, &post);
	if (root->extra_branch) {
		if (root->extra_branch->left_branch) {
			ret += std::to_string((int)(root->extra_support * 100));
		}

		char buf[20];
		snprintf(buf, sizeof(buf), "%1.4e", root->extra_dist);
		ret += std::string(":") + buf;
	}
	ret += ");";

	return ret;
}

enum { SET_D = 0, SET_A, SET_B, SET_C };

/** Colorize according to the following scheme.
 *
 *  A -left--             -right- C
 *           \           /
 *            --left-- self
 *           /           \
 *  B -right-             -extra- D
 *
 * self is the node of the caller. It has two branches pointing down (left
 * and right). The extra branch points to the parent.
 *
 */

void quartet_left(tree_node *self, const matrix &distance)
{
	if (!self->left_branch || !self->left_branch->left_branch) return;

	auto buffer = std::vector<uint8_t>(distance.get_size(), SET_D);

	colorize(self->left_branch->left_branch, buffer, SET_A);
	colorize(self->left_branch->right_branch, buffer, SET_B);
	colorize(self->right_branch, buffer, SET_C);

	self->left_support = support_fn(distance, buffer);
}

void quartet_right(tree_node *self, const matrix &distance)
{
	if (!self->left_branch || !self->right_branch->left_branch) return;

	auto buffer = std::vector<uint8_t>(distance.get_size(), SET_D);

	colorize(self->right_branch->left_branch, buffer, SET_A);
	colorize(self->right_branch->right_branch, buffer, SET_B);
	colorize(self->left_branch, buffer, SET_C);

	self->right_support = support_fn(distance, buffer);
}

void colorize(tree_node *self, std::vector<uint8_t> &buffer, uint8_t color)
{
	if (!self) return;

	self->traverse([color = color, &buffer](const tree_node *self) {
		if (self->left_branch == nullptr) {
			buffer[self->index] = color;
		}
	});
}

void quartet_all(tree &baum, const matrix &distance)
{
	// iterate over all nodes
	size_t size = distance.get_size();
	tree_node *inner_nodes = baum.pool.data() + size;

	// #pragma omp parallel for schedule(dynamic) num_threads(THREADS)
	for (size_t i = 0; i < size - 2; i++) {
		quartet_left(&inner_nodes[i], distance);
		quartet_right(&inner_nodes[i], distance);
	}

	tree_root *root = &baum.root;
	quartet_left(root, distance);
	quartet_right(root, distance);

	if (root->extra_branch->left_branch) {
		// Support Value for Root→Extra
		auto buffer = std::vector<uint8_t>(distance.get_size(), SET_D);

		colorize(root->extra_branch->left_branch, buffer, SET_A);
		colorize(root->extra_branch->right_branch, buffer, SET_B);
		colorize(root->left_branch, buffer, SET_C);

		root->extra_support = support_fn(distance, buffer);
	}
}

struct quartet {
	size_t A, B, C, D;
	bool operator<(struct quartet other) const noexcept
	{
		return memcmp(this, &other, sizeof(*this)) < 0;
	}
};

double support_sample(const matrix &distance,
					  const std::vector<uint8_t> &buffer)
{
	const size_t size = distance.get_size();

	size_t set_sizes[4] = {0};
	for (size_t i = 0; i < size; i++) {
		set_sizes[buffer[i]]++;
	}

	size_t quartet_number = set_sizes[SET_A] * set_sizes[SET_B] *
							set_sizes[SET_C] * set_sizes[SET_D];

	if (quartet_number < sample_size) {
		return support_full(distance, buffer);
	}

	// sample `sample_size` quartets and use them to compute the support
	auto quartets = std::set<struct quartet>();

	if (seed == 0) {
		std::random_device rd;
		seed = rd();
	}

	// todo: seeding from a single int is insufficient.
	auto engine = std::default_random_engine{seed};
	std::uniform_int_distribution<size_t> dists[4];
	for (size_t i = 0; i < 4; i++) {
		dists[i] = std::uniform_int_distribution<size_t>{
			0, set_sizes[i] - 1}; // both values are inclusive
	}

	// convert number in set to index into matrix
	auto indices_A = std::vector<size_t>(set_sizes[SET_A]);
	auto indices_B = std::vector<size_t>(set_sizes[SET_B]);
	auto indices_C = std::vector<size_t>(set_sizes[SET_C]);
	auto indices_D = std::vector<size_t>(set_sizes[SET_D]);
	auto i_A = indices_A.begin();
	auto i_B = indices_B.begin();
	auto i_C = indices_C.begin();
	auto i_D = indices_D.begin();
	for (size_t i = 0; i < buffer.size(); i++) {
		switch (buffer[i]) {
			case SET_A: *i_A++ = i; break;
			case SET_B: *i_B++ = i; break;
			case SET_C: *i_C++ = i; break;
			case SET_D: *i_D++ = i; break;
		}
	}

	// generate quartets
	while (quartets.size() < sample_size) {
		struct quartet q = {
			indices_A[dists[SET_A](engine)], indices_B[dists[SET_B](engine)],
			indices_C[dists[SET_C](engine)], indices_D[dists[SET_D](engine)]};
		quartets.insert(q);
	}

	size_t non_supporting_counter = 0;
	auto q = quartets.begin();

	for (size_t i = 0; i < sample_size; i++, q++) {
		double D_abcd = distance.entry(q->A, q->B) + distance.entry(q->C, q->D);

		if (distance.entry(q->A, q->C) + distance.entry(q->B, q->D) < D_abcd ||
			distance.entry(q->A, q->D) + distance.entry(q->B, q->C) < D_abcd) {

			non_supporting_counter++;
		}
	}

	return 1 - (static_cast<double>(non_supporting_counter) / sample_size);
}

double support_full(const matrix &distance, const std::vector<uint8_t> &buffer)
{
	const size_t size = distance.get_size();

	size_t non_supporting_counter = 0;
	size_t quartet_counter = 0;

	size_t A = 0, B, C, D;
	for (; A < size; A++) {
		if (buffer[A] != SET_A) continue;

		for (B = 0; B < size; B++) {
			if (buffer[B] != SET_B) continue;

			for (C = 0; C < size; C++) {
				if (buffer[C] != SET_C) continue;

				for (D = 0; D < size; D++) {
					if (buffer[D] != SET_D) continue;

					quartet_counter++;

					double D_abcd = distance.entry(A, B) + distance.entry(C, D);
					if (distance.entry(A, C) + distance.entry(B, D) < D_abcd ||
						distance.entry(A, D) + distance.entry(B, C) < D_abcd) {

						// printf("%zu %zu %zu %zu\n", A, B, C, D);
						non_supporting_counter++;
					}
				}
			}
		}
	}

	return 1 - (static_cast<double>(non_supporting_counter) / quartet_counter);
}

static void mat_nj_usage(int);

/**
 * @brief The main function of `mat nj`.
 *
 * @param argc - It's argc, stupid. (Advanced by one from global argc.)
 * @param argv - It's argv, stupid. (Advanced by one from global argv.)
 * @returns 0 iff successful.
 */
int mat_nj(int argc, char **argv)
{
	static struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"sample-size", required_argument, 0, 0},
		{"seed", required_argument, 0, 0},
		{"no-support", no_argument, 0, 0},
		{0, 0, 0, 0} // hack
	};

	support_fn = &support_full;
	auto support = true;

	while (true) {
		int long_index;
		int c = getopt_long(argc, argv, "h", long_options, &long_index);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: {
				auto option_str = std::string(long_options[long_index].name);
				if (option_str == "no-support") {
					support = false;
					break;
				}
				if (option_str == "sample-size") {
					sample_size = std::stod(optarg);
					support = true;
					support_fn = &support_sample;
					break;
				}
				if (option_str == "seed") {
					seed = std::stod(optarg);
					break;
				}
				break;
			}
			case 'h': mat_nj_usage(EXIT_SUCCESS);
			case '?': // intentional fall-through
			default: mat_nj_usage(EXIT_FAILURE);
		}
	}

	argc -= optind, argv += optind;

	auto matrices = parse_all(argv);

	for (const auto &mat : matrices) {
		if (mat.get_size() < 4) {
			errx(1, "expected at least four species");
		}

		auto t = nj(mat);
		if (support) {
			quartet_all(t, mat);
		}
		std::cout << to_newick(t, mat) << std::endl;
	}

	return 0;
}

static void mat_nj_usage(int status)
{
	static const char str[] = {
		"usage: mat nj [OPTIONS] [FILE...]\n"
		"Build a tree via neighbor joining.\n\n"
		"Available options:\n"
		"  -h, --help           print this help\n"
		"      --no-support     do not compute support values\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}
