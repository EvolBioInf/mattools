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
#include <err.h>
#include <stdio.h>
#include <string>

int mat_compare(int, char **);
int mat_grep(int, char **);
int mat_nj(int, char **);
int mat_format(int, char **);
static void usage(int status);
static void version();

/** @brief The main function of the mat tools. It dispatches into the sub
 * commands.
 *
 * @param argc - argc.
 * @param argv - argv.
 * @returns 0 iff successful.
 */
int main(int argc, char *argv[])
{
	if (argc < 2) {
		usage(EXIT_FAILURE);
	}

	auto first_arg = std::string{argv[1]};
	if (first_arg == "--version") {
		version();
		exit(EXIT_SUCCESS);
	}

	if (first_arg == "--help") {
		usage(EXIT_SUCCESS);
	}

	// strip binary
	argc -= 1, argv += 1;

	auto command = first_arg;
	if (command == "compare") {
		return mat_compare(argc, argv);
	}

	if (command == "grep") {
		return mat_grep(argc, argv);
	}

	if (command == "nj") {
		return mat_nj(argc, argv);
	}

	if (command == "format") {
		return mat_format(argc, argv);
	}

	warnx("unknown command '%s'.", command.c_str());
	usage(EXIT_FAILURE);

	return 1;
}

static void usage(int status)
{
	static const char str[] = {
		"usage: mat [--version] [--help] <command> [<args>]\n\n"
		"The available commands are:\n"
		" compare     Compute the distance between two matrices\n"
		" format      Format the distance matrix\n"
		" grep        Print submatrix for names matching a pattern\n"
		" nj          Convert to a tree by neighbor joining\n"
		"\n"
		"Use 'mat <command> --help' to get guidance on the usage of a "
		"command.\n"};

	fprintf(status == EXIT_SUCCESS ? stdout : stderr, str);
	exit(status);
}

static void version()
{
	static const char str[] = {PACKAGE_STRING "\n"};
	printf(str);
}
