/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2024-04-10 12:05:53
 * @modify date 2024-04-10 12:05:53
 * @desc "throughput benchmark"
 */

#include <iostream>
#include <chrono>
#include <boost/program_options.hpp>

#include "dada.h"

using namespace boost::program_options;

int main(int argc, const char *argv[])
{
	/* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("num,n", value<int>()->default_value(100), "Make statistics for every n blocks")
			("key_output,k", value<std::string>()->default_value("1111"), "Output dada key");

	command_line_parser parser{argc, argv};
	parser.options(desc).style(command_line_style::default_style | command_line_style::allow_short);
	parsed_options parsed_options = parser.run();

	variables_map vm;
	store(parsed_options, vm);
	notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << '\n';
		return 0;
	}
	if (vm.count("verbose"))
	{
		verbose = 1;
	}

	// write output header
	nlohmann::json output_header = {
		{"telescope", "benchmark"},
		{"source_name", "benchmark"},
		{"ra", "00:00:00"},
		{"dec", "00:00:00"},
		{"beam", "0"},
		{"tstart", "0"},
		{"nifs", 1},
		{"nbits", 8},
		{"nchans", 2048},
		{"tsamp", 100e-6},
		{"fch1", 1600.},
		{"foff", -0.4}
	};

	std::string key_output = vm["key_output"].as<std::string>();
	PSRDADA::Writer writer(key_output);
	writer.prepare(output_header);

	// create databuffer
	long int bufsz = writer.get_bufsz();
	std::vector<unsigned char> buf(bufsz, 0);

	int n = vm["num"].as<int>();

	while (true)
	{
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

		for (size_t k=0; k<n; k++)
		{
			writer.run((char *)(buf.data()), buf.size());
		}

		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

		double time_elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();

		std::cout << "throughput = " << n * bufsz / 1024 / 1024 / time_elapsed << " MB/s" << std::endl;
	}
}