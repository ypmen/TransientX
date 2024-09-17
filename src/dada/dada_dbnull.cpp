/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2024-04-10 12:36:37
 * @modify date 2024-04-10 12:36:37
 * @desc [description]
 */

#include <iostream>
#include <boost/program_options.hpp>

#include "dada.h"

using namespace boost::program_options;

unsigned int num_threads;

int main(int argc, const char *argv[])
{
	/* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("key_input,k", value<std::string>()->default_value("2222"), "Input dada key");

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

	nlohmann::json input_header;

	PSRDADA::Reader reader(vm["key_input"].as<std::string>());
	reader.prepare(input_header);

	long int bufsz = reader.get_bufsz();

	std::vector<char> tim(bufsz, 0);

	while (true)
	{
		reader.run((char *)(tim.data()), tim.size());
	}
}