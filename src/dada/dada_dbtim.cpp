/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2024-04-05 09:34:29
 * @modify date 2024-04-05 09:34:29
 * @desc "dump psrdada buffer to disk"
 */

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>

#include "dada.h"
#include "logging.h"

using namespace boost::program_options;

unsigned int num_threads;

int main(int argc, const char *argv[])
{
	init_logging();

	/* options */
	int verbose = 0;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("rootname,o", value<std::string>()->default_value("fake"), "Rootname of output tim files")
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

	num_threads = vm["threads"].as<unsigned int>();
	std::string rootname = vm["rootname"].as<std::string>();

	nlohmann::json input_header;

	PSRDADA::Reader reader(vm["key_input"].as<std::string>());
	reader.prepare(input_header);

	double tsamp = input_header["tsamp"];
	double dms = input_header["dms"];
	double ddm = input_header["ddm"];
	double ndm = input_header["ndm"];

	long int ndump = reader.get_bufsz() / sizeof(float) / ndm;

	std::vector<std::ofstream> outfiles;

	for (long int k=0; k<ndm; k++)
	{
		double dm = dms + k * ddm;
			
		std::stringstream ss_dm;
		ss_dm << "DM" << std::setprecision(2) << std::fixed << std::setfill('0') << dm;
		std::string s_dm = ss_dm.str();

		std::string fname = rootname + "_" + s_dm + ".dat";
		std::ofstream outfile;
		outfile.open(fname, std::ios::binary);
		outfiles.push_back(std::move(outfile));
	}

	std::vector<float> tim(ndm * ndump);

	while (true)
	{
		reader.run((char *)(tim.data()), tim.size() * sizeof(float));
		for (size_t k=0; k<ndm; k++)
		{
			outfiles[k].write((char *)(tim.data() + k * ndump), ndump * sizeof(float));
		}
	}
}