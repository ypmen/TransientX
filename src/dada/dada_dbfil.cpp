/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2025-09-22 18:16:25
 * @modify date 2025-09-22 18:16:25
 * @desc [description]
 */

#include <iostream>
#include <boost/program_options.hpp>

#include "dadareader.h"
#include "filterbank.h"
#include "filterbankwriter.h"
#include "logging.h"

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
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("key_input,k", value<std::string>()->default_value("2222"), "Input dada key")
			("config,c", value<std::string>(), "Config json file");

	positional_options_description pos_desc;
	pos_desc.add("config", -1);
	command_line_parser parser{argc, argv};
	parser.options(desc).style(command_line_style::default_style | command_line_style::allow_short);
	parser.options(desc).positional(pos_desc);
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

	// read config
	std::ifstream config_f(vm["config"].as<std::string>());
	nlohmann::json config = nlohmann::json::parse(config_f);
	std::string rootname = config["rootname"];
	float outmean = config["mean"];
	float outstd = config["std"];
	int outnbits = config["nbits"];

	// read dada header
	PSRDADA::DADAreader reader(vm["key_input"].as<std::string>());
	reader.read_header();

	long int nchans = reader.nchans;
	double tsamp = reader.tsamp;
	long int ndump = reader.get_bufsz() / (reader.nbits / 8) / nchans;
	
	DataBuffer<float> databuf(ndump, nchans);
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], reader.frequencies.data(), sizeof(double)*nchans);

	// setup filterbank writer
	FilterbankWriter filwriter;
	filwriter.fil.filename = rootname + ".fil";
	std::strcpy(filwriter.fil.source_name, reader.source_name.c_str());
	filwriter.fil.telescope_id = get_telescope_id(reader.telescope);
	std::string ra = reader.ra;
	std::string dec = reader.dec;
	ra.erase(remove(ra.begin(), ra.end(), ':'), ra.end());
	dec.erase(remove(dec.begin(), dec.end(), ':'), dec.end());
	filwriter.fil.src_raj = std::stod(ra);
	filwriter.fil.src_dej = std::stod(dec);
	filwriter.fil.ibeam = std::stoi(reader.beam);
	filwriter.fil.nbits = outnbits;
	filwriter.outmean = outmean;
	filwriter.outstd = outstd;
	filwriter.prepare(databuf);

	while (true)
	{
		if (reader.read_data(databuf, ndump) != ndump) 
		{
			BOOST_LOG_TRIVIAL(error)<<"data read failed: size not correct";
			break;
		}
		
		filwriter.run(databuf);
	}
}