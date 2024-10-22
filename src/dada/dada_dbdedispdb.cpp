/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2024-04-04 13:51:29
 * @modify date 2024-04-04 13:51:29
 * @desc [description]
 */

#include <iostream>
#include <boost/program_options.hpp>

#include "preprocesslite.h"
#include "pipeline.h"
#include "subdedispersion.h"
#include "dada.h"
#include "dadareader.h"
#include "logging.h"
#include "json.hpp"

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
			("key_input,i", value<std::string>()->default_value("1111"), "Input dada key")
			("key_output,o", value<std::string>()->default_value("2222"), "Output dada key")
			("key_output2", value<std::string>(), "Output subband dada key")
			("dada_args", value<std::string>()->default_value(""), "Extra dada args")
			("nblock,n", value<int>()->default_value(16), "Number of output dada block")
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
	if (vm.count("config") == 0)
	{
		cerr<<"Error: no config file"<<endl;
		return -1;
	}

	std::ifstream config_f(vm["config"].as<std::string>());
	nlohmann::json config = nlohmann::json::parse(config_f);

	nlohmann::json config_prep = config["preprocesslite"];
	nlohmann::json config_ds = config["downsample"];
	nlohmann::json config_eq;
	nlohmann::json config_bs = config["baseline"];
	nlohmann::json config_rfi = config["rfi"];
	nlohmann::json config_dedisp = config["subdedispersion"];

	num_threads = vm["threads"].as<unsigned int>();

	PSRDADA::DADAreader reader(vm["key_input"].as<std::string>());
	reader.read_header();

	long int nchans = reader.nchans;
	double tsamp = reader.tsamp;

	long int ndump = reader.get_bufsz() / (reader.nbits / 8) / nchans;

	DataBuffer<float> databuf(ndump, nchans);
	databuf.tsamp = tsamp;
	memcpy(&databuf.frequencies[0], reader.frequencies.data(), sizeof(double)*nchans);

	PreprocessLite prep(config_prep);
	prep.prepare(databuf);

	XLIBS::Pipeline pipeline(config_ds, config_eq, config_bs, config_rfi);
	pipeline.prepare(prep);

	RealTime::SubbandDedispersion dedisp(config_dedisp);
	dedisp.ndump = pipeline.nsamples;
	dedisp.prepare(pipeline);

	// create output dada buffer
	std::string cmd = "dada_db " + vm["dada_args"].as<std::string>() + " ";
	cmd += "-b " + std::to_string(dedisp.ndm * dedisp.ndump * sizeof(float)) + " ";
	cmd += "-n " + std::to_string(vm["nblock"].as<int>()) + " ";
	cmd += "-k " + vm["key_output"].as<std::string>();	
	system(cmd.c_str());

	// write output header
	nlohmann::json output_header = {
		{"dms", dedisp.dms},
		{"ddm", dedisp.ddm},
		{"ndm", dedisp.ndm},
		{"obsinfo", {
			{"telescope", reader.telescope},
			{"source_name", reader.source_name},
			{"ra", reader.ra},
			{"dec", reader.dec},
			{"beam", reader.beam},
			{"tstart", reader.start_mjd.to_day()},
			{"tsamp", reader.tsamp},
			{"fch1", reader.frequencies.front()},
			{"foff", reader.frequencies[1] - reader.frequencies[0]},
			{"nchans", reader.nchans}
		}}
	};

	PSRDADA::Writer writer(vm["key_output"].as<std::string>());
	writer.prepare(output_header);

	PSRDADA::Writer writer_subband;
	if (vm.count("key_output2"))
	{
		std::string cmd = "dada_db ";
		cmd += "-b " + std::to_string(dedisp.sub.nsub * dedisp.sub.nchans * dedisp.sub.nsamples * sizeof(float)) + " ";
		cmd += "-n " + std::to_string(vm["nblock"].as<int>()) + " ";
		cmd += "-k " + vm["key_output2"].as<std::string>();	
		system(cmd.c_str());

		writer_subband.setup(vm["key_output2"].as<std::string>());
		writer_subband.prepare(output_header);
	}
	
	while (true)
	{
		if (reader.read_data(databuf, ndump) != ndump) 
		{
			BOOST_LOG_TRIVIAL(error)<<"data read failed: size not correct";
			break;
		}

		databuf.counter += ndump;

		DataBuffer<float> *data = prep.run(databuf);

		data = pipeline.run(*data);

		dedisp.run(*data, data->nsamples);

		if (dedisp.counter < dedisp.offset-dedisp.noverlap+dedisp.ndump) continue;

		writer.run((char *)dedisp.sub.buffertim.data(), dedisp.ndm * dedisp.ndump * sizeof(float));

		if (vm.count("key_output2"))
		{
			writer_subband.run((char *)dedisp.sub.bufferT.data(), dedisp.sub.nsub * dedisp.sub.nchans * dedisp.sub.nsamples * sizeof(float));
		}
	}
}