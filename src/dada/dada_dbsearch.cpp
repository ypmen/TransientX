/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2024-10-06 21:06:17
 * @modify date 2024-10-06 21:06:17
 * @desc [description]
 */

#include <iostream>
#include <boost/program_options.hpp>

#include "boxcar.h"
#include "cluster.h"
#include "candplot.h"

#include "dada.h"
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
			("key_input", value<std::string>()->default_value("2222"), "Input dada tim key")
			("key_input2", value<std::string>()->default_value("3333"), "Input dada sub key")
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

	nlohmann::json config_boxcar = config["boxcar"];
	nlohmann::json config_clustering = config["clustering"];
	nlohmann::json config_candplot = config["candplot"];

	num_threads = vm["threads"].as<unsigned int>();

	nlohmann::json dedisp_header;
	nlohmann::json dedisp_header2;

	PSRDADA::Reader reader(vm["key_input"].as<std::string>());
	reader.prepare(dedisp_header);

	PSRDADA::Reader reader_sub(vm["key_input2"].as<std::string>());
	reader_sub.prepare(dedisp_header2);

	// reconstruct dedisp
	RealTime::SubbandDedispersion dedisp;
	dedisp.dms = dedisp_header["dms"];
	dedisp.ddm = dedisp_header["ddm"];
	dedisp.ndm = dedisp_header["ndm"];
	dedisp.ndump = reader.get_bufsz()/ sizeof(float) / dedisp.ndm;

	nlohmann::json obsinfo_header = dedisp_header["obsinfo"];

	DataBuffer<float> databuf(dedisp.ndump, obsinfo_header["nchans"]);
	databuf.tsamp = obsinfo_header["tsamp"];
	double fch1 = obsinfo_header["fch1"];
	double foff = obsinfo_header["foff"];
	for (size_t j=0; j<obsinfo_header["nchans"]; j++)
	{
		databuf.frequencies[j] = fch1 + j * foff;
	}

	dedisp.prepare(databuf);

	// configure boxcar
	Boxcar boxcar(config_boxcar);
	boxcar.prepare(dedisp);

	// configure clustering
	Cluster<double> cluster(config_clustering);

	// configure candplot
	CandPlot candplot(config_candplot);

	/** form obsinfo*/
	std::map<std::string, std::string> obsinfo;
	obsinfo["Source_name"] = obsinfo_header["source_name"];
	obsinfo["Telescope"] = obsinfo_header["telescope"];
	stringstream ss_tstart;
	ss_tstart << setprecision(13) << fixed << obsinfo_header["tstart"];
	string s_tstart = ss_tstart.str();
	obsinfo["Tstart"] = s_tstart;
	obsinfo["RA"] = obsinfo_header["ra"];
	obsinfo["DEC"] = obsinfo_header["dec"];
	stringstream ss_beam;
	ss_beam << setw(5) << setfill('0') << std::stoi(obsinfo_header["beam"].get<std::string>());
	obsinfo["Beam"] = "Beam" + ss_beam.str();

	double gl = 0., gb = 0.;
#ifdef HAVE_SOFA
	get_gl_gb(gl, gb, obsinfo["RA"], obsinfo["DEC"]);
#endif
	obsinfo["GL"] = to_string(gl);
	obsinfo["GB"] = to_string(gb);

	double ymw16_maxdm = 0.;
	ymw16_maxdm = get_maxdm_ymw16(gl, gb);

	obsinfo["MaxDM_YMW16"] = to_string(ymw16_maxdm);

	// run pipeline
	std::string fname = "";

	while (true)
	{
		dedisp.counter += dedisp.ndump;
		if (dedisp.counter < dedisp.offset-dedisp.noverlap+dedisp.ndump) continue;

		reader.run((char *)(dedisp.sub.buffertim.data()), dedisp.ndm * dedisp.ndump * sizeof(float));
		reader_sub.run((char *)(dedisp.sub.bufferT.data()), dedisp.sub.nsub * dedisp.sub.nchans * dedisp.sub.nsamples * sizeof(float));

		if (boxcar.run(dedisp))
		{
			if (cluster.run(boxcar))
			{
				candplot.plot(cluster, boxcar, dedisp, 0, fname, obsinfo);
			}
		}
	}
}