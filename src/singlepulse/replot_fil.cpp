/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-04-13 20:57:43
 * @modify date 2021-04-13 20:57:43
 * @desc [description]
 */

#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "logging.h"
#include "dedisperse.h"
#include "downsample.h"
#include "candidate.h"
#include "psrfits.h"
#include "psrfitsreader.h"
#include "filterbankreader.h"
#include "ringbuffer.h"

#include "json.hpp"
using json = nlohmann::json;

using namespace boost::program_options;

unsigned int num_threads;

template<typename InputIterator, typename ValueType>
ValueType closestdist(InputIterator first, InputIterator last, ValueType value)
{
	ValueType temp = *std::min_element(first, last, [&](ValueType x, ValueType y){return std::abs(x - value) < std::abs(y - value);});
	return std::abs(temp-value);
}

int main(int argc, char *argv[])
{
	init_logging();

	/* options */
	int verbose = 0;
	bool outbest = false;

	options_description desc{"Options"};
	desc.add_options()
			("help,h", "Help")
			("verbose,v", "Print debug information")
			("threads,t", value<unsigned int>()->default_value(1), "Number of threads")
			("td", value<int>()->default_value(1), "Time downsample")
			("fd", value<int>()->default_value(1), "Frequency downsample")
			("nw", value<int>()->default_value(20), "Chop data by nw*width")
			("snrloss", value<float>()->default_value(1e-1), "S/N loss")
			("dm", value<double>(), "Update dm")
			("coherent", "Apply coherent dedispersion")
			("zapthre", value<float>()->default_value(3), "Threshold in IQR for zapping channels")
			("zap", value<std::vector<double>>()->multitoken()->zero_tokens()->composing(), "Zap channels, e.g. --zap 1000 1100 1200 1300")
			("zdot", "Perform zero-DM matched filter")
			("clip", value<std::vector<int>>()->multitoken(), "Perform clip filter [td, fd, threshold]")
			("kadane", value<std::vector<int>>()->multitoken(), "Perform kadane filter [td, fd, threshold]")
			("candfile", value<std::string>(), "Input cand file")
			("template", value<std::string>(), "Input fold template file")
			("srcname", value<std::string>(), "Souce name")
			("ibeam,i", value<int>(), "Beam number")
			("telescope", value<std::string>(), "Telescope name")
			("ra", value<std::string>(), "RA (hhmmss.s)")
			("dec", value<std::string>(), "DEC (ddmmss.s)")
			("incoherent", "The beam is incoherent (ifbf). Coherent beam by default (cfbf)")
			("arch,a", "Generate archive file")
			("pow2bin", "choose nbin of power 2")
			("nopt", "Don't perform optimization")
			("dmcutoff", value<float>()->default_value(0), "DM cutoff")
			("ddmcutoff", value<float>()->default_value(5), "DM change cutoff (unit of pulse width)")
			("widthcutoff", value<float>()->default_value(0.1), "Pulse width cutoff (s)")
			("snrcutoff", value<float>()->default_value(7), "S/N cutoff")
			("clean,c", "Remove candidates by dm width snr cutoff")
			("nosumif", "Do not sum polarizations")
			("cont", "Input files are contiguous")
			("wts", "Apply DAT_WTS")
			("scloffs", "Apply DAT_SCL and DAT_OFFS")
			("zero_off", "Apply ZERO_OFF")
			("psrfits", "Input psrfits format data")
			("input,f", value<std::vector<std::string>>()->multitoken()->composing(), "Input files");

	positional_options_description pos_desc;
	pos_desc.add("input", -1);
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

	if (vm.count("arch") != 0)
	{
		if (vm.count("template") == 0)
		{
			BOOST_LOG_TRIVIAL(error)<<"Error: no template file"<<endl;
			return -1;
		}
		else
		{
			std::ifstream tmp(vm["template"].as<std::string>().c_str());
			if (!tmp.good())
			{
				BOOST_LOG_TRIVIAL(error)<<"Error: can not open file "<<vm["template"].as<std::string>()<<endl;
				return -1;
			}
		}
	}

	if (vm.count("candfile") == 0)
	{
		BOOST_LOG_TRIVIAL(error)<<"Error: no candfile file"<<endl;
		return -1;
	}
	else
	{
		std::ifstream tmp(vm["candfile"].as<std::string>().c_str());
		if (!tmp.good())
		{
			BOOST_LOG_TRIVIAL(error)<<"Error: can not open file "<<vm["candfile"].as<std::string>()<<endl;
			return -1;
		}
	}

	if (vm.count("input") == 0)
	{
		BOOST_LOG_TRIVIAL(error)<<"Error: no input file"<<std::endl;
		return -1;
	}

	bool contiguous = vm.count("cont");
	std::string rootname = vm["candfile"].as<std::string>().substr(vm["candfile"].as<std::string>().find_last_of("/") + 1);
	rootname.erase(rootname.find(".cands"), 6);
	num_threads = vm["threads"].as<unsigned int>();
	float snr_threhold = vm["snrcutoff"].as<float>();

	bool apply_wts = false;
	bool apply_scloffs = false;
	bool apply_zero_off = false;

	if (vm.count("wts"))
		apply_wts = true;
	if (vm.count("scloffs"))
		apply_scloffs = true;
	if (vm.count("zero_off"))
		apply_zero_off = true;

	vector<std::string> fnames = vm["input"].as<std::vector<std::string>>();

	PSRDataReader * reader;

	if (vm.count("psrfits"))
		reader= new PsrfitsReader;
	else
		reader= new FilterbankReader;

	if (vm.count("srcname"))
		reader->source_name = vm["srcname"].as<std::string>();
	if (vm.count("ibeam"))
		reader->beam = std::to_string(vm["ibeam"].as<int>());
	if (vm.count("telescope"))
		reader->telescope = vm["telescope"].as<std::string>();
	if (vm.count("ra"))
		reader->ra = vm["ra"].as<std::string>();
	if (vm.count("dec"))
		reader->dec = vm["dec"].as<std::string>();

	reader->fnames = fnames;
	reader->sumif = vm.count("nosumif") ? false : true;
	reader->contiguous = contiguous;
	reader->verbose = verbose;
	reader->apply_scloffs = apply_scloffs;
	reader->apply_wts = apply_wts;
	reader->apply_zero_off = apply_zero_off;
	reader->check();
	reader->read_header();

	BOOST_LOG_TRIVIAL(info)<<"initialize...";

	/** downsample */
	int td = vm["td"].as<int>();
	int fd = vm["fd"].as<int>();

	/** rfi */
	std::vector<std::pair<double, double>> zaplist;
	if (vm.count("zap"))
	{
		std::vector<double> zap_opts = vm["zap"].as<std::vector<double>>();
		for (auto opt=zap_opts.begin(); opt!=zap_opts.end(); ++opt)
		{
			zaplist.push_back(pair<double, double>(*(opt+0), *(opt+1)));
			advance(opt, 1);
		}
	}

	long int nchans = reader->nchans;
	double tsamp = reader->tsamp * td;
	int nifs = reader->sumif ? 1 : reader->nifs;
	long int nsblk = reader->nsblk / td * td;

	long double start_mjd = reader->start_mjd.to_day();

	int nwidth = vm["nw"].as<int>()/2;

	/* read candidates from  */
	std::vector<Candidate> cands;

	string filename = vm["candfile"].as<string>();
	string line;
	ifstream candfile(filename);

	BOOST_LOG_TRIVIAL(info)<<"read candfile "<<filename;
	
	Candidate cand;
	cand.mjd_start = start_mjd;
	cand.obsinfo["Date"] = std::to_string(start_mjd);
	stringstream ss_ibeam;
	if (vm.count("incoherent"))
		ss_ibeam << "ifbf" << setw(5) << setfill('0') << reader->beam;
	else
		ss_ibeam << "cfbf" << setw(5) << setfill('0') << reader->beam;
	std::string s_ibeam = ss_ibeam.str();
	cand.obsinfo["Beam"] = s_ibeam;
	cand.obsinfo["Source_name"] = reader->source_name;
	cand.obsinfo["RA"] = reader->ra;
	cand.obsinfo["DEC"] = reader->dec;
	cand.obsinfo["Telescope"] = reader->telescope;
	double gl = 0., gb = 0.;
#ifdef HAVE_SOFA
	get_gl_gb(gl, gb, reader->ra, reader->dec);
#endif
	cand.obsinfo["GL"] = to_string(gl);
	cand.obsinfo["GB"] = to_string(gb);
	double ymw16_maxdm = 0.;
	ymw16_maxdm = get_maxdm_ymw16(gl, gb);
	cand.obsinfo["MaxDM_YMW16"] = to_string(ymw16_maxdm);
	cand.tbin = tsamp;
	cand.frequencies.resize(nchans);
	memcpy(cand.frequencies.data(), reader->frequencies.data(), sizeof(double)*nchans);

	BOOST_LOG_TRIVIAL(info)<<"create meta file";
	json meta;
	meta = {
		{"Telescope", cand.obsinfo["Telescope"]},
		{"Beam", cand.obsinfo["Beam"]},
		{"RA", cand.obsinfo["RA"]},
		{"DEC", cand.obsinfo["DEC"]},
		{"Source_name", cand.obsinfo["Source_name"]},
		{"GL", cand.obsinfo["GL"]},
		{"GB", cand.obsinfo["GB"]},
		{"MaxDM_YMW16", cand.obsinfo["MaxDM_YMW16"]}
	};

	std::string s_meta = meta.dump();
	std::ofstream metafile;
	metafile.open(rootname+"_replot.json");
	metafile << s_meta << std::endl;
	metafile.close();

	double dm1delay = std::abs(Candidate::dmdelay(1., cand.frequencies.front(), cand.frequencies.back()));

	while (getline(candfile, line))
	{
		boost::trim(line);
		if (line.rfind("#", 0) == 0) continue;
		
		std::vector<std::string> parameters;
		boost::split(parameters, line, boost::is_any_of("\t "), boost::token_compress_on);

		cand.s_beam = parameters[0];
		cand.id = std::stoi(parameters[1]);
		cand.mjd = std::stold(parameters[2]);
		if (vm.count("dm") == 0)
			cand.dm = std::stod(parameters[3]);
		else
			cand.dm = vm["dm"].as<double>();
		cand.width_orig = cand.width = std::stod(parameters[4])/1000.;
		cand.snr = std::stof(parameters[5]);
		cand.fl = std::stod(parameters[6]);
		cand.fh = std::stod(parameters[7]);
		cand.pngname = parameters[8];
		cand.planid = std::stoi(parameters[9]);
		cand.rawdata = parameters[10];

		double dmdelay = std::abs(Candidate::dmdelay(cand.dm, cand.frequencies.front(), cand.frequencies.back()));
		int nbin = (dmdelay+2*nwidth*cand.width)/tsamp;
		if (vm.count("pow2bin"))
			nbin = std::pow(2, std::ceil(std::log2(nbin)));
		cand.npol = nifs;
		cand.nchan = nchans;
		cand.nbin = nbin;

		if (cand.mjd+cand.nbin*cand.tbin/86400. > start_mjd)
			cands.push_back(cand);
	}

	long int count = 0;
	long int intcnt = 0;

	std::vector<int> candscnt(cands.size(), 0);
	for (long int k=0; k<cands.size(); k++)
	{
		candscnt[k] = cands[k].nbin;
	}

	std::vector<long double> candsmjd;
	candsmjd.push_back(0.);

	std::vector<bool> candsisfirst(cands.size(), true);

	RingBuffer ringbuffer;
	ringbuffer.resize(*std::max_element(candscnt.begin(), candscnt.end())+1, nifs*nchans);
	
	BOOST_LOG_TRIVIAL(info)<<"start processing... ";

	DataBuffer<float> databuffer(nsblk, nifs*nchans);
	Downsample downsample(td, 1);
	downsample.prepare(databuffer);

	while (!reader->is_end)
	{
		bool hit = false;
		for (long int k=0; k<cands.size(); k++)
		{
			double left = std::max(count*tsamp, double(cands[k].mjd-start_mjd)*86400.-2*nwidth*cands[k].width_orig);
			double right = std::min((count+nsblk/td)*tsamp, double(cands[k].mjd-start_mjd)*86400.+cands[k].nbin*tsamp+2*nwidth*cands[k].width_orig);

			if (left < right) hit = true;
		}

		if (hit)
		{
			reader->read_data(databuffer, nsblk);
			DataBuffer<float> *data = downsample.run(databuffer);
			
			for (long int i=0; i<nsblk/td; i++)
			{
				count++;
				ringbuffer.append(data->buffer.begin() + i * nifs * nchans, data->buffer.begin() + (i + 1) * nifs * nchans);

				for (long int k=0; k<cands.size(); k++)
				{
					if (!cands[k].captured)
					{
						if (candscnt[k] && count*tsamp > ((cands[k].mjd-start_mjd)*86400.-nwidth*cands[k].width))
						{
							if (candsisfirst[k])
							{
								candsisfirst[k] = false;
								cands[k].mjd_start = start_mjd+(count-1)*tsamp/86400.;
							}

							candscnt[k]--;
						}
						else if (candscnt[k] == 0)
						{
							cands[k].resize(cands[k].npol, cands[k].nchan, cands[k].nbin);

							std::vector<float> temp;
							ringbuffer.read(temp, count-cands[k].nbin-1, cands[k].nbin);

							transpose_pad<float>(cands[k].data.data(), temp.data(), cands[k].nbin, cands[k].npol*cands[k].nchan);

							temp.clear();
							temp.shrink_to_fit();

							if (vm.count("zdot"))
							{
								BOOST_LOG_TRIVIAL(debug)<<"zero-DM matched filter...";
								cands[k].zdot();
							}

							BOOST_LOG_TRIVIAL(debug)<<"calculate spectra stats...";
							cands[k].get_stats();

							BOOST_LOG_TRIVIAL(debug)<<"dedisperse at DM="<<cands[k].dm<<"...";
							cands[k].dedisperse(vm.count("coherent"));
							if (vm.count("kadane"))
								cands[k].shrink_to_fit(2*nwidth, vm["kadane"].as<std::vector<int>>()[0]);
							else
								cands[k].shrink_to_fit(2*nwidth);

							BOOST_LOG_TRIVIAL(debug)<<"remove rfi...";
							cands[k].azap(vm["zapthre"].as<float>());

							BOOST_LOG_TRIVIAL(debug)<<"manually zap channels...";
							cands[k].zap(zaplist);

							BOOST_LOG_TRIVIAL(debug)<<"calculate spectra stats...";
							cands[k].get_stats();

							BOOST_LOG_TRIVIAL(debug)<<"normalize...";
							cands[k].normalize();

							if (vm.count("clip"))
							{
								BOOST_LOG_TRIVIAL(debug)<<"clip outliers...";
								cands[k].clip(vm["clip"].as<std::vector<int>>()[0], vm["clip"].as<std::vector<int>>()[1], vm["clip"].as<std::vector<int>>()[2]);
							}

							if (vm.count("kadane"))
							{
								BOOST_LOG_TRIVIAL(debug)<<"perform kadane filter...";
								cands[k].kadaneF(vm["kadane"].as<std::vector<int>>()[0], vm["kadane"].as<std::vector<int>>()[1], vm["kadane"].as<std::vector<int>>()[2]);
							}

							BOOST_LOG_TRIVIAL(debug)<<"downsample...";
							cands[k].downsample(1, vm["fd"].as<int>());
							
							if (!vm.count("nopt"))
							{
								BOOST_LOG_TRIVIAL(debug)<<"hough transform...";
								cands[k].optimize();

								BOOST_LOG_TRIVIAL(debug)<<"boxcar filter...";
								cands[k].matched_filter(vm["snrloss"].as<float>());
							}
							else
							{
								BOOST_LOG_TRIVIAL(debug)<<"peak search...";
								cands[k].peak_search(vm["snrloss"].as<float>());
							}

							if (vm.count("clean") && (cands[k].snr<snr_threhold || closestdist(candsmjd.begin(), candsmjd.end(), cands[k].mjd)*86400. < cands[k].width || cands[k].dm_maxsnr < vm["dmcutoff"].as<float>() || cands[k].width > vm["widthcutoff"].as<float>() || std::abs(cands[k].dm_maxsnr-cands[k].dm)*dm1delay > vm["ddmcutoff"].as<float>()*cands[k].width))
							{
								BOOST_LOG_TRIVIAL(info)<<"candidate removed id="<<cands[k].id;
							}
							else
							{
								BOOST_LOG_TRIVIAL(info)<<"candidate preserved id="<<cands[k].id;

								candsmjd.push_back(cands[k].mjd);

								BOOST_LOG_TRIVIAL(debug)<<"de-dedisperse...";
								cands[k].dededisperse(vm.count("coherent"));

								if (vm.count("arch"))
								{
									BOOST_LOG_TRIVIAL(debug)<<"save to archive...";
									cands[k].save2ar(vm["template"].as<std::string>());
								}

								cands[k].dm = cands[k].dm_maxsnr;

								BOOST_LOG_TRIVIAL(debug)<<"dedisperse at DM="<<cands[k].dm<<"...";
								cands[k].dedisperse(vm.count("coherent"));

								BOOST_LOG_TRIVIAL(debug)<<"save to png...";
								cands[k].save2png(rootname, snr_threhold);
							}

							cands[k].close();

							cands[k].captured = true;
						}
					}
				}
			}
		}
		else
		{
			reader->read_data(databuffer, nsblk, true);
			for (long int i=0; i<nsblk/td; i++)
			{
				count++;
				ringbuffer.append();
			}
		}
	}

	BOOST_LOG_TRIVIAL(info)<<"done!";

	return 0;
}