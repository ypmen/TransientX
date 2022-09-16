/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-09-01 08:40:21
 * @modify date 2022-09-01 08:40:21
 * @desc [description]
 */

#include "filterbankreader.h"
#include "logging.h"
#include "mjd.h"
#include "utils.h"

FilterbankReader::FilterbankReader()
{
	nsblk = 1024;

	ntot = 0;
	count = 0;
	ns_filn = 0;

	ifile_cur = 0;
	isubint_cur = 0;
	isample_cur = 0;

	update_file = true;
	update_subint = true;
}

FilterbankReader::~FilterbankReader()
{
}

void FilterbankReader::check()
{
	size_t nfil = fnames.size();
	BOOST_LOG_TRIVIAL(info)<<"scan "<<nfil<<" filterbank files";

	fil.resize(nfil);
	for (long int i=0; i<nfil; i++)
	{
		fil[i].filename = fnames[i];
	}

	for (long int i=0; i<nfil; i++)
	{
		fil[i].read_header();
		nsamples += fil[i].nsamples;
		MJD tstart(fil[i].tstart);
		mjd_starts.push_back(tstart);
		mjd_ends.push_back(tstart+fil[i].nsamples*fil[i].tsamp);
	}

	// check continuity
	idmap = argsort(mjd_starts);

	for (long int i=0; i<nfil-1; i++)
	{
		if (abs((mjd_ends[idmap[i]]-mjd_starts[idmap[i+1]]).to_second())>0.5*fil[idmap[i]].tsamp)
		{
			if (contiguous)
			{
				BOOST_LOG_TRIVIAL(warning)<<"time not contiguous";
			}
			else
			{
				BOOST_LOG_TRIVIAL(error)<<"time not contiguous";
				exit(-1);
			}
		}
	}

	// update ifile, isubint, isample
	for (size_t idxn=ifile_cur; idxn<nfil; idxn++)
	{
		size_t n = idmap[idxn];

		long int nseg = ceil(1.*fil[n].nsamples/nsblk);

		size_t ns_filn_tmp = 0;

		for (size_t s=0; s<nseg; s++)
		{
			for (size_t i=0; i<nsblk; i++)
			{
				if (count == skip_start)
				{
					ifile_cur = idxn;
					isubint_cur = s;
					isample_cur = i;

					return;
				}

				count++;

				if (++ns_filn_tmp == fil[n].nsamples)
				{
					goto next;
				}
			}
		}
		next:;
	}
}

void FilterbankReader::read_header()
{
	BOOST_LOG_TRIVIAL(info)<<"read header";

	if (beam.empty())
	{
		if (fil[0].ibeam != 0)
		{
			BOOST_LOG_TRIVIAL(info)<<"read beam_id from file";
			beam = std::to_string(fil[0].ibeam);
		}
	}

	if (source_name.empty())
	{
		if (strcmp(fil[0].source_name, "") != 0)
		{
			BOOST_LOG_TRIVIAL(info)<<"read source name from file";
			source_name = fil[0].source_name;
		}
	}

	if (telescope.empty())
	{
		BOOST_LOG_TRIVIAL(info)<<"read telescope from file";
		get_telescope_name(fil[0].telescope_id, telescope);
	}

	if (ra.empty() or dec.empty())
	{
		get_s_radec(fil[0].src_raj, fil[0].src_dej, ra, dec);
	}

	nchans = fil[0].nchans;
	tsamp = fil[0].tsamp;
	nifs = fil[0].nifs;

	start_mjd = MJD(fil[idmap[0]].tstart);

	frequencies.resize(nchans, 0.);
	std::memcpy(frequencies.data(), fil[0].frequency_table, sizeof(double) * nchans);
}

size_t FilterbankReader::read_data(DataBuffer<float> &databuffer, size_t ndump, bool virtual_reading)
{
	assert(databuffer.buffer.size() > 0);

	size_t bcnt1 = 0;

	size_t nfil = fil.size();
	for (size_t idxn=ifile_cur; idxn<nfil; idxn++)
	{
		size_t n = idmap[idxn];
		size_t nseg = ceil(1.*fil[n].nsamples/nsblk);

		if (update_file)
		{
			ns_filn = 0;
		}
		update_file = false;

		double zero_off = 0.;

		for (size_t s=isubint_cur; s<nseg; s++)
		{
			if (verbose)
			{
				cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
				cerr<<"("<<100.*count/nsamples<<"%)";
			}

			if (!virtual_reading)
			{
				if (update_subint)
				{
					fil[n].read_data(ns_filn, nsblk);
				}
				update_subint = false;
			}

			for (size_t i=isample_cur; i<nsblk; i++)
			{
				if (!virtual_reading)
				{
						if (!sumif or nifs == 1)
						{
							for (size_t k=0; k<nifs; k++)
							{
								for (size_t j=0; j<nchans; j++)
								{
									databuffer.buffer[bcnt1*nifs*nchans+k*nchans+j] = ((unsigned char *)(fil[n].data))[i*nifs*nchans+k*nchans+j] - zero_off;
								}
							}
						}
						else
						{
							for (size_t j=0; j<nchans; j++)
							{
								float xx = ((unsigned char *)(fil[n].data))[i*nifs*nchans+0*nchans+j] - zero_off;
								float yy = ((unsigned char *)(fil[n].data))[i*nifs*nchans+1*nchans+j] - zero_off;

								databuffer.buffer[bcnt1*nchans+j] = xx + yy;
							}
						}
				}

				bcnt1++;
				count++;
				ntot++;
				ns_filn++;

				if (count == nsamples - skip_end)
				{
					if (ns_filn == fil[n].nsamples)
					{
						fil[n].free();
						update_file = true;
						update_subint = true;
					}

					is_end = true;
					if (verbose)
					{
						cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
						cerr<<"("<<100.*count/nsamples<<"%)";
					}
					return bcnt1;
				}

				if (bcnt1 == ndump)
				{
					ifile_cur = idxn;
					isubint_cur = s;
					isample_cur = i+1;

					if (isample_cur == nsblk)
					{
						isample_cur = 0;
						isubint_cur++;
						update_subint = true;
					}

					if (ns_filn == fil[n].nsamples)
					{
						isample_cur = 0;
						isubint_cur = 0;
						ifile_cur++;

						fil[n].free();
						update_file = true;
						update_subint = true;
					}

					return bcnt1;
				}

				if (ns_filn == fil[n].nsamples)
				{
					goto next;
				}
			}
			isample_cur = 0;
			update_subint = true;
		}
		next:
		isample_cur = 0;
		isubint_cur = 0;
		fil[n].free();
		update_file = true;
		update_subint = true;
	}

	is_end = true;
	if (verbose)
	{
		cerr<<"\r\rfinish "<<setprecision(2)<<fixed<<tsamp*count<<" seconds ";
		cerr<<"("<<100.*count/nsamples<<"%)";
	}
	return bcnt1;
}

void FilterbankReader::get_filterbank_template(Filterbank &filtem)
{
	filtem = fil[idmap[0]];
}