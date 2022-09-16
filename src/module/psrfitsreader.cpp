/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-08-20 23:10:30
 * @modify date 2022-08-20 23:10:30
 * @desc [description]
 */

#include "psrfitsreader.h"
#include "logging.h"
#include "mjd.h"
#include "utils.h"

PsrfitsReader::PsrfitsReader()
{
	ntot = 0;
	count = 0;
	ns_psfn = 0;

	ifile_cur = 0;
	isubint_cur = 0;
	isample_cur = 0;

	update_file = true;
	update_subint = true;
}

PsrfitsReader::~PsrfitsReader()
{
}

void PsrfitsReader::check()
{
	size_t npsf = fnames.size();
	BOOST_LOG_TRIVIAL(info)<<"scan "<<npsf<<" psrfits files";

	psf.resize(npsf);
	for (size_t i=0; i<npsf; i++)
	{
		psf[i].filename = fnames[i];
	}

	for (size_t i=0; i<npsf; i++)
	{
		psf[i].open();
		psf[i].primary.load(psf[i].fptr);
		psf[i].load_mode();
		psf[i].subint.load_header(psf[i].fptr);
		nsamples += psf[i].subint.nsamples;
		mjd_starts.push_back(psf[i].primary.start_mjd);
		mjd_ends.push_back(psf[i].primary.start_mjd+psf[i].subint.nsamples*psf[i].subint.tbin);
		psf[i].close();
	}

	// check continuity
	idmap = argsort(mjd_starts);

	for (size_t i=0; i<npsf-1; i++)
	{
		if (abs((mjd_ends[idmap[i]]-mjd_starts[idmap[i+1]]).to_second())>0.5*psf[idmap[i]].subint.tbin)
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
	for (size_t idxn=ifile_cur; idxn<npsf; idxn++)
	{
		size_t n = idmap[idxn];

		size_t ns_psfn_tmp = 0;

		psf[n].open();
		psf[n].primary.load(psf[n].fptr);
		psf[n].load_mode();
		psf[n].subint.load_header(psf[n].fptr);

		for (size_t s=0; s<psf[n].subint.nsubint; s++)
		{
			for (size_t i=0; i<psf[n].subint.nsblk; i++)
			{
				if (count == skip_start)
				{
					ifile_cur = idxn;
					isubint_cur = s;
					isample_cur = i;

					return;
				}

				count++;

				if (++ns_psfn_tmp == psf[n].subint.nsamples)
				{
					goto next;
				}
			}
		}
		next:
		psf[n].close();
	}
}

void PsrfitsReader::read_header()
{
	BOOST_LOG_TRIVIAL(info)<<"read header";

	psf[0].open();
	psf[0].primary.load(psf[0].fptr);
	psf[0].load_mode();
	psf[0].subint.load_header(psf[0].fptr);

	if (psf[0].mode != Integration::SEARCH)
	{
		BOOST_LOG_TRIVIAL(error)<<"mode is not SEARCH";
		exit(-1);
	}

	if (beam.empty())
	{
		if (strcmp(psf[0].primary.ibeam, "") != 0)
		{
			BOOST_LOG_TRIVIAL(info)<<"read beam_id from file";
			beam = psf[0].primary.ibeam;
		}
	}
	
	if (source_name.empty())
	{
		if (strcmp(psf[0].primary.src_name, "") != 0)
		{
			BOOST_LOG_TRIVIAL(info)<<"read source name from file";
			source_name = psf[0].primary.src_name;
		}
	}

	if (telescope.empty())
	{
		if (strcmp(psf[0].primary.telesop, "") != 0)
		{
			BOOST_LOG_TRIVIAL(info)<<"read telescope from file";
			telescope = psf[0].primary.telesop;
		}
	}

	if (ra.empty())
	{
		if (strcmp(psf[0].primary.ra, "") != 0)
		{
			BOOST_LOG_TRIVIAL(info)<<"read ra from file";
			ra = psf[0].primary.ra;
		}
	}

	if (dec.empty())
	{
		if (strcmp(psf[0].primary.dec, "") != 0)
		{
			BOOST_LOG_TRIVIAL(info)<<"read dec from file";
			dec = psf[0].primary.dec;
		}
	}

	psf[0].subint.load_integration(psf[0].fptr, 0, it);

	nchans = it.nchan;
	tsamp = psf[0].subint.tbin;
	nifs = it.npol;
	nsblk = it.nsblk;

	start_mjd = psf[idmap[0]].primary.start_mjd;

	frequencies.resize(nchans, 0.);
	std::memcpy(frequencies.data(), it.frequencies, sizeof(double)*nchans);

	psf[0].close();
}

size_t PsrfitsReader::read_data(DataBuffer<float> &databuffer, size_t ndump, bool virtual_reading)
{
	assert(databuffer.buffer.size() > 0);

	size_t bcnt1 = 0;

	size_t npsf = psf.size();
	for (size_t idxn=ifile_cur; idxn<npsf; idxn++)
	{
		size_t n = idmap[idxn];

		if (update_file)
		{
			psf[n].open();
			psf[n].primary.load(psf[n].fptr);
			psf[n].load_mode();
			psf[n].subint.load_header(psf[n].fptr);

			ns_psfn = 0;
		}
		update_file = false;

		double zero_off = 0.;
		if (apply_zero_off)
			zero_off = psf[n].subint.zero_off;

		for (size_t s=isubint_cur; s<psf[n].subint.nsubint; s++)
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
					if (apply_scloffs or apply_wts)
						psf[n].subint.load_integration(psf[n].fptr, s, it);
					else
						psf[n].subint.load_integration_data(psf[n].fptr, s, it);
				}
				update_subint = false;
			}

			for (size_t i=isample_cur; i<it.nsblk; i++)
			{
				if (!virtual_reading)
				{
					if (it.dtype == Integration::UINT8)
					{
						if (!sumif or nifs == 1)
						{
							for (size_t k=0; k<nifs; k++)
							{
								if (apply_wts)
								{
									if (apply_scloffs)
									{
										for (size_t j=0; j<nchans; j++)
										{
											databuffer.buffer[bcnt1*nifs*nchans+k*nchans+j] = it.weights[j] * ((((unsigned char *)(it.data))[i*nifs*nchans+k*nchans+j] - zero_off) * it.scales[k*nchans+j] + it.offsets[k*nchans+j]);
										}
									}
									else
									{
										for (size_t j=0; j<nchans; j++)
										{
											databuffer.buffer[bcnt1*nifs*nchans+k*nchans+j] = it.weights[j] * (((unsigned char *)(it.data))[i*nifs*nchans+k*nchans+j] - zero_off);
										}
									}
								}
								else
								{
									if (apply_scloffs)
									{
										for (size_t j=0; j<nchans; j++)
										{
											databuffer.buffer[bcnt1*nifs*nchans+k*nchans+j] = (((unsigned char *)(it.data))[i*nifs*nchans+k*nchans+j] - zero_off) * it.scales[k*nchans+j] + it.offsets[k*nchans+j];
										}
									}
									else
									{
										for (size_t j=0; j<nchans; j++)
										{
											databuffer.buffer[bcnt1*nifs*nchans+k*nchans+j] = ((unsigned char *)(it.data))[i*nifs*nchans+k*nchans+j] - zero_off;
										}
									}
								}
							}
						}
						else
						{
							if (apply_wts)
							{
								if (apply_scloffs)
								{
									for (size_t j=0; j<nchans; j++)
									{
										float xx = it.weights[j] * ((((unsigned char *)(it.data))[i*nifs*nchans+0*nchans+j] - zero_off) * it.scales[0*nchans+j] + it.offsets[0*nchans+j]);
										float yy = it.weights[j] * ((((unsigned char *)(it.data))[i*nifs*nchans+1*nchans+j] - zero_off) * it.scales[1*nchans+j] + it.offsets[1*nchans+j]);
										
										databuffer.buffer[bcnt1*nchans+j] = xx + yy;
									}
								}
								else
								{
									for (size_t j=0; j<nchans; j++)
									{
										float xx = it.weights[j] * (((unsigned char *)(it.data))[i*nifs*nchans+0*nchans+j] - zero_off);
										float yy = it.weights[j] * (((unsigned char *)(it.data))[i*nifs*nchans+1*nchans+j] - zero_off);

										databuffer.buffer[bcnt1*nchans+j] = xx + yy;
									}
								}
							}
							else
							{
								if (apply_scloffs)
								{
									for (size_t j=0; j<nchans; j++)
									{
										float xx = (((unsigned char *)(it.data))[i*nifs*nchans+0*nchans+j] - zero_off) * it.scales[0*nchans+j] + it.offsets[0*nchans+j];
										float yy = (((unsigned char *)(it.data))[i*nifs*nchans+1*nchans+j] - zero_off) * it.scales[1*nchans+j] + it.offsets[1*nchans+j];

										databuffer.buffer[bcnt1*nchans+j] = xx + yy;
									}
								}
								else
								{
									for (size_t j=0; j<nchans; j++)
									{
										float xx = ((unsigned char *)(it.data))[i*nifs*nchans+0*nchans+j] - zero_off;
										float yy = ((unsigned char *)(it.data))[i*nifs*nchans+1*nchans+j] - zero_off;

										databuffer.buffer[bcnt1*nchans+j] = xx + yy;
									}
								}
							}
						}
					}
					else
					{
						BOOST_LOG_TRIVIAL(error)<<"data type is not supported"<<endl;
					}
				}

				bcnt1++;
				count++;
				ntot++;
				ns_psfn++;

				if (count == nsamples - skip_end)
				{
					if (ns_psfn == psf[n].subint.nsamples)
					{
						psf[n].close();
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

					if (isample_cur == it.nsblk)
					{
						isample_cur = 0;
						isubint_cur++;
						update_subint = true;
					}

					if (ns_psfn == psf[n].subint.nsamples)
					{
						isample_cur = 0;
						isubint_cur = 0;
						ifile_cur++;

						psf[n].close();
						update_file = true;
						update_subint = true;
					}

					return bcnt1;
				}

				if (ns_psfn == psf[n].subint.nsamples)
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
		psf[n].close();
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

void PsrfitsReader::get_filterbank_template(Filterbank &filtem)
{
	filtem.use_frequence_table = false;
	filtem.data_type = 1;
	strcpy(filtem.rawdatafile, psf[idmap[0]].filename.c_str());
	filtem.tstart =start_mjd.to_day();
	filtem.tsamp = tsamp;
	filtem.nifs = nifs;
	if (sumif) filtem.nifs = 1;
	filtem.nchans = nchans;
	
	get_fch1_foff(filtem.fch1, filtem.foff);

	if (filtem.frequency_table != NULL) delete [] filtem.frequency_table;
	filtem.frequency_table = new double [16320];
	memcpy(filtem.frequency_table, frequencies.data(), sizeof(double)*nchans);
}