/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-11-02 16:24:33
 * @modify date 2020-11-02 16:24:33
 * @desc [description]
 */

#include "assert.h"

#include "archivewriter_tx.h"

using namespace TransientX;

ArchiveWriter::ArchiveWriter()
{
	mode = Integration::FOLD;
	ibeam = 1;
	npol = 0;
	nchan = 0;
	nbin = 0;
	nsubint = 0;
	tbin = 0.;
	dm = 0.;
}

ArchiveWriter::~ArchiveWriter()
{
	if (fits.fptr != NULL)
	{
		fits.close();
	}
}

void ArchiveWriter::close()
{
	if (fits.fptr != NULL)
	{
		fits.close();
	}    
}

void ArchiveWriter::prepare()
{
	assert(mode == Integration::FOLD);

	fits.primary.start_mjd = start_mjd;
	strcpy(fits.primary.src_name, src_name.c_str());
	strcpy(fits.primary.ra, ra.c_str());
	strcpy(fits.primary.dec, dec.c_str());

	fits.filename = rootname + ".ar";
	strcpy(fits.primary.obs_mode, "PSR");
	strcpy(fits.primary.ibeam, to_string(ibeam).c_str());

	fits.primary.chan_dm = dm;

	fits.subint.mode = Integration::FOLD;
	fits.subint.dtype = Integration::SHORT;
	fits.subint.npol = npol;
	fits.subint.nchan = nchan;
	fits.subint.nbin = nbin;
	fits.subint.tbin = tbin;
	fits.subint.dm = dm;

	fits.parse_template(template_file);
	fits.primary.unload(fits.fptr);
	fits.subint.unload_header(fits.fptr);

	it.mode = Integration::FOLD;
	it.dtype = Integration::SHORT;
}

void ArchiveWriter::run(vector<float> &profile, int np, int nc, int nb, double fold_period = 0., double tsubint = 0., double offs_sub = 0.)
{
	assert(np == npol);
	assert(nc == nchan);
	assert(nb == nbin);

	it.load_data(&profile[0], npol, nchan, nbin);
	it.load_frequencies(&frequencies[0], nchan);
	it.folding_period = fold_period;
	it.tsubint = tsubint;
	it.offs_sub = offs_sub;
	fits.subint.unload_integration(fits.fptr, it);
	nsubint++;
}

void ArchiveWriter::run(vector<float> &profiles, int ns, int np, int nc, int nb, vector<double> &fold_periods, vector<double> &tsubints, vector<double> &offs_subs)
{
	assert(np == npol);
	assert(nc == nchan);
	assert(nb == nbin);

	for (long int k=0; k<ns; k++)
	{
		it.load_data(&profiles[0]+k*npol*nchan*nbin, npol, nchan, nbin);
		it.load_frequencies(&frequencies[0], nchan);
		it.folding_period = fold_periods[k];
		it.tsubint = tsubints[k];
		it.offs_sub = offs_subs[k];
		fits.subint.unload_integration(fits.fptr, it);
		nsubint++;
	}
}
