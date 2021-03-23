/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-03-19 19:26:17
 * @modify date 2021-03-19 19:26:17
 * @desc [description]
 */

#include "preprocess.h"
#include "utils.h"
#include "dedisperse.h"

void Preprocess::prepare(DataBuffer<short> &databuffer)
{
    nsamples = databuffer.nsamples/td;
    nchans = databuffer.nchans/fd;

    resize(nsamples, nchans);
    tsamp = databuffer.tsamp*td;

    fill(frequencies.begin(), frequencies.end(), 0.);
    for (long int j=0; j<nchans; j++)
    {
        for (long int k=0; k<fd; k++)
        {
            frequencies[j] += databuffer.frequencies[j*fd+k];
        }
        frequencies[j] /= fd;
    }

    closable = false;
}

DataBuffer<float> * Preprocess::run(DataBuffer<short> &databuffer)
{
    if (closable) open();

    std::vector<short> temp(databuffer.nchans*databuffer.nsamples, 0);
    
    transpose_pad<short>(temp.data(), databuffer.buffer.data(), databuffer.nsamples, databuffer.nchans);

    std::vector<float> bufferT(nchans*nsamples, 0.);

    int nwin = ceil(width/databuffer.tsamp);
    nwin = nwin/2*2+1;

    for (long int l=0; l<fd; l++)
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_threads)
#endif
        for (long int j=0; j<nchans; j++)
        {
            std::vector<size_t> hist1(512, 0);
            std::vector<size_t> hist2(1024, 0);
            std::vector<short> med(databuffer.nsamples, 0);
            /* remove baseline */
            if (nwin >= 3 && nwin<=databuffer.nsamples)
            {
                int cnt = 0;
                for (long int i=0; i<(nwin-1)/2; i++)
                {
                    hist1[temp[(j*fd+l)*databuffer.nsamples+i]]++;
                    cnt++;
                }

                int median = -1;
                int res = 0;
                int tmpcnt = 0;
                for (long int k=0; k<512; k++)
                {
                    tmpcnt += hist1[k];
                    if (tmpcnt > cnt/2 && median < 0)
                    {
                        median = k;
                        res = tmpcnt-cnt/2-1;
                    }
                }

                int left = 0;
                int right = (nwin-1)/2;
                for (long int i=0; i<nwin/2+1; i++)
                {
                    hist1[temp[(j*fd+l)*databuffer.nsamples+right]]++;
                    cnt++;
                    /* add right */
                    if (temp[(j*fd+l)*databuffer.nsamples+right] > median)
                    {
                        if (cnt/2-(cnt-1)/2 == 1)
                        {
                            forward(hist1, median, res);
                        }
                    }
                    else if (temp[(j*fd+l)*databuffer.nsamples+right] < median)
                    {
                        if (cnt/2-(cnt-1)/2 == 0)
                        {
                            backward(hist1, median, res);
                        }
                    }
                    else
                    {
                        if (cnt/2-(cnt-1)/2 == 0)
                        {
                            backward(hist1, median, res);
                        }
                    }

                    med[i] = median;

                    right++;
                }
                for (long int i=nwin/2+1; i<databuffer.nsamples-(nwin-1)/2; i++)
                {
                    hist1[temp[(j*fd+l)*databuffer.nsamples+right]]++;
                    cnt++;
                    /* add right */
                    if (temp[(j*fd+l)*databuffer.nsamples+right] > median)
                    {
                        if (cnt/2-(cnt-1)/2 == 1)
                        {
                            forward(hist1, median, res);
                        }
                    }
                    else if (temp[(j*fd+l)*databuffer.nsamples+right] < median)
                    {
                        if (cnt/2-(cnt-1)/2 == 0)
                        {
                            backward(hist1, median, res);
                        }
                    }
                    else
                    {
                        if (cnt/2-(cnt-1)/2 == 0)
                        {
                            backward(hist1, median, res);
                        }
                    }

                    hist1[temp[(j*fd+l)*databuffer.nsamples+left]]--;
                    cnt--;
                    /* remove left */
                    if (temp[(j*fd+l)*databuffer.nsamples+left] > median)
                    {
                        if ((cnt+1)/2-cnt/2 == 1)
                        {
                            backward(hist1, median, res);
                        }
                    }
                    else if (temp[(j*fd+l)*databuffer.nsamples+left] < median)
                    {
                        if ((cnt+1)/2-cnt/2 == 0)
                        {
                            forward(hist1, median, res);
                        }
                    }
                    else
                    {
                        if (res == 0)
                        {
                            if (hist1[median] > 0)
                            {
                                if ((cnt+1)/2-cnt/2 == 0)
                                {
                                    forward(hist1, median, res);
                                }
                            }
                            else
                            {
                                if ((cnt+1)/2-cnt/2 == 0)
                                {
                                    forward(hist1, median, res);
                                }
                                else
                                {
                                    backward(hist1, median, res);
                                }
                            }
                        }
                        else
                        {
                            res--;
                            if ((cnt+1)/2-cnt/2 == 1)
                            {
                                backward(hist1, median, res);
                            }
                        }
                    }

                    left++;
                    right++;

                    med[i] = median;
                }
                for (long int i=databuffer.nsamples-(nwin-1)/2; i<databuffer.nsamples; i++)
                {
                    hist1[temp[(j*fd+l)*databuffer.nsamples+left]]--;
                    cnt--;
                    /* remove left */
                    if (temp[(j*fd+l)*databuffer.nsamples+left] > median)
                    {
                        if ((cnt+1)/2-cnt/2 == 1)
                        {
                            backward(hist1, median, res);
                        }
                    }
                    else if (temp[(j*fd+l)*databuffer.nsamples+left] < median)
                    {
                        if ((cnt+1)/2-cnt/2 == 0)
                        {
                            forward(hist1, median, res);
                        }
                    }
                    else
                    {
                        if (res == 0)
                        {
                            if (hist1[median] > 0)
                            {
                                if ((cnt+1)/2-cnt/2 == 0)
                                {
                                    forward(hist1, median, res);
                                }
                            }
                            else
                            {
                                if ((cnt+1)/2-cnt/2 == 0)
                                {
                                    forward(hist1, median, res);
                                }
                                else
                                {
                                    backward(hist1, median, res);
                                }
                            }
                        }
                        else
                        {
                            res--;
                            if ((cnt+1)/2-cnt/2 == 1)
                            {
                                backward(hist1, median, res);
                            }
                        }
                    }

                    left++;

                    med[i] = median;
                }
            }

            /* zapping based on skewness and kurtosis */
            double mean = 0.;
            for (long int i=0; i<databuffer.nsamples; i++)
            {
                temp[(j*fd+l)*databuffer.nsamples+i] -= med[i];
                temp[(j*fd+l)*databuffer.nsamples+i] += 512;

                mean += temp[(j*fd+l)*databuffer.nsamples+i];

                hist2[temp[(j*fd+l)*databuffer.nsamples+i]]++;
            }
            mean /= databuffer.nsamples;

            double mean4 = 0.;
            double mean3 = 0.;
            double var = 0.;
            long int q1 = -1;
            long int q2 = -1;
            long int q3 = -1;
            int cnt = 0;
            for (long int k=0; k<1024; k++)
            {
                double tmp = k-mean;
                double tmp2 = tmp*tmp;

                mean4 += hist2[k]*tmp2*tmp2;
                mean3 += hist2[k]*tmp*tmp2;
                var += hist2[k]*tmp2;
                cnt += hist2[k];

                if (cnt > databuffer.nsamples/4 && q1 < 0)
                {
                    q1 = k;
                }
                if (cnt > databuffer.nsamples/2 && q2 < 0)
                {
                    q2 = k;
                }
                if (cnt > databuffer.nsamples*3/4 && q3 < 0)
                {
                    q3 = k;
                }
            }
            mean4 /= databuffer.nsamples;
            mean3 /= databuffer.nsamples;
            var /= databuffer.nsamples;

            double chkurtosis = mean4/(var*var)-3.;
            double chskewness = mean3/(var*std::sqrt(var));
            double chmean = q2;
            double chstd = (q3-q1)/1.349;
            chstd = chstd==0. ? 1.:chstd;

            if (std::abs(chkurtosis) < thresig*std::sqrt(24./databuffer.nsamples) && std::abs(chskewness) < thresig*std::sqrt(6./databuffer.nsamples))
            {
                for (long int m=0; m<td; m++)
                {
                    for (long int i=0; i<nsamples; i++)
                    {
                        bufferT[j*nsamples+i] += (temp[(j*fd+l)*databuffer.nsamples+(i*td+m)]-chmean)/chstd;
                    }
                }
            }
            else
            {
                for (long int m=0; m<td; m++)
                {
                    for (long int i=0; i<nsamples; i++)
                    {
                        bufferT[j*nsamples+i] = 0.;
                    }
                }
            }
        }
    }

    transpose_pad<float>(buffer.data(), bufferT.data(), nchans, nsamples);

    if (td == 1 && fd == 1)
        equalized = true;
    else
        equalized = false;

    counter += nsamples;

    databuffer.isbusy = false;
    isbusy = true;

    if (databuffer.closable) databuffer.close();

    return this;
}