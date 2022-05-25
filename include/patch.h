/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2022-05-13 11:59:35
 * @modify date 2022-05-13 11:59:35
 * @desc [description]
 */

#ifndef PATCH_H
#define PATCH_H

#include "databuffer.h"

class Patch : public DataBuffer<float>
{
public:
	Patch();
	~Patch();
	void prepare(DataBuffer<float> &databuffer);
	DataBuffer<float> * filter(DataBuffer<float> &databuffer);
	DataBuffer<float> * filter2(DataBuffer<float> &databuffer);
	DataBuffer<float> * get(){return this;}
public:
	string filltype;
	float width;
	float threshold;
	float killrate;
};

#endif /* PATCH_H */
