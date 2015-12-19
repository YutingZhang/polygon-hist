#include "zpolyhist.h"
#include "zutils.h"
#include <opencv/cxcore.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>
#include <cstdio>
#include <string>
#include <cstring>
#include <cmath>
#include <boost/timer.hpp>
using namespace std;
using namespace zPolyHist;

#ifdef _WIN_X64
	#pragma warning(disable:4267)
#endif

class piFun4Test : public HistFunBase {
public:
	piFun4Test() {}
	virtual output_type operator () (bin_type b, input_type n) = 0;
	virtual ~piFun4Test() {};
};

class piFun4TestLinear : public piFun4Test {
private:
	std::vector<FloatType> w;
public:
	piFun4TestLinear(size_t binnum) : w(binnum) {
		for ( size_t i=0; i<binnum; ++i ) {
			//w[i] = (double)rand()/RAND_MAX;
			w[i]=(double)i+1;
		}
	}
	virtual output_type operator () (bin_type b, input_type n) {
		//return n;
		return n*w[b];
	}
	virtual ~piFun4TestLinear() {};
};

class piFun4TestChi : public piFun4Test {
private:
	std::vector<FloatType> w;
public:
	piFun4TestChi(size_t binnum) : w(binnum) {
		for ( size_t i=0; i<binnum; ++i )
			w[i] = ((double)rand()/RAND_MAX)*10.0;
	}
	virtual output_type operator () (bin_type b, input_type n) {
		FloatType a=(n+w[b]);
		return a*a/(n-w[b]);
	}
	virtual ~piFun4TestChi() {};
};

class piFun4TestEntropy : public piFun4Test {
private:
	std::vector<FloatType> w;
public:
	piFun4TestEntropy(){}
	virtual output_type operator () (bin_type b, input_type n) {
		return (n+1)*log(n+1);
	}
	virtual ~piFun4TestEntropy() {};
};

class piFun4TestSVMChiKernel : public piFun4Test {
private:
	std::vector<FloatType> w;
	FloatArray2D	m;
public:
	piFun4TestSVMChiKernel(size_t binnum, size_t sv = 50) : w(sv) {
		m = FloatArray2D::Create(binnum,sv);
		for ( size_t k=0; k<sv; ++k ) {
			w[k] = ((FloatType)rand()/RAND_MAX);
			for ( size_t i=0; i<binnum; ++i )
				m[k][i] = ((FloatType)rand()/RAND_MAX);
		}
	}
	virtual output_type operator () (bin_type b, input_type n) {
		FloatType s = 0.0;
		for ( size_t k=0; k<w.size(); ++k ) {
			FloatType a=n+m[k][b];
			s += w[k]*(a*a/(n-m[k][b]));
		}
		return s;
	}
	virtual ~piFun4TestSVMChiKernel() {
		FloatArray2D::Release(m);
	};
};

int main_time(const int argc, const char* argv[]) {
	/*
	*	polyint	
			-m	method
			-f	input image path
			-b	binnum base, binnum = N^2
			-p	polygon group	(Can use multiple times)
			-z	polygon group zoom factor
			-t	function type
		-m:
	*	-2 - integral image
	*	-1 - naive
	*	 0 - fast auto
	*	 1 - fast edge
	*	 2 - fast naive

		-t
		1 - linear
		2 - Chi
		3 - Entropy
	*/

	int method = 0;
	int binnum_base1 = 16;
	int binnum_base2 = 16;
	int zoom = 1;
	int fun_type = 1;
	int strideX = 1, strideY = 1;
	const char _input_file[]  = "test.jpg";
	const char _output_file[] = "tmp.hist";
	const char* input_file = _input_file;
	const char* output_file = _output_file;
	bool	is_labelMap   = false;
	bool	is_textOutput = false;
	list<string> poly_files;

	for ( int i=1; i<argc; ++i ) {
		if ( argv[i][0] == '-' ) {
			const char mc = argv[i][1];
			const char* str = argv[i+1];
			++i;
			switch (mc) {
				case 'm':
					sscanf(str, "%d", &method);
					cerr << "Method: " << method << endl;
					break;
				case 'f':
					input_file = str;
					cerr << "Input: " << input_file << endl;
					break;
				case 'b':
					sscanf(str, "%d", &binnum_base1);
					binnum_base2 = binnum_base1;
					cerr << "Bin Number Base: " << binnum_base1 << '\t' << binnum_base2 << endl;
					break;
				case 'B':
					sscanf(str, "%d", &binnum_base1);
					++i;
					sscanf(argv[i], "%d", &binnum_base2);
					cerr << "Bin Number Base: " << binnum_base1 << '\t' << binnum_base2 << endl;
					break;
				case 'p':
					poly_files.push_back(string(str));
					cerr << "Add polygon group: " << poly_files.back() << endl;
					break;
				case 'z':
					sscanf(str, "%d", &zoom);
					cerr << "Zoom factor: " << zoom << endl;
					break;
				case 't':
					sscanf(str, "%d", &fun_type);
					cerr << "Function type: " << fun_type << endl;
					break;
				case 'o':
					output_file = str;
					cerr << "Output: " << output_file << endl;
					break;
				case 'd':	//Data type 'L' - label map; 'I' - image (default)
					is_labelMap = (str[0]=='L');
					cerr << "Input format: " << ((is_labelMap)?("label map"):("color image")) << endl;
					break;
				case 's':	//Stride
					sscanf(str, "%d", &strideX);
					strideY = strideX;
					cerr << "Stride: " << strideX << "\t" << strideY << endl;
					break;
				case 'S':	//Stride X Y
					sscanf(str, "%d", &strideX);
					++i;
					sscanf(argv[i], "%d", &strideY);
					cerr << "Stride: " << strideX << "\t" << strideY << endl;
					break;
				case 'O':	//Data type 'T' - text; 'B' - binary (default)
					is_textOutput = (str[0]=='T');
					cerr << "Output format: " << ((is_textOutput)?("text"):("binary")) << endl;
					break;
			}
		}
	}

	Shapes	shapes;
	{
		//Load polygons
		for ( list<string>::iterator iter = poly_files.begin();
			iter != poly_files.end(); ++iter ) {
				ifstream in(iter->c_str());
				if (!in) {
					cerr << "Can't open " << *iter << endl;
					continue;
				}
				IntPoint p;
				shapes.NewShape();
				while ( in >> p.x >> p.y ) {
					if ( p.x== -1 && p.y== -1 ) {		//use (-1,-1) as terminator
						shapes.NewPolygon();
					} else {
						shapes.AddVertex2Current(IntPoint(p.x*zoom,p.y*zoom));
					}
				}
		}
		shapes.Complete();
	}


	if ( !shapes.size() ) {
		cerr << "No polygroup file." << endl;
		exit(1);
	}



	int binnum;
	typedef UIntArray2D		binArr;
	typedef ArrayOfUIntArray2D		binArrs;
	typedef FloatArray2D	weightArr;
	typedef ArrayOfFloatArray2D		weightArrs;
	binArr arrp;
	weightArr warrp;
	if (!is_labelMap)
	{
		//Read image
		IplImage* img = cvLoadImage(input_file, CV_LOAD_IMAGE_COLOR);
		if(!img) {
			cerr << "Can't open " << input_file << endl;
			exit(1);
		}
		CvSize imgSize = cvGetSize(img);
		IplImage* imgc;

		imgc = cvCreateImage( imgSize, img->depth, img->nChannels );
		cvCvtColor( img, imgc, CV_RGB2HLS);


		int colorDepthFrac1 = 256/binnum_base1;
		int colorDepthFrac2 = 256/binnum_base2;
		int true_base1 = 256/colorDepthFrac1;
		int true_base2 = 256/colorDepthFrac2;
		binnum = true_base1*true_base2;
		arrp = binArr::Create( imgSize.width, imgSize.height ) ;
		warrp= weightArr::Create( imgSize.width, imgSize.height );

		for ( int y=0; y<imgSize.height; ++y ) {
			for ( int x=0; x<imgSize.width; ++x ) {
				CvScalar c = cvGet2D( imgc, y, x );
				arrp[y][x] = ((int)c.val[0]/colorDepthFrac1) * true_base2 +
					(int)c.val[2]/colorDepthFrac2;
				warrp[y][x]= (double)c.val[1]/255.0;
			}
		}
	}
	else
	{
		binnum = binnum_base1;
		//Read bin map
		FILE* in = fopen(input_file,"r");
		if(!in) {
			cerr << "Can't open " << input_file << endl;
			exit(1);
		}
		const size_t buffer_size = 1024*1024;
		char* lbuffer = new char[buffer_size];
		char seps[]=",";
		char* token;
		typedef std::vector<int>	intVec;
		ptrvec< intVec > content;
		while ( fgets(lbuffer,buffer_size,in) ) {
			intVec* lvec = new intVec;
			content.push_back(lvec);
			token = strtok( lbuffer, seps );
			while ( token ) {
				int b;
				if (sscanf(token, "%d", &b)>0)
					lvec->push_back(b);
				token = strtok( NULL, seps );
			}
		}
		fclose(in);

		arrp = binArr::Create( content[0].size(), content.size() ) ;
		warrp= weightArr::Create( content[0].size(), content.size() );
		for ( int y=0; y<(int)arrp.height; ++y ) {
			for ( int x=0; x<(int)arrp.width; ++x ) {
				arrp[y][x] = content[y][x];
				warrp[y][x]= 1.0;
			}
		}
	}

	cerr << "Bin num: " << binnum << endl;

	binArrs		arr  = binArrs::Create(2);
	weightArrs	warr = weightArrs::Create(2);
	arr(0)  = arrp;  arr(1)  = arrp;
	warr(0) = warrp; warr(1) = warrp;

	typedef piFun4Test FunType;

	typedef auto_ptr<FunType> FunPtr;
	FunPtr	hfun;
	if ( fun_type == 1 ) {
		hfun = FunPtr( new piFun4TestLinear(binnum) );
	} else if ( fun_type == 2 ) {
		hfun = FunPtr( new piFun4TestChi(binnum) );
	} else if ( fun_type == 3 ) {
		hfun = FunPtr( new piFun4TestEntropy() );
	} else if ( fun_type == 4 ) {
		hfun = FunPtr( new piFun4TestSVMChiKernel(binnum, 50) );
	}

	OutputFactory fac( shapes, *hfun, binnum, IntPoint(strideX,strideY) );

	OutputArr ra;
	SparseHistOutputArr shra;

	double tc = 0;
	switch (method) {
		case -2:	//integral image
			{
				boost::timer timer;
				ra = fac.ComputeResponses( arr, warr, OutputFactory::METHOD_INTEGRAL );
				shra = fac.ComputeHistogram( arr, warr, OutputFactory::METHOD_INTEGRAL, true );
				tc = timer.elapsed();
			}
			break;
		case -1:	//naive
			{
				boost::timer timer;
				ra = fac.ComputeResponses( arr, warr, OutputFactory::METHOD_NORMAL );
				shra = fac.ComputeHistogram( arr, warr, OutputFactory::METHOD_NORMAL, true );
				tc = timer.elapsed();
			}
			break;
		case 0:		//fast auto
			{
				boost::timer timer;
				ra = fac.ComputeResponses( arr, warr, OutputFactory::METHOD_FAST_AUTO );
				shra = fac.ComputeHistogram( arr, warr, OutputFactory::METHOD_FAST_AUTO, true );
				tc = timer.elapsed();
			}
			break;
		case 1:		//fast edge
			{
				boost::timer timer;
				ra = fac.ComputeResponses( arr, warr, OutputFactory::METHOD_FAST_EDGE );
				shra = fac.ComputeHistogram( arr, warr, OutputFactory::METHOD_FAST_EDGE, true );
				tc = timer.elapsed();
			}
			break;
		case 2:		//fast naive
			{
				boost::timer timer;
				ra = fac.ComputeResponses( arr, warr, OutputFactory::METHOD_FAST_NAIVE );
				shra = fac.ComputeHistogram( arr, warr, OutputFactory::METHOD_FAST_NAIVE, true );
				tc = timer.elapsed();
			}
			break;
	}

	binArrs::Release(arr);
	weightArrs::Release(warr);

	binArr::Release(arrp);
	weightArr::Release(warrp);

	/*for ( size_t y=0; y<rp[0].size().height; ++y ) {
		for ( size_t x=0; x<rp[0].size().width; ++x) {
			cout << rp[0][y][x] << '\t';
		}
		cout << endl;
	}*/

	double winCount = 0;
	{
		struct piBlockHeader {
			int elementSize;
			unsigned int width;
			unsigned int height;
		};

		FILE* out, * hout;
		string	hist_output_file = string(output_file) + "-h";
		if (is_textOutput) {
			out = fopen( output_file, "w"  );
			hout= fopen( hist_output_file.c_str(), "w" );
		} else {
			out = fopen( output_file, "wb" );
			hout= fopen( hist_output_file.c_str(), "wb" );
		}
		
		for ( size_t i=0; i<ra.size(); ++i ) {
			FloatArray2D&	ro = ra(i);
			SparseHistArr&	shro = shra(i);
			piBlockHeader header;
			header.elementSize = sizeof(FloatType);
			header.width       = ro.width;
			header.height	   = ro.height;
			if (is_textOutput) {
				//fprintf(out, "Element size:\t%d\n", header.elementSize);
				//fprintf(out, "Size:\t%u\t%u\n", header.width, header.height );
				FloatType* dptr = ro[0];
				for ( unsigned int y=0; y<header.height; ++y ) {
					for ( unsigned int x=0; x<header.width; ++x) {
						if (x>0)
							fprintf(out, ", ");
						fprintf(out, "%.10lf",  *dptr++);

						fprintf( hout, "(%d,%d):", x, y );
						SparseHistType psh = shro(x,y);
						for ( size_t k=0; k<psh.size() ; ++k ) {
							const SparseHistElement& she = psh[k];
							fprintf( hout, "\t%d:%.5lf", she.idx, she.count );
						}
						fprintf(hout, "\n");
					}
					fprintf(out, "\n");
				}
			} else {
				fwrite( &header, sizeof(header), 1, out );
				fwrite( ro[0], sizeof(FloatType), ro.width*ro.height, out );
				//the histogram output is ABSENT here
			}
			winCount += (double)header.width*header.height;
		}
		fclose(out);
		fclose(hout);
	
	}

	printf( "%.6lf\t%.0lf\n", tc, winCount );

	return 0;
}
