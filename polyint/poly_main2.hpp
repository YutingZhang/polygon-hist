


#include "sliding_hist.h"
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

class piFun4Test : public piHistFunConcept<piFloat, piFloat, size_t> {
public:
	piFun4Test() {}
	virtual output_type operator () (bin_type b, input_type n) = 0;
	output_type zero() { return 0.0; }
	virtual ~piFun4Test() {};
};

class piFun4TestLinear : public piFun4Test {
private:
	std::vector<piFloat> w;
public:
	piFun4TestLinear(size_t binnum) : w(binnum) {
		for ( size_t i=0; i<binnum; ++i ) {
			//w[i] = (double)rand()/RAND_MAX;
			w[i]=i+1;
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
	std::vector<piFloat> w;
public:
	piFun4TestChi(size_t binnum) : w(binnum) {
		for ( size_t i=0; i<binnum; ++i )
			w[i] = ((double)rand()/RAND_MAX)*10.0;
	}
	virtual output_type operator () (bin_type b, input_type n) {
		piFloat a=(n+w[b]);
		return a*a/(n-w[b]);
	}
	virtual ~piFun4TestChi() {};
};

class piFun4TestEntropy : public piFun4Test {
private:
	std::vector<piFloat> w;
public:
	piFun4TestEntropy(){}
	virtual output_type operator () (bin_type b, input_type n) {
		return (n+1)*log(n+1);
	}
	virtual ~piFun4TestEntropy() {};
};

class piFun4TestSVMChiKernel : public piFun4Test {
private:
	std::vector<piFloat> w;
	array2d<piFloat>	m;
public:
	piFun4TestSVMChiKernel(size_t binnum, size_t sv = 50) : w(sv), m(uISize(binnum, sv)) {
		for ( size_t k=0; k<sv; ++k ) {
			w[k] = ((piFloat)rand()/RAND_MAX);
			for ( size_t i=0; i<binnum; ++i )
				m[k][i] = ((piFloat)rand()/RAND_MAX);
		}
	}
	virtual output_type operator () (bin_type b, input_type n) {
		piFloat s = 0.0;
		for ( size_t k=0; k<w.size(); ++k ) {
			piFloat a=n+m[k][b];
			s += w[k]*(a*a/(n-m[k][b]));
		}
		return s;
	}
	virtual ~piFun4TestSVMChiKernel() {};
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
	bool	is_labelMap = false;
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

	piPolyGroupArr	shapes;
	{
		//Load polygons
		for ( list<string>::iterator iter = poly_files.begin();
			iter != poly_files.end(); ++iter ) {
				ifstream in(iter->c_str());
				if (!in) {
					cerr << "Can't open " << *iter << endl;
					continue;
				}
				iPoint p;
				piPolygonGroup* pg = new piPolygonGroup;
				list<iPoint> vertice;
				while ( in >> p.x >> p.y ) {
					if ( p.x== -1 && p.y== -1 ) {		//use (-1,-1) as terminator
						array1d<iPoint> parr(vertice.begin(),vertice.end(),array1d<iPoint>::FromIterator());
						piPolygon* poly = new piPolygon(parr);
						pg->push_back(poly);
						vertice.clear();
					} else {
						vertice.push_back(p*zoom);
					}
				}
				if (!vertice.empty()) {
					array1d<iPoint> parr(vertice.begin(),vertice.end(),array1d<iPoint>::FromIterator());
					piPolygon* poly = new piPolygon(parr);
					pg->push_back(poly);
				}
				pg->locateOrig(); //Better or not?
				shapes.push_back(pg);
		}
	}


	if ( !shapes.size() ) {
		cerr << "No polygroup file." << endl;
		exit(1);
	}


	
	int binnum;
	typedef array2d<unsigned> binArr;
	typedef std::auto_ptr<binArr> binArrPtr;
	binArrPtr arrp;
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
		arrp = binArrPtr( new binArr( uISize( imgSize.width, imgSize.height ) ) );
		for ( int y=0; y<imgSize.height; ++y ) {
			for ( int x=0; x<imgSize.width; ++x ) {
				CvScalar c = cvGet2D( imgc, y, x );
				(*arrp)[y][x] = ((int)c.val[0]/colorDepthFrac1) * true_base2 +
					(int)c.val[2]/colorDepthFrac2;
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

		arrp = binArrPtr( new binArr( uISize( content[0].size(), content.size() ) ) );
		for ( int y=0; y<(int)arrp->size().height; ++y ) {
			for ( int x=0; x<(int)arrp->size().width; ++x ) {
				(*arrp)[y][x] = content[y][x];
			}
		}
	}

	cerr << "Bin num: " << binnum << endl;

	binArr& arr = *arrp;
	swConf4BetterSlidingWindows bsConf;
	bsConf.stride = iPoint(strideX,strideY);
	swBetterSlidingWindows windows(shapes, bsConf);
	typedef piFun4Test FunType;
	typedef piBinMap<binArr> binmap_t;
	typedef piHistSlider< FunType, binmap_t > sliding_t;
	binmap_t bmap(arr, arr.size(), binnum);

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

	sliding_t hs(bmap, hfun.get(), windows);

	double tc = 0;
	int incM = INT_MIN;
	sliding_t::ResultOutputPVec rp;
	switch (method) {
		case -2:	//integral image
			{
				boost::timer timer;
				hs.generateIntegralResults(rp, 0);
				tc = timer.elapsed();
			}
			break;
		case -1:	//naive
			{
				boost::timer timer;
				hs.generateNaiveResults(rp, 0);
				tc = timer.elapsed();
			}
			break;
		case 0:		//fast auto
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Auto;
		case 1:		//fast edge
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Edge;
		case 2:		//fast naive
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Naive;
			windows.ChangeIncMethod( incM );
			{
				boost::timer timer;
				hs.generateFastResults(rp, 0);
				tc = timer.elapsed();
			}
	}

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

		FILE* out;
		if (is_textOutput)
			out = fopen( output_file, "w" );
		else
			out = fopen( output_file, "wb" );
		
		for ( size_t i=0; i<rp.size(); ++i ) {
			sliding_t::ResultOutput& ro = rp[i];
			piBlockHeader header;
			header.elementSize = sizeof(sliding_t::result_t);
			header.width       = ro.size().width;
			header.height	   = ro.size().height;
			if (is_textOutput) {
				//fprintf(out, "Element size:\t%d\n", header.elementSize);
				//fprintf(out, "Size:\t%u\t%u\n", header.width, header.height );
				sliding_t::result_t* dptr = ro[0];
				for ( unsigned int y=0; y<header.height; ++y ) {
					for ( unsigned int x=0; x<header.width; ++x) {
						if (x>0)
							fprintf(out, ", ");
						fprintf(out, "%.10lf",  *dptr++);
					}
					fprintf(out, "\n");
				}
			} else {
				fwrite( &header, sizeof(header), 1, out );
				fwrite( ro[0], sizeof(sliding_t::result_t), ro.size().width*ro.size().height, out );
			}
			winCount += (double)header.width*header.height;
		}
		fclose(out);
	
	}

	printf( "%.6lf\t%.0lf\n", tc, winCount );

	return 0;
}

