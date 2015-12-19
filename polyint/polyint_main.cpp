/*
 * polyint_main.cpp
 *
 *  Created on: 2010-10-21
 *      Author: zhangyuting
 */

#include<cstdlib>
#include<opencv/cxcore.h>
#include<opencv/cv.h>
#include<opencv/highgui.h>
#include<list>
#include<vector>
#include<cstdio>
#include<iostream>
#include<fstream>
#include"zutils.h"
#include"polyint.h"
#include"hashspecial.h"
#include"sliding_hist.h"
#include<boost/rational.hpp>
#include<boost/progress.hpp>
using namespace std;


int test_consistency( swBetterSlidingWindows& sw );

int test_consistency_fast( swBetterSlidingWindows& );

int test_consistency_integral( swBetterSlidingWindows& sw );

int test_fast( swBetterSlidingWindows& sw );

int test_integral( swBetterSlidingWindows& sw );

int test_naive( swBetterSlidingWindows& sw );

const char def_imgfile[]= "test.jpg";
const char* imgfile;
int use_fast = 1;

int main_time(const int argc, const char* argv[]);
int main_test(const int argc, const char* argv[]);

int main(const int argc, const char* argv[]) {
	return main_time(argc, argv);
}

int main_test(const int argc, const char* argv[]) {

	if ( argc>1 ) {
		sscanf(argv[1], "%d", &use_fast);
	};

	if ( argc>2 )
		imgfile = argv[2];
	else
		imgfile = def_imgfile;

#if 0

	{
		void basic_test();
		basic_test();
	}

#else
	{
		void fun();
		fun();
	}

#endif


	//printf("Hello world!\n");
	return 0;
}

struct byte3 {
	unsigned char a[3];
};

void fun() {
	piPolyGroupArr shapes;
	piPolygonGroup* pg = new piPolygonGroup();

	//-------------------------------------------
	vector<iPoint> vva;

	/*vva.push_back(iPoint(20,0));
	vva.push_back(iPoint(40,0));
	vva.push_back(iPoint(50,10));
	vva.push_back(iPoint(50,20));
	vva.push_back(iPoint(40,30));
	vva.push_back(iPoint(30,10));
	vva.push_back(iPoint(20,30));
	vva.push_back(iPoint(30,40));
	vva.push_back(iPoint(0,40));*/


	vva.push_back(iPoint(0,0));
	vva.push_back(iPoint(30,10));
	vva.push_back(iPoint(25,30));
	vva.push_back(iPoint(10,20));


	piPolygon::VertexArr va( vva.begin(), vva.end(), piPolygon::VertexArr::FromIterator() );
	pg->push_back(new piPolygon(va));
	pg->locateOrig();
	shapes.push_back(pg);

	//piSlidingWindows windows(shapes);
	swBetterSlidingWindows windows(shapes);
	//--------------------------------------------

	//sliding_t hs(arr, hfun, windows);
	
	windows.ChangeIncMethod( swBetterSlidingWindows::Conf::Inc_Auto );

	if ( use_fast==1 )
		test_fast(windows);
	else if ( use_fast==0 )
		test_naive(windows);
	else if ( use_fast==-1 )
		test_consistency(windows);
	else if ( use_fast==2 )
		test_consistency_fast(windows);
	else if ( use_fast==3 )
		test_integral(windows);
	else if ( use_fast==4 )
		test_consistency_integral(windows);

}

int test_consistency_fast( swBetterSlidingWindows& sw ) {

	typedef piFun4TestSum FunType;
	typedef array2d<unsigned char> binmap_o;
	typedef piBinMap<binmap_o> binmap_t;
	typedef piHistSlider< FunType, binmap_t > sliding_t;

	FunType hfun;

	IplImage* img = cvLoadImage(imgfile, CV_LOAD_IMAGE_GRAYSCALE);
	if(!img) {
		cerr << "Can't open " << imgfile << endl;
		exit(1);
	}

	CvSize imgsize = cvGetSize(img);

	binmap_o  arr1(img->imageData, uISize(imgsize.width,imgsize.height), img->widthStep );
	binmap_t arr(arr1, arr1.size(), 256);
	swBetterSlidingWindows sw1(sw);
	sliding_t hs(arr, &hfun, sw);
	sliding_t hs1(arr, &hfun, sw1);

	sliding_t::ResultOutputPVec rp;
	hs1.generateFastResults(rp);

	size_t err_n = 0;
	while ( hs.stepSlide() ) {
		cout << '(' << hs.current().x << ',' << hs.current().y << ")\t";
		for ( size_t i=0; i<hs.Results.size(); ++i ) {
			if (!hs.ResultMask[i]) {
				cout << "NONE\t";
				continue;
			}
			piFloat rt = rp[i][hs.current().y][hs.current().x];
			if ( hs.Results[i] == rt ) {
				cout << "OK\t" << hs.Results[i] << '\t' << rt << "\t";
			} else {
				cout << "WRONG: " << hs.Results[i] << '\t' << rt << "\t";
				++err_n;
			}
		}
		cout << endl;
	}

	cvReleaseImage(&img);

	return err_n;
}

int test_consistency( swBetterSlidingWindows& sw ) {

	typedef piFun4Hist FunType;
	typedef array2d<unsigned char> binmap_o;
	typedef piBinMap<binmap_o> binmap_t;
	typedef piHistSlider< FunType, binmap_t > sliding_t;

	FunType hfun(256);

	IplImage* img = cvLoadImage(imgfile, CV_LOAD_IMAGE_GRAYSCALE);
	if(!img) {
		cerr << "Can't open " << imgfile << endl;
		exit(1);
	}

	CvSize imgsize = cvGetSize(img);

	binmap_o  arr1(img->imageData, uISize(imgsize.width,imgsize.height), img->widthStep );
	binmap_t arr(arr1, arr1.size(), 256);
	sliding_t hs(arr, &hfun, sw);

	size_t err_n = 0;
	while ( hs.stepSlide() ) {
		cout << '(' << hs.current().x << ',' << hs.current().y << ")\t";
		hs.histAt(hs.current());
		for ( size_t i=0; i<hs.Results.size(); ++i ) {
			if (!hs.ResultMask[i]) {
				cout << "NONE\t";
				continue;
			}
			if ( hs.Results[i].is_equal2arr(hs.HistFixed[i]) ) {
				cout << "OK\t";
			} else {
				cout << "WRONG" << "\t";
				for ( size_t j=0; j<arr.size().width; ++j) {
					cout <<"@\t"<< hs.Results[i][j] << '\t' << hs.HistFixed[i][j] << '\t'
						<< hs.Results[i][j]-hs.HistFixed[i][j] << endl;
				}
				++err_n;
			}
		}
		cout << endl;
	}

	cvReleaseImage(&img);

	return err_n;
}

int test_fast( swBetterSlidingWindows& sw ) {
	typedef piFun4TestSum FunType;
	typedef array2d<unsigned char> binmap_o;
	typedef piBinMap<binmap_o> binmap_t;
	typedef piHistSlider< FunType, binmap_t > sliding_t;
	typedef piBuffer4HistSlider< FunType > buffer_t;


	IplImage* img = cvLoadImage(imgfile, CV_LOAD_IMAGE_GRAYSCALE);
	if(!img) {
		cerr << "Can't open " << imgfile << endl;
		exit(1);
	}

	CvSize imgsize = cvGetSize(img);

	binmap_o arr1(img->imageData, uISize(imgsize.width,imgsize.height), img->widthStep );
	FunType hfun;
	buffer_t buffer(sw);
	binmap_t arr(arr1, arr1.size(), 256);
#if 0
	buffer.initialize(2000, 256);
	sliding_t hs(arr, hfun, sw, &buffer);
#else
	sliding_t hs(arr, &hfun, sw);
#endif
	

	int n=0;
#if 0
	{
		boost::progress_timer timer;
		while ( hs.stepSlide() ) {
			++n;
			//cout << hs.Results[0] << endl;
		}
	}
#else
	{
		boost::progress_timer timer;
		sliding_t::ResultOutputPVec rp;
		hs.generateFastResults(rp);
	}
#endif
	
	cvReleaseImage(&img);

	return n;
}

int test_naive( swBetterSlidingWindows& sw ) {
	typedef piFun4TestSum FunType;
	typedef array2d<unsigned char> binmap_o;
	typedef piBinMap<binmap_o> binmap_t;
	typedef piHistSlider< FunType, binmap_t > sliding_t;

	FunType hfun;

	IplImage* img = cvLoadImage(imgfile, CV_LOAD_IMAGE_GRAYSCALE);
	if(!img) {
		cerr << "Can't open " << imgfile << endl;
		exit(1);
	}

	CvSize imgsize = cvGetSize(img);

	binmap_o  arr1(img->imageData, uISize(imgsize.width,imgsize.height), img->widthStep );
	binmap_t arr(arr1, arr1.size(), 256);
	

	sliding_t hs(arr, &hfun, sw);

	int n=0;
	{
		boost::progress_timer timer;
		//cout << hs.SlidingRange[0].width << ',' << hs.SlidingRange[0].height << endl;
		for (int y=hs.SlidingRange[0].height; y>=0; --y) {
			for (int x=hs.SlidingRange[0].width; x>=0; --x) {
				hs.funAt(x,y);
				cout << hs.ResultsFixed[0] << endl;
				++n;
			}
		}
	}

	cvReleaseImage(&img);

	return n;
}

int test_integral( swBetterSlidingWindows& sw ) {
	typedef piFun4TestSum FunType;
	typedef array2d<unsigned char> binmap_o;
	typedef piBinMap<binmap_o> binmap_t;
	typedef piHistSlider< FunType, binmap_t > sliding_t;
	typedef piBuffer4HistSlider< FunType > buffer_t;

	IplImage* img = cvLoadImage(imgfile, CV_LOAD_IMAGE_GRAYSCALE);
	if(!img) {
		cerr << "Can't open " << imgfile << endl;
		exit(1);
	}

	CvSize imgsize = cvGetSize(img);

	binmap_o arr1(img->imageData, uISize(imgsize.width,imgsize.height), img->widthStep );
	FunType hfun;
	binmap_t arr(arr1, arr1.size(), 256);
	sliding_t hs(arr, &hfun, sw);
	
	size_t n = 0;
	{
		boost::progress_timer timer;
		hs.InitializeIntegralImage();

		while ( hs.IntegralSlide() ) {
			++n;
		}
	}

	cvReleaseImage(&img);

	return n;
}

int test_consistency_integral( swBetterSlidingWindows& sw ) {
	typedef piFun4Hist FunType;
	typedef array2d<unsigned char> binmap_o;
	typedef piBinMap<binmap_o> binmap_t;
	typedef piHistSlider< FunType, binmap_t > sliding_t;
	typedef piBuffer4HistSlider< FunType > buffer_t;


	IplImage* img = cvLoadImage(imgfile, CV_LOAD_IMAGE_GRAYSCALE);
	if(!img) {
		cerr << "Can't open " << imgfile << endl;
		exit(1);
	}

	CvSize imgsize = cvGetSize(img);

	binmap_o arr1(img->imageData, uISize(imgsize.width,imgsize.height), img->widthStep );
	FunType hfun(8);
	
	for ( size_t y=0; y<arr1.size().height; ++y ) {
		for ( size_t x=0; x<arr1.size().width; ++x ) {
			arr1[y][x] = arr1[y][x] % 8;
		}
	}
	binmap_t arr(arr1, arr1.size(), 8);
	sliding_t hs(arr, &hfun, sw);
	swBetterSlidingWindows& sw1(sw);
	sliding_t hs1(arr, &hfun, sw1);
	
	hs.InitializeIntegralImage();

	size_t err_n = 0;
	while ( hs.IntegralSlide() ) {
		hs1.stepSlide();
		cout << '(' << hs.current().x << ',' << hs.current().y << ")\t";
		for ( size_t i=0; i<hs.Results.size(); ++i ) {
			if (!hs.ResultMask[i]) {
				cout << "NONE\t";
				continue;
			}
			if ( hs1.Results[i].is_equal2arr(hs.HistFixed[i]) ) {
				cout << "OK\t";
				for ( size_t j=0; j<hs.HistFixed[i].size(); ++j ) {
					cout << hs.HistFixed[i][j] << '\t';
				}
			} else {
				cout << "WRONG" << "\t";
				/*for ( size_t j=0; j<arr.size().width; ++j) {
					cout <<"@\t"<< hs1.Results[i][j] << '\t' << hs.HistFixed[i][j] << '\t'
						<< hs1.Results[i][j]-hs.HistFixed[i][j] << endl;
				}*/
				++err_n;
			}
		}
		cout << endl;
	}

	cvReleaseImage(&img);

	return err_n;
}

void basic_test() {
	{
		cout << "piPolygon test --------------" << endl;

		vector<iPoint> vva;

		
		vva.push_back(iPoint(20,0));
		vva.push_back(iPoint(40,0));
		vva.push_back(iPoint(50,10));
		vva.push_back(iPoint(50,20));
		vva.push_back(iPoint(40,30));
		vva.push_back(iPoint(30,10));
		vva.push_back(iPoint(20,30));
		vva.push_back(iPoint(30,40));
		vva.push_back(iPoint(0,40));

		/*
		vva.push_back(iPoint(0,0));
		vva.push_back(iPoint(30,0));
		vva.push_back(iPoint(30,20));
		vva.push_back(iPoint(0,20));
		*/
		
		piPolygon::VertexArr va( vva.begin(), vva.end(), piPolygon::VertexArr::FromIterator() );


		piPolygon poly(va);
		cout << "LeftTop: " << poly.LeftTop.x << '\t' << poly.LeftTop.y << endl;
		cout << "Size: " << poly.Size.width << '\t' << poly.Size.height << endl;
		
		for (size_t i=0;i<poly.RightPPR.size();++i) {
			cout << poly.RightPPR[i] << "\t" << poly.TopPPR[i] << "\t" << poly.Boundary[i].Slope <<
					"\t(" << poly.NormEdges[i].sp.x << ',' << poly.NormEdges[i].sp.y <<
					")+("  << poly.NormEdges[i].vec.x << ',' << poly.NormEdges[i].vec.y << ")" << endl;
		}
		cout << "--------------" << endl;
		char _polyMark[] = {'.','X','O'};
		char* polyMark = _polyMark + 1;
		for (size_t y=0;y<poly.Size.height;++y) {
			for (size_t x=0;x<poly.Size.width;++x) {
				cout << polyMark[poly.GetMask()[y][x]];
			}
			cout << endl ;
		}
		
		cout << "////--------------" << endl;
	}

	{
		cout << "piUnitRasterizedLineSeg test --------------" << endl;
		piUnitRasterizedLineSeg urls(iPoint(30,10));
		cout << "Vector: " << urls.Vec.x << "\t" << urls.Vec.y << endl;
		cout << "Slope: " << urls.Slope << endl;
		for (size_t i=0;i<urls.PixelCount;++i) {
			cout    << urls.Pixels[i].x << "\t"
					<< urls.Pixels[i].y << "\t"
					<< urls.LeftPortions[i] << "\t" << urls.RightPortions[i] << "\t"
					<< urls.TopPortions[i] << "\t" << urls.BottomPortions[i]
					<< endl;
		}
		cout << "--------------" << endl;
		for (size_t i=0;i<urls.InnerYs.size();++i) {
			cout    << urls.InnerYs[i] << "\t"
					<< urls.TopYs[i] << endl;
		}
		cout << "////--------------" << endl;
	}

	{
		cout << "piRasterizedLineSeg test --------------" << endl;
		piRasterizedLineSeg rls(lineseg_t(iPoint(3,2),iPoint(6,8)));

		cout << "Pixel count: " << rls.Pixels.size() << endl;

		for (size_t i=0;i<rls.PixelCount;++i) {
			cout    << rls.Pixels[i].x << "\t"
					<< rls.Pixels[i].y << "\t"
					<< rls.LeftPortions[i] << "\t" << rls.RightPortions[i] << "\t"
					<< rls.TopPortions[i] << "\t" << rls.BottomPortions[i]
					<< endl;
		}
		cout << "--------------" << endl;
		for (size_t i=0;i<rls.InnerYs.size();++i) {
			cout    << rls.InnerYs[i] << "\t"
					<< rls.TopYs[i] << endl;
		}

		cout << "////--------------" << endl;

		cout << "piPixelIntersection test --------------" << endl;
		piRasterizedLineSeg rls_y(lineseg_t(iPoint(3,3),iPoint(6,8)));
		piRasterizedLineSeg rls_x(lineseg_t(iPoint(4,2),iPoint(6,8)));
		piRasterizedLineSeg rls_xy(lineseg_t(iPoint(4,3),iPoint(6,8)));
		piPixelIntersection pinter( rls.Pixels, rls_xy.Pixels );
		cout << "Intersection: " << endl;
		for (size_t i=0;i<pinter.Intersection.size();++i) {
			cout    << pinter.Intersection[i].x << "\t"
					<< pinter.Intersection[i].y << endl;
		}
		cout << "Residual1: " << endl;
		for (size_t i=0;i<pinter.Residual1.size();++i) {
			cout    << pinter.Residual1[i].x << "\t"
					<< pinter.Residual1[i].y << endl;
		}
		cout << "Residual2: " << endl;
		for (size_t i=0;i<pinter.Residual2.size();++i) {
			cout    << pinter.Residual2[i].x << "\t"
					<< pinter.Residual2[i].y << endl;
		}
		cout << "////--------------" << endl;
	}

	{
		cout << "piHashSpecial --------------" << endl;
		typedef piZeroScheme::Epsilon<piFloat> FloatZeroScheme;
		typedef piHashSpecial<piFloat,FloatZeroScheme> floatHashSpecial;
		floatHashSpecial hs(100, 10 );
		hs.WriteTable(31,5233);
		hs.WriteTable(33,5132);
		hs.clear();
		hs.WriteTable(4,1246);
		hs.WriteTable(4,0);
		hs.WriteTable(2,6);
		hs.WriteTable(1,13);
		hs.WriteTable(7,73);
		hs.WriteTable(98,552);
		hs.WriteTable(7,hs[7]+0.4f);
		hs.WriteTable(5,0.546f);
		hs.WriteTable(2,0);
		hs.WriteTable(11,1);
		hs.WriteTable(12,2);
		hs.WriteTable(13,3);
		hs.WriteTable(14,4);
		hs.WriteTable(15,5);
		hs.WriteTable(13,0);
		hs.WriteTable(16,6);
		hs.WriteTable(17,7);
		hs.WriteTable(17,0);
		hs.WriteTable(18,8);
		cout << "Nonzero Size: " << hs.nonzero_size() << endl;
		for (floatHashSpecial::iterator iter=hs.begin();
			iter!=hs.end(); ++iter) {
				cout << iter->idx << '\t' << iter->val << endl;
		}
		cout << "--------------" << endl;
		for (size_t i=0; i<hs.size(); ++i) {
			cout << i << '\t' << hs.ReadTable(i) << endl;
		}

		cout << "////--------------" << endl;
	}

	{
		cout << "zArray2D test --------------" << endl;
		array1d<int> arr1( 4*4 );
		for ( size_t i=0; i<arr1.size(); ++i) {
			arr1[i] = (int)i+1;
		}

		const array2d<int> arr2( arr1.begin(), uISize(3,4), sizeof(int)*4 );
		for ( size_t y=0; y<arr2.size().height; ++y) {
			for ( size_t x=0; x<arr2.size().width; ++x) {
				cout << arr2[y][x] << '(' << arr2(x,y) << ")\t";
			}
			cout << endl;
		}

		cout << "////--------------" << endl;
	}
}
