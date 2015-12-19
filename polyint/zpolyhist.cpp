#include "zpolyhist.h"
#include "polyint.h"
#include "sliding_hist.h"
#include "zutils.h"
#include <memory>
#include <list>
#include <algorithm>

namespace zPolyHist {
	using namespace ::std;

	//Shapes =============================

	struct ShapesContent {
		typedef auto_ptr<piPolyGroupArr> ShapesPtr;
		enum ObjectState {
			SC_EDITTING  = 0,
			SC_COMPLETED = 1
		};

		int			state;
		ShapesPtr	shapes;
		piPolygonGroup* tmp_pg;
		list<iPoint>	tmp_vertice;
		int			margin_left;
		int			margin_right;
		int			margin_top;
		int			margin_bottom;
	};

	Shapes::Shapes() {
		_content = new ShapesContent;
		_content->state = ShapesContent::SC_EDITTING;
		_content->shapes = ShapesContent::ShapesPtr( new piPolyGroupArr );
		_content->tmp_pg = new piPolygonGroup;
		_content->margin_top	= 0;
		_content->margin_right	= 0;
		_content->margin_bottom	= 0;
		_content->margin_left	= 0;
	}

	Shapes::Shapes( const Shapes& r ) {
		this->operator =( r );
	}

	Shapes::~Shapes() {
		delete _content;
	}

	const Shapes& Shapes::operator = ( const Shapes& r ) {
		if ( r.IsCompleted() )
			throw Exception_OperateOnNoncompletedShapes();	//Bad style, throw during Construction
		_content = new ShapesContent;
		const ShapesContent* sc2 = (const ShapesContent*)(r._content);
		_content->shapes = ShapesContent::ShapesPtr( new piPolyGroupArr(*(sc2->shapes)) );
		_content->state  = sc2->state;
		return *this;
	}

	void Shapes::AddVertex2Current( IntPoint p ) {
		if ( this->IsCompleted() ) {
			throw Exception_EditCompletedShapes();
		}
		
		list<iPoint>& vertice = _content->tmp_vertice;

		vertice.push_back( iPoint(p.x,p.y) );
	}

	void _shapesEndPolygon( piPolygonGroup*& pg, list<iPoint>& vertice ) {
		array1d<iPoint> parr(vertice.begin(),vertice.end(),array1d<iPoint>::FromIterator());
		piPolygon* poly = new piPolygon(parr);
		pg->push_back(poly);
		vertice.clear();
	}

	void Shapes::NewPolygon() {
		if ( this->IsCompleted() ) {
			throw Exception_EditCompletedShapes();
		}
		
		list<iPoint>& vertice = _content->tmp_vertice;
		piPolygonGroup*& pg = _content->tmp_pg;

		if (!vertice.empty())
			_shapesEndPolygon( pg, vertice );
	}

	void _shapesEndShape( piPolyGroupArr& shapes, piPolygonGroup*& pg, list<iPoint>& vertice, const ShapesContent* _content ) {
		if (!vertice.empty())
			_shapesEndPolygon( pg, vertice );
		//pg->locateOrig();	//Never manually locate at orig ********
		pg->transit( iPoint( _content->margin_left, _content->margin_top ) );
		iPoint& rb = const_cast<iPoint&>(pg->RightBottom);
		rb = pg->RightBottom + iPoint( _content->margin_right, _content->margin_bottom );
		shapes.push_back(pg);
	}

	void Shapes::NewShape() {
		if ( this->IsCompleted() ) {
			throw Exception_EditCompletedShapes();
		}
		
		list<iPoint>& vertice = _content->tmp_vertice;
		piPolygonGroup*& pg = _content->tmp_pg;

		if ( pg->size() || !vertice.empty() ) {
			//EndShape
			_shapesEndShape( *(_content->shapes), pg, vertice, _content );

			//NewShape
			_content->margin_top	= 0;
			_content->margin_right	= 0;	//Currently useless
			_content->margin_bottom	= 0;	//Currently useless
			_content->margin_left	= 0;
			_content->tmp_pg = new piPolygonGroup;
		}
	}

	void Shapes::Complete() {
		if ( this->IsCompleted() ) {
			throw Exception_EditCompletedShapes();
		}
		
		list<iPoint>& vertice = _content->tmp_vertice;
		piPolygonGroup*& pg = _content->tmp_pg;

		if ( pg->size() || !vertice.empty() ) {
			//EndShape
			_shapesEndShape( *(_content->shapes), pg, vertice, _content );
		} else {
			delete _content->tmp_pg;
		}

		_content->state = ShapesContent::SC_COMPLETED;
	}

	bool Shapes::IsCompleted() const {	
		return (_content->state == ShapesContent::SC_COMPLETED);
	}

	size_t Shapes::size() const {
		return _content->shapes->size();
	}
	size_t Shapes::size( size_t idx ) const {
		return (*(_content->shapes))[idx].size();
	}

	IntPoint Shapes::lefttop( size_t idx ) const {
		iPoint p = (*_content->shapes)[idx].LeftTop;
		return IntPoint( p.x, p.y );
	}

	IntPoint Shapes::rightbottom( size_t idx ) const {
		iPoint p = (*_content->shapes)[idx].RightBottom;
		return IntPoint( p.x, p.y );
	}

	void Shapes::SetMargin4CurrentShape( int top, int right, int bottom, int left ) {
		_content->margin_top	= top;
		_content->margin_right	= right;
		_content->margin_bottom	= bottom;
		_content->margin_left	= left;
	}


	//Histogram function =========================
	class FunHull 
		: public piHistFunConcept< HistFunBase::output_type, HistFunBase::input_type, HistFunBase::bin_type> {
	private:
		HistFunBase* _hist_fun;
	public:
		FunHull( HistFunBase& hist_fun ) : _hist_fun(&hist_fun) {}
		output_type operator () (bin_type b, input_type n) {
			return _hist_fun->operator () (b,n);
		}
		output_type zero() { return 0.0; }
		void SetHistFun( HistFunBase& hist_fun ) {
			_hist_fun = &hist_fun;
		}
		~FunHull() {};
	};


	//Output Factory =============================
	
	struct OutputFactoryContent {

		typedef auto_ptr<FunHull>			FunHullPtr;

		typedef piBinMap<UIntArray2D>		BinMapType;
		typedef piBinMap<UIntArray2D,FloatArray2D>		BinWeightedMapType;

		typedef swConf4BetterSlidingWindows WindowConfigType;

		typedef swBetterSlidingWindows		WindowType;
		typedef auto_ptr<WindowType>		WindowPtr;

		typedef piBuffer4HistSlider< FunHull >			SliderBufferType;
		typedef auto_ptr<SliderBufferType>				SliderBufferPtr;

		typedef piHistSlider< FunHull, BinMapType >				SliderType;
		typedef piHistSlider< FunHull, BinWeightedMapType >		WeightedSliderType;

		typedef SliderType::ResultOutputPVec		ResultPool;
		typedef auto_ptr<ResultPool>				ResultPoolPtr;

		typedef piFun4Hist							Fun4Hist;
		typedef piHashHist<piFloat>					HistResultType;
		typedef piBuffer4HistSlider< Fun4Hist >			Slider4HistBufferType;
		typedef piHistSlider< Fun4Hist, BinMapType >	Slider4HistType;
		typedef piHistSlider< Fun4Hist, BinWeightedMapType >	WeightedSlider4HistType;
		typedef Slider4HistType::ResultOutputPVec		HistResultPool;

		size_t			binnum;

		const Shapes*	shapes;
		FunHullPtr		fun;

		WindowPtr		window;

		SliderBufferPtr	slider_buffer;
		SliderBufferPtr	slider_buffer_fake;

		ResultPoolPtr	result_pool;
		OutputArr		outputs;
		SparseHistOutputArr	hist_outputs;
		vector<size_t>	hist_nonzero_count;
	};

	OutputFactory::OutputFactory( const Shapes& _shapes, HistFunBase& _fun, size_t _bin_num, IntPoint _stride ) 
		: _content(new OutputFactoryContent()), HistgramNonzeroCount(_content->hist_nonzero_count) {

		_content->binnum = _bin_num;

		_content->shapes = &_shapes;
		SetHistFun( _fun );

		OutputFactoryContent::WindowConfigType bsConf;
		bsConf.stride = iPoint(_stride.x,_stride.y);
		_content->window = OutputFactoryContent::WindowPtr(
			new OutputFactoryContent::WindowType(*(_content->shapes->_content->shapes), bsConf) );
		
		_content->slider_buffer = OutputFactoryContent::SliderBufferPtr(
			new OutputFactoryContent::SliderBufferType(*(_content->window), _content->binnum ) );
		_content->slider_buffer_fake = OutputFactoryContent::SliderBufferPtr(
			new OutputFactoryContent::SliderBufferType(*(_content->window), _content->binnum ) );
		//_content->slider_buffer->initialize( 5000, _content->binnum );

		_content->outputs = OutputArr::Create( _content->window->win.size() );

		_content->hist_outputs = SparseHistOutputArr::Create( _content->window->win.size() );
		for ( size_t i=0; i<_content->window->win.size(); ++i ) {
			_content->hist_outputs(i).data = NULL;
		}
		_content->hist_nonzero_count.resize( _content->window->win.size() );

	}

	OutputFactory::~OutputFactory() {
		OutputArr::Release( _content->outputs );
		for ( size_t i=0; i<_content->hist_outputs.size(); ++i ) {
			SparseHistArr ha = _content->hist_outputs(i);
			if ( ha.data )
				SparseHistArr::Release( ha );
		}
		SparseHistOutputArr::Release( _content->hist_outputs );
		delete _content;
	}

	void OutputFactory::SetHistFun( HistFunBase& _fun ) {
		_content->fun = OutputFactoryContent::FunHullPtr(
			new FunHull(_fun) );
	}

	void OutputFactory::SetBinNumber( size_t _bin_num ) {
		_content->binnum = _bin_num;
	}

	// Compute responses

	template <class HIST_SLIDER_T>
	const OutputArr _OutputFactory_ComputeResponses
		( OutputFactoryContent* _content, const typename HIST_SLIDER_T::BinMapArr& bmaps, int method ) {
		typedef HIST_SLIDER_T	TheSliderType;

		OutputFactoryContent::SliderBufferType* slider_buffer;
		if ( method == OutputFactory::METHOD_INTEGRAL ) {
			slider_buffer = _content->slider_buffer.get();
		} else {
			slider_buffer = _content->slider_buffer_fake.get();
		}

		TheSliderType hs( bmaps, _content->fun.get(), *(_content->window), slider_buffer );

		int incM = INT_MIN;
		_content->result_pool = OutputFactoryContent::ResultPoolPtr( new OutputFactoryContent::ResultPool );
		OutputFactoryContent::ResultPool& rp = *(_content->result_pool);
		switch (method) {
		case OutputFactory::METHOD_INTEGRAL:	//integral image
			{
				hs.generateIntegralResults(rp, 0);
			}
			break;
		case OutputFactory::METHOD_NORMAL:	//naive
			{
				hs.generateNaiveResults(rp, 0);
			}
			break;
		case OutputFactory::METHOD_FAST_AUTO:		//fast auto
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Auto;
		case OutputFactory::METHOD_FAST_EDGE :		//fast edge
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Edge;
		case OutputFactory::METHOD_FAST_NAIVE:		//fast naive
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Naive;
			_content->window->ChangeIncMethod( incM );
			{
				hs.generateFastResults(rp, 0);
			}
		}

		for ( size_t i=0; i<_content->window->win.size(); ++i) {
			typename TheSliderType::ResultOutput& ro = rp[i];
			FloatType* ptr = (ro[0]);
			_content->outputs(i) = FloatArray2D(ptr, ro.size().width, ro.size().height);
		}

		return _content->outputs;

	}

	const OutputArr OutputFactory::ComputeResponses( const UIntArray2D label_map, int method ) const {
		ArrayOfUIntArray2D  label_maps ( const_cast<UIntArray2D*>(&label_map), 1, 1);
		return ComputeResponses( label_maps, method );
	}

	

	const OutputArr OutputFactory::ComputeResponses( const ArrayOfUIntArray2D label_maps, int method ) const {
		typedef OutputFactoryContent::SliderType::BinMapArr TheBinMapArr;
		TheBinMapArr bmaps(label_maps.size(), 0, TheBinMapArr::WithoutConstruction() );
		for ( size_t i=0; i<label_maps.size(); ++i ) {
			const UIntArray2D& label_map = label_maps(i);
			new (&(bmaps[i])) OutputFactoryContent::BinMapType(label_map, uISize(label_map.width, label_map.height), _content->binnum);
		}

		return _OutputFactory_ComputeResponses<OutputFactoryContent::SliderType>( _content, bmaps, method );
	}

	const OutputArr OutputFactory::ComputeResponses( const UIntArray2D label_map, const FloatArray2D weight_map, int method ) const {
		ArrayOfUIntArray2D  label_maps ( const_cast<UIntArray2D*>(&label_map), 1, 1);
		ArrayOfFloatArray2D weight_maps( const_cast<FloatArray2D*>(&weight_map), 1, 1);
		return ComputeResponses( label_maps, weight_maps, method );
	}

	const OutputArr OutputFactory::ComputeResponses( const ArrayOfUIntArray2D label_maps, const ArrayOfFloatArray2D weight_maps, int method ) const {
		typedef OutputFactoryContent::WeightedSliderType::BinMapArr TheBinMapArr;
		TheBinMapArr bmaps(label_maps.size(), 0, TheBinMapArr::WithoutConstruction() );
		for ( size_t i=0; i<label_maps.size(); ++i ) {
			const UIntArray2D&	label_map  = label_maps(i);
			const FloatArray2D&	weight_map = weight_maps(i);
			new (&(bmaps[i])) OutputFactoryContent::BinWeightedMapType(label_map, weight_map, uISize(label_map.width, label_map.height), _content->binnum);
		}

		return _OutputFactory_ComputeResponses<OutputFactoryContent::WeightedSliderType>( _content, bmaps, method );
	}


	// Compute SparseHistHistogram ------------------------------

	int _cmp_SparseHistElement( const SparseHistElement& a, const SparseHistElement& b ) {
		return (a.idx<b.idx);
	}

	template <class HIST_SLIDER_T>
	const SparseHistOutputArr _OutputFactory_ComputeHistogram
		( OutputFactoryContent* _content, const typename HIST_SLIDER_T::BinMapArr& bmaps, int method, bool need_sorting ) {
		
		typedef HIST_SLIDER_T	TheSliderType;

		OutputFactoryContent::Fun4Hist hfun(_content->binnum);
		TheSliderType hs( bmaps, &hfun, *(_content->window) );
		
		//
		int incM = INT_MIN;
		OutputFactoryContent::HistResultPool	rp;
		switch (method) {
		case OutputFactory::METHOD_INTEGRAL:	//integral image
			{
				hs.generateIntegralResults(rp);
			}
			break;
		case OutputFactory::METHOD_NORMAL:		//naive
			{
				hs.generateNaiveResults(rp);
			}
			break;
		case OutputFactory::METHOD_FAST_AUTO:		//fast auto
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Auto;
		case OutputFactory::METHOD_FAST_EDGE :		//fast edge
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Edge;
		case OutputFactory::METHOD_FAST_NAIVE:		//fast naive
			if (incM == INT_MIN)
				incM = swBetterSlidingWindows::Conf::Inc_Naive;
			_content->window->ChangeIncMethod( incM );
			{
				hs.generateFastResults(rp);
			}
		}

		SparseHistOutputArr ho = _content->hist_outputs;
		for ( size_t i=0; i<_content->window->win.size(); ++i) {
			OutputFactoryContent::Slider4HistType::ResultOutput& ro = rp[i];
			SparseHistArr& ha = ho(i);
			if ( ha.data )
				SparseHistArr::Release( ha );
			ha = SparseHistArr::Create( ro.size().width, ro.size().height );
			size_t& cnt = _content->hist_nonzero_count[i];
			cnt = 0;
			for ( size_t y=0; y<ha.height; ++y ) {
				for ( size_t x=0; x<ha.width; ++x ) {
					SparseHistType& l = ha[y][x];
					OutputFactoryContent::Slider4HistType::ColumnHist& r = ro[y][x];

					l.resize(r.nonzero_size());
					SparseHistType::iterator iter_dst = l.begin();
					for ( OutputFactoryContent::Slider4HistType::ColumnHist::const_iterator iter = r.begin();
						iter != r.end(); ++iter ) {
							iter_dst->idx = iter->idx;
							iter_dst->count = iter->val;
							++iter_dst;
					}
					cnt += l.size();
					if ( need_sorting ) {
						sort( l.begin(),l.end(), _cmp_SparseHistElement );
					}
				}
			}
		}


		return _content->hist_outputs;

	}

	const SparseHistOutputArr OutputFactory::ComputeHistogram( const UIntArray2D label_map, int method, bool need_sorting  ) const {
		ArrayOfUIntArray2D  label_maps ( const_cast<UIntArray2D*>(&label_map), 1, 1);
		return ComputeHistogram( label_maps, method, need_sorting );
	}

	

	const SparseHistOutputArr OutputFactory::ComputeHistogram( const ArrayOfUIntArray2D label_maps, int method, bool need_sorting ) const {
		typedef OutputFactoryContent::SliderType::BinMapArr TheBinMapArr;
		TheBinMapArr bmaps(label_maps.size(), 0, TheBinMapArr::WithoutConstruction() );
		for ( size_t i=0; i<label_maps.size(); ++i ) {
			const UIntArray2D& label_map = label_maps(i);
			new (&(bmaps[i])) OutputFactoryContent::BinMapType(label_map, uISize(label_map.width, label_map.height), _content->binnum);
		}

		return _OutputFactory_ComputeHistogram<OutputFactoryContent::Slider4HistType>( _content, bmaps, method, need_sorting );
	}

	const SparseHistOutputArr OutputFactory::ComputeHistogram( const UIntArray2D label_map, const FloatArray2D weight_map, int method, bool need_sorting ) const {
		ArrayOfUIntArray2D  label_maps ( const_cast<UIntArray2D*>(&label_map), 1, 1);
		ArrayOfFloatArray2D weight_maps( const_cast<FloatArray2D*>(&weight_map), 1, 1);
		return ComputeHistogram( label_maps, weight_maps, method, need_sorting );
	}



	const SparseHistOutputArr OutputFactory::ComputeHistogram( const ArrayOfUIntArray2D label_maps, const ArrayOfFloatArray2D weight_maps, int method, bool need_sorting ) const {

		typedef OutputFactoryContent::WeightedSliderType::BinMapArr TheBinMapArr;
		TheBinMapArr bmaps(label_maps.size(), 0, TheBinMapArr::WithoutConstruction() );
		for ( size_t i=0; i<label_maps.size(); ++i ) {
			const UIntArray2D&	label_map  = label_maps(i);
			const FloatArray2D&	weight_map = weight_maps(i);
			new (&(bmaps[i])) OutputFactoryContent::BinWeightedMapType(label_map, weight_map, uISize(label_map.width, label_map.height), _content->binnum);
		}

		return _OutputFactory_ComputeHistogram<OutputFactoryContent::WeightedSlider4HistType>( _content, bmaps, method, need_sorting );

	}

}
