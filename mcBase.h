

#pragma once

//  Main file in the simulation library
//  All the base classes and template algorithms
//  See chapters 6, 7, 12 and 14

#include "AAD.h"

#include <vector>
#include <memory>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>

using namespace std;

#include "matrix.h"
#include "ThreadPool.h"

using Time = double;
extern Time systemTime;

//Data to be simulated
struct SampleStructure
{
    bool            needNumeraire = true;
    vector<Time>    forwardMats;
    vector<Time>    discountMats;

    struct RateDef
    {
        Time    start;
        Time    end;
        string  curve;

        RateDef(const Time s, const Time e, const string& c) :
            start(s), end(e), curve(c) {};
    };

    vector<RateDef> liborDefs;
};

//--------------------------------------------------------------------------------------------
template <class T>
struct ValuesForSampleStruct
{
    T           needNumeraire;
    vector<T>   forwards;
    vector<T>   discounts;
    vector<T>   libors;

    //  Allocate given SampleDef
    void allocate(const SampleStructure& data)
    {
        forwards.resize(data.forwardMats.size());
        discounts.resize(data.discountMats.size());
        libors.resize(data.liborDefs.size());
    }

    void initialize()
    {
        needNumeraire = T(1.0);
		fill(forwards.begin(), forwards.end(), T(100.0));
        fill(discounts.begin(), discounts.end(), T(1.0));
        fill(libors.begin(), libors.end(), T(0.0));
    }
};
//--------------------------------------------------------------------------------------------
template <class T>
using Scenario = vector<ValuesForSampleStruct<T>>;

template <class T>
inline void allocatePath(const vector<SampleStructure>& getDeflineRef, Scenario<T>& path)
{
    path.resize(getDeflineRef.size());
    for (size_t i = 0; i < getDeflineRef.size(); ++i)
    {
        path[i].allocate(getDeflineRef[i]);
    }
}
//--------------------------------------------------------------------------------------------
template <class T>
inline void initializePath(Scenario<T>& path)
{
    for (auto& scen : path) scen.initialize();
}
//--------------------------------------------------------------------------------------------
template <class T>
class Product
{
public:
    virtual const vector<Time>& getTimelineRef() const = 0;
    virtual const vector<SampleStructure>& getDeflineRef() const = 0;
    virtual const vector<string>& getPayoffLabelsRef() const = 0;

    virtual void computePayoffs(const Scenario<T>&  path,     
                         vector<T>&       computePayoffs) const = 0;

	virtual unique_ptr<Product<T>> getClonePtr() const = 0;
    virtual ~Product() {}
};
//--------------------------------------------------------------------------------------------
template <class T>
class Model
{
public:
    virtual void allocate(const vector<Time>&         prdTimeline, 
                          const vector<SampleStructure>&    prdDefline) = 0;
    virtual void init( const vector<Time>&         prdTimeline, 
                       const vector<SampleStructure>&    prdDefline) = 0;

    virtual size_t simDim() const = 0;

    virtual void generatePath(const vector<double>& gaussVec, 
                              Scenario<T>&              path) const = 0;

	virtual unique_ptr<Model<T>> getClonePtr() const = 0;
    virtual ~Model() {}

    virtual const vector<T*>& getParameters() = 0;
    virtual const vector<string>& getParameterLabels() const = 0;
    size_t getNumParams() const { return const_cast<Model*>(this)->getParameters().size();}

	void putParametersOnTape() { putParametersOnTapeT<T>(); }// only valid for T = Number

private:
    template<class U> 
    void putParametersOnTapeT(){}

    //  If T = Number : put on tape
    template <>
    void putParametersOnTapeT<Number>() { for (Number* param : getParameters()) param->putOnTape();}
};
//-------------------------------------------------------------------------------------------
class RNG
{
public:
    
    //  Initialise with dimension simDim
    virtual void init(const size_t simDim) = 0;

    //  Compute the next vector[simDim] of independent Uniforms or Gaussians
    //  The vector is filled by the function and must be pre-allocated
	virtual void nextU(vector<double>& uVec) = 0;
	virtual void nextG(vector<double>& gaussVec) = 0;

    virtual unique_ptr<RNG> getClonePtr() const = 0;

    virtual ~RNG() {}

    //  Skip ahead
    virtual void skipTo(const unsigned b) = 0;
};

//  Template algorithms
//  ===================

//  Serial valuation, chapter 6

//	MC simulator: free function that conducts simulations 
//      and returns a matrix (as vector of vectors) of payoffs 
//          (0..nPath-1 , 0..nPay-1) 
inline vector<vector<double>> mcSimul(
    const Product<double>&      prd,
    const Model<double>&        mdl,
    const RNG&                  rng,			            
    const size_t                nPath)                      
{
    //  Work with copies of the model and RNG
    //      which are modified when we set up the simulation
    //  Copies are OK at high level
    auto cMdl = mdl.getClonePtr();
    auto cRng = rng.getClonePtr();

    //	Allocate results
    const size_t nPay = prd.getPayoffLabelsRef().size();
    vector<vector<double>> results(nPath, vector<double>(nPay));
    //  Init the simulation timeline
    cMdl->allocate(prd.getTimelineRef(), prd.getDeflineRef());
    cMdl->init(prd.getTimelineRef(), prd.getDeflineRef());              
    //  Init the RNG
    cRng->init(cMdl->simDim());                        
    //  Allocate Gaussian vector
    vector<double> gaussVec(cMdl->simDim());           
    //  Allocate path
    Scenario<double> path;
    allocatePath(prd.getDeflineRef(), path);
    initializePath(path);

    //	Iterate through paths	
    for (size_t i = 0; i<nPath; i++)
    {
        //  Next Gaussian vector, dimension D
        cRng->nextG(gaussVec);                        
        //  Generate path, consume Gaussian vector
        cMdl->generatePath(gaussVec, path);     
        //	Compute result
        prd.computePayoffs(path, results[i]);
    }

    return results;	//	C++11: move
}

//  Parallel valuation, chapter 7

#define BATCHSIZE 64
//	Parallel equivalent of mcSimul()
inline vector<vector<double>> mcParallelSimul(
    const Product<double>&      prd,
    const Model<double>&        mdl,
    const RNG&                  rng,
    const size_t                nPath)
{
    auto cMdl = mdl.getClonePtr();

    const size_t nPay = prd.getPayoffLabelsRef().size();
    vector<vector<double>> results(nPath, vector<double>(nPay));

    cMdl->allocate(prd.getTimelineRef(), prd.getDeflineRef());
    cMdl->init(prd.getTimelineRef(), prd.getDeflineRef());

    //  Allocate space for Gaussian vectors and paths, 
    //      one for each thread
    ThreadPool *pool = ThreadPool::getInstance();
    const size_t nThread = pool->getTotalNumerOfThreads();
    vector<vector<double>> gaussVecs(nThread+1);    //  +1 for main
    vector<Scenario<double>> paths(nThread+1);
    for (auto& vec : gaussVecs) vec.resize(cMdl->simDim());
    for (auto& path : paths)
    {
        allocatePath(prd.getDeflineRef(), path);
        initializePath(path);
    }
    
    //  One RNG per thread
    vector<unique_ptr<RNG>> rngs(nThread + 1);
    for (auto& random : rngs)
    {
        random = rng.getClonePtr();
        random->init(cMdl->simDim());
    }

    //  Reserve memory for futures
    vector<TaskHandle> futures;
    futures.reserve(nPath / BATCHSIZE + 1); 

    //  Start
    //  Same as mcSimul() except we send tasks to the pool 
    //  instead of executing them

    size_t firstPath = 0;
    size_t pathsLeft = nPath;
    while (pathsLeft > 0)
    {
        size_t pathsInTask = min<size_t>(pathsLeft, BATCHSIZE);

        futures.push_back( pool->spawnTask ( [&, firstPath, pathsInTask]()
        {
            //  Inside the parallel task, 
            //      pick the right pre-allocated vectors
            const size_t getThreadNumer = pool->getThreadNumer();
            vector<double>& gaussVec = gaussVecs[getThreadNumer];
            Scenario<double>& path = paths[getThreadNumer];

            //  Get a RNG and position it correctly
            auto& random = rngs[getThreadNumer];
            random->skipTo(firstPath);

            //  And conduct the simulations, exactly same as sequential
            for (size_t i = 0; i < pathsInTask; i++)
            {
                //  Next Gaussian vector, dimension D
                random->nextG(gaussVec);
                //  Path
                cMdl->generatePath(gaussVec, path);       
                //  Payoff
                prd.computePayoffs(path, results[firstPath + i]);
            }

            //  Remember tasks must return bool
            return true;
        }));

        pathsLeft -= pathsInTask;
        firstPath += pathsInTask;
    }

    //  Wait and help
    for (auto& future : futures) pool->activeWait(future);

    return results;	//	C++11: move
}

//  AAD instrumentation of mcSimul(), chapter 12

//  returns the following results:
struct AADSimulResults
{
    AADSimulResults(const size_t nPath, const size_t nPay, const size_t nParam) :
        computePayoffs(nPath, vector<double>(nPay)),
        aggregated(nPath),
        risks(nParam)
    {}

    //  matrix(0..nPath - 1, 0..nPay - 1) of payoffs, same as mcSimul()
    vector<vector<double>>  computePayoffs;

    //  vector(0..nPath) of aggregated payoffs
    vector<double>          aggregated;

    //  vector(0..nParam - 1) of risk sensitivities
    //  of aggregated payoff, averaged over paths
    vector<double>          risks;
};

//  Default aggregator = 1st payoff = payoff[0]
const auto defaultAggregator = [](const vector<Number>& v) {return v[0]; };

template<class F = decltype(defaultAggregator)>
inline AADSimulResults
mcSimulAAD(
    const Product<Number>&  prd,
    const Model<Number>&    mdl,
    const RNG& rng,
    const size_t            nPath,
    const F&                aggFun = defaultAggregator)
{
    //  Work with copies of the model and RNG
    //      which are modified when we set up the simulation
    //  Copies are OK at high level
    auto cMdl = mdl.getClonePtr();
    auto cRng = rng.getClonePtr();

    //  Allocate path and model
	Scenario<Number> path;
    allocatePath(prd.getDeflineRef(), path);
	cMdl->allocate(prd.getTimelineRef(), prd.getDeflineRef());

    //  Dimensions
    const size_t nPay = prd.getPayoffLabelsRef().size();
    const vector<Number*>& params = cMdl->getParameters();
    const size_t nParam = params.size();

    //  AAD - 1
    //  Access to tape
    Tape& tape = *Number::tape;
    //  Clear and initialise tape
    tape.clear();
	auto resetter = setNumResultsForAAD();
	//  Put parameters on tape
    //  note that also initializes all adjoints
    cMdl->putParametersOnTape();
    //  Init the simulation timeline
    //  CAREFUL: simulation timeline must be on tape
    //  Hence moved here
    cMdl->init(prd.getTimelineRef(), prd.getDeflineRef());
    //  Initialize path
    initializePath(path);
    //  Mark the tape straight after initialization
    tape.mark();
    //

    //  Init the RNG
    cRng->init(cMdl->simDim());                         
                                                            
    //  Allocate workspace
    vector<Number> nPayoffs(nPay);
    //  Gaussian vector
    vector<double> gaussVec(cMdl->simDim());            

    //  Results
    AADSimulResults results(nPath, nPay, nParam);

    //	Iterate through paths	
    for (size_t i = 0; i<nPath; i++)
    {
        //  AAD - 2
        //  Rewind tape to mark
        //  parameters stay on tape but the rest is wiped
        tape.rewindToMark();
        //

        //  Next Gaussian vector, dimension D
        cRng->nextG(gaussVec);
        //  Generate path, consume Gaussian vector
        cMdl->generatePath(gaussVec, path);     
        //	Compute result
        prd.computePayoffs(path, nPayoffs);
        //  Aggregate
        Number result = aggFun(nPayoffs);

        //  AAD - 3
        //  Propagate adjoints
        result.propagateToMark();
        //  Store results for the path
        results.aggregated[i] = double(result);
        convertCollection(
            nPayoffs.begin(), 
            nPayoffs.end(), 
            results.computePayoffs[i].begin());
		//
    }

    //  AAD - 4
    //  Mark = limit between pre-calculations and path-wise operations
    //  Operations above mark have been propagated and accumulated
    //  We conduct one propagation mark to start
    Number::propagateMarkToStart();
    //

    //  Pick sensitivities, summed over paths, and normalize
    transform(
        params.begin(),
        params.end(),
        results.risks.begin(),
        [nPath](const Number* p) {return p->adjoint() / nPath; });

    //  Clear the tape
    tape.clear();

    return results;
}

//  Parallel AAD, chapter 12

//  Init model and out on tape
inline void initModel4ParallelAAD(
    //  Inputs
    const Product<Number>&      prd,
    //  Cloned model, must have been allocated prior
    Model<Number>&              clonedMdl,
    //  Path, also allocated prior
    Scenario<Number>&           path)
{
    //  Access to tape
    Tape& tape = *Number::tape;
    //  Rewind tape
    tape.rewind();
    //  Put parameters on tape
    //  note that also initializes all adjoints
    clonedMdl.putParametersOnTape();
    //  Init the simulation timeline
    //  CAREFUL: simulation timeline must be on tape
    //  Hence moved here
    clonedMdl.init(prd.getTimelineRef(), prd.getDeflineRef());
    //  Path
    initializePath(path);
    //  Mark the tape straight after parameters
    tape.mark();
    //
}

//  Parallel version of mcSimulAAD()
template<class F = decltype(defaultAggregator)>
inline AADSimulResults
mcParallelSimulAAD(
    const Product<Number>&  prd,
    const Model<Number>&    mdl,
    const RNG& rng,
    const size_t            nPath,
    const F&                aggFun = defaultAggregator)
{
    const size_t nPay = prd.getPayoffLabelsRef().size();
    const size_t nParam = mdl.getNumParams();

    //  Allocate results
    AADSimulResults results(nPath, nPay, nParam);

    //  Clear and initialise tape
	Number::tape->clear();
	auto resetter = setNumResultsForAAD();
	
    //  We need one of all these for each thread
    //  0: main thread
    //  1 to n : worker threads

    ThreadPool *pool = ThreadPool::getInstance();
    const size_t nThread = pool->getTotalNumerOfThreads();

    //  Allocate workspace

    //  One model clone per thread
    vector<unique_ptr<Model<Number>>> models(nThread + 1);
    for (auto& model : models)
    {
        model = mdl.getClonePtr();
        model->allocate(prd.getTimelineRef(), prd.getDeflineRef());
    }

    //  One scenario per thread
    vector<Scenario<Number>> paths(nThread + 1);
    for (auto& path : paths)
    {
        allocatePath(prd.getDeflineRef(), path);
    }

    //  One vector of payoffs per thread
    vector<vector<Number>> computePayoffs(nThread + 1, vector<Number>(nPay));

    //  ~workspace

    //  Tapes for the worker threads
    //  The main thread has one of its own
    vector<Tape> tapes(nThread);

    //  Model initialized?
    //  Note we don't use vector<bool>
    //      because vector<bool> is not thread safe
    vector<int> mdlInit(nThread + 1, false);

    //  Initialize main thread
    initModel4ParallelAAD(prd, *models[0], paths[0]);

    //  Mark main thread as initialized
    mdlInit[0] = true;

    //  Init the RNGs, one pet thread
    //  One RNG per thread
    vector<unique_ptr<RNG>> rngs(nThread + 1);
    for (auto& random : rngs)
    {
        random = rng.getClonePtr();
        random->init(models[0]->simDim());
    }

    //  One Gaussian vector per thread
    vector<vector<double>> gaussVecs
        (nThread + 1, vector<double>(models[0]->simDim()));

    //  Reserve memory for futures
    vector<TaskHandle> futures;
    futures.reserve(nPath / BATCHSIZE + 1);

    //  Start
    //  Same as mcSimul() except we send tasks to the pool 
    //  instead of executing them

    size_t firstPath = 0;
    size_t pathsLeft = nPath;
    while (pathsLeft > 0)
    {
        size_t pathsInTask = min<size_t>(pathsLeft, BATCHSIZE);

        futures.push_back(pool->spawnTask([&, firstPath, pathsInTask]()
        {
            const size_t getThreadNumer = pool->getThreadNumer();

            //  Use this thread's tape
            //  Thread local magic: each thread its own pointer
            //  Note main thread = 0 is not reset
            if (getThreadNumer > 0) Number::tape = &tapes[getThreadNumer - 1];

            //  Initialize once on each thread
            if (!mdlInit[getThreadNumer])
            {
                //  Initialize
                initModel4ParallelAAD(prd, *models[getThreadNumer], paths[getThreadNumer]);

                //  Mark as initialized
                mdlInit[getThreadNumer] = true;
            }

            //  Get a RNG and position it correctly
            auto& random = rngs[getThreadNumer];
            random->skipTo(firstPath);

            //  And conduct the simulations, exactly same as sequential
            for (size_t i = 0; i < pathsInTask; i++)
            {
                //  Rewind tape to mark
                //  Notice : this is the tape for the executing thread

                Number::tape->rewindToMark();
                //  Next Gaussian vector, dimension D
                random->nextG(gaussVecs[getThreadNumer]);
                //  Path
                models[getThreadNumer]->generatePath(
                    gaussVecs[getThreadNumer], 
                    paths[getThreadNumer]);
                //  Payoff
                prd.computePayoffs(paths[getThreadNumer], computePayoffs[getThreadNumer]);

                //  Propagate adjoints
                Number result = aggFun(computePayoffs[getThreadNumer]);
                result.propagateToMark();
                //  Store results for the path
                results.aggregated[firstPath + i] = double(result);
                convertCollection(
                    computePayoffs[getThreadNumer].begin(), 
                    computePayoffs[getThreadNumer].end(),
                    results.computePayoffs[firstPath + i].begin());
            }

            //  Remember tasks must return bool
            return true;
        }));

        pathsLeft -= pathsInTask;
        firstPath += pathsInTask;
    }

    //  Wait and help
    for (auto& future : futures) pool->activeWait(future);
    
    //  Mark = limit between pre-calculations and path-wise operations
    //  Operations above mark have been propagated and accumulated
    //  We conduct one propagation mark to start
    //  On the main thread's tape
    Number::propagateMarkToStart();
    //  And on the worker thread's tapes
    Tape* mainThreadPtr = Number::tape;
    for (size_t i = 0; i < nThread; ++i)
    {
        if (mdlInit[i + 1])
        {
            //  Set tape pointer
            Number::tape = &tapes[i];
            //  On that tape, propagate
            Number::propagateMarkToStart();
        }
    }
    //  Reset tape to main thread's
    Number::tape = mainThreadPtr;

    //  Sum sensitivities over threads
    for (size_t j = 0; j < nParam; ++j)
    {
        results.risks[j] = 0.0;
        for (size_t i = 0; i < models.size(); ++i)
        {
            if (mdlInit[i]) results.risks[j] += models[i]->getParameters()[j]->adjoint();
        }
        results.risks[j] /= nPath;
    }

	//  Clear the main thread's tape
    //  The other tapes are cleared on the destruction of the vector of tapes
    Number::tape->clear();

    return results;
}

//  Multi-dimensional AAD, chapter 14
//	Rewrite code for the risk reports of multiple payoffs for clarity

struct AADMultiSimulResults
{
	AADMultiSimulResults(const size_t nPath, const size_t nPay, const size_t nParam) :
		computePayoffs(nPath, vector<double>(nPay)),
		risks(nParam, nPay)
	{}

	//  matrix(0..nPath - 1, 0..nPay - 1) of payoffs, same as mcSimul()
	vector<vector<double>>  computePayoffs;

	//  matrix(0..nParam - 1, 0..nPay - 1) of risk sensitivities
	//		of all payoffs, averaged over paths
	matrix<double>          risks;
};

//  Serial

inline AADMultiSimulResults
mcSimulAADMulti(
	const Product<Number>&  prd,
	const Model<Number>&    mdl,
	const RNG&              rng,
	const size_t            nPath)
{
	auto cMdl = mdl.getClonePtr();
	auto cRng = rng.getClonePtr();

	Scenario<Number> path;
	allocatePath(prd.getDeflineRef(), path);
    cMdl->allocate(prd.getTimelineRef(), prd.getDeflineRef());

	const size_t nPay = prd.getPayoffLabelsRef().size();
	const vector<Number*>& params = cMdl->getParameters();
	const size_t nParam = params.size();

	Tape& tape = *Number::tape;
	tape.clear();

    //  Set the AAD environment to multi-dimensional with dimension nPay
    //  Reset to 1D is automatic when resetter exits scope
	auto resetter = setNumResultsForAAD(true, nPay);

    cMdl->putParametersOnTape();
	cMdl->init(prd.getTimelineRef(), prd.getDeflineRef());
	initializePath(path);
	tape.mark();

	cRng->init(cMdl->simDim());

	vector<Number> nPayoffs(nPay);
	vector<double> gaussVec(cMdl->simDim());

    //  Allocate multi-dimensional results
    //      including a matrix(0..nParam - 1, 0..nPay - 1) of risk sensitivities
	AADMultiSimulResults results(nPath, nPay, nParam);

	for (size_t i = 0; i<nPath; i++)
	{
		tape.rewindToMark();

		cRng->nextG(gaussVec);
		cMdl->generatePath(gaussVec, path);
		prd.computePayoffs(path, nPayoffs);

        //  Multi-dimensional propagation
        //      client code seeds the tape with the correct boundary conditions 
		for (size_t j = 0; j < nPay; ++j)
		{
			nPayoffs[j].adjoint(j) = 1.0;
		}
        //      multi-dimensional propagation over simulation, end to mark
		Number::propagateAdjointsMulti(prev(tape.end()), tape.markIt());

		convertCollection(
            nPayoffs.begin(), 
            nPayoffs.end(), 
            results.computePayoffs[i].begin());
	}

    //  Multi-dimensional propagation over initialization, mark to start
	Number::propagateAdjointsMulti(tape.markIt(), tape.begin());

    //  Pack results 
	for (size_t i = 0; i < nParam; ++i)
	{
		for (size_t j = 0; j < nPay; ++j)
		{
			results.risks[i][j] = params[i]->adjoint(j) / nPath;
		}
	}

	tape.clear();

	return results;
}

//  Parallel

inline AADMultiSimulResults
mcParallelSimulAADMulti(
	const Product<Number>&  prd,
	const Model<Number>&    mdl,
	const RNG& rng,
	const size_t            nPath)
{
	const size_t nPay = prd.getPayoffLabelsRef().size();
	const size_t nParam = mdl.getNumParams();

	Number::tape->clear();
	auto resetter = setNumResultsForAAD(true, nPay);

	ThreadPool *pool = ThreadPool::getInstance();
	const size_t nThread = pool->getTotalNumerOfThreads();

	vector<unique_ptr<Model<Number>>> models(nThread + 1);
	for (auto& model : models)
	{
		model = mdl.getClonePtr();
		model->allocate(prd.getTimelineRef(), prd.getDeflineRef());
	}

	vector<Scenario<Number>> paths(nThread + 1);
	for (auto& path : paths)
	{
		allocatePath(prd.getDeflineRef(), path);
	}

	vector<vector<Number>> computePayoffs(nThread + 1, vector<Number>(nPay));

	vector<Tape> tapes(nThread);

	vector<int> mdlInit(nThread + 1, false);

	initModel4ParallelAAD(prd, *models[0], paths[0]);

	mdlInit[0] = true;

	vector<unique_ptr<RNG>> rngs(nThread + 1);
	for (auto& random : rngs)
	{
		random = rng.getClonePtr();
		random->init(models[0]->simDim());
	}

	vector<vector<double>> gaussVecs
	(nThread + 1, vector<double>(models[0]->simDim()));

	AADMultiSimulResults results(nPath, nPay, nParam);

	vector<TaskHandle> futures;
	futures.reserve(nPath / BATCHSIZE + 1);

	size_t firstPath = 0;
	size_t pathsLeft = nPath;
	while (pathsLeft > 0)
	{
		size_t pathsInTask = min<size_t>(pathsLeft, BATCHSIZE);

		futures.push_back(pool->spawnTask([&, firstPath, pathsInTask]()
		{
			const size_t getThreadNumer = pool->getThreadNumer();

			if (getThreadNumer > 0) Number::tape = &tapes[getThreadNumer - 1];

			if (!mdlInit[getThreadNumer])
			{
				initModel4ParallelAAD(prd, *models[getThreadNumer], paths[getThreadNumer]);
				mdlInit[getThreadNumer] = true;
			}

			auto& random = rngs[getThreadNumer];
			random->skipTo(firstPath);

			for (size_t i = 0; i < pathsInTask; i++)
			{

				Number::tape->rewindToMark();
				random->nextG(gaussVecs[getThreadNumer]);
				models[getThreadNumer]->generatePath(
					gaussVecs[getThreadNumer],
					paths[getThreadNumer]);
				prd.computePayoffs(paths[getThreadNumer], computePayoffs[getThreadNumer]);

				const size_t n = computePayoffs[getThreadNumer].size();
				for (size_t j = 0; j < n; ++j)
				{
					computePayoffs[getThreadNumer][j].adjoint(j) = 1.0;
				}
				Number::propagateAdjointsMulti(prev(Number::tape->end()), Number::tape->markIt());

				convertCollection(
					computePayoffs[getThreadNumer].begin(),
					computePayoffs[getThreadNumer].end(),
					results.computePayoffs[firstPath + i].begin());
			}

			return true;
		}));

		pathsLeft -= pathsInTask;
		firstPath += pathsInTask;
	}

	for (auto& future : futures) pool->activeWait(future);

	Number::propagateAdjointsMulti(Number::tape->markIt(), Number::tape->begin());
	for (size_t i = 0; i < nThread; ++i)
	{
		if (mdlInit[i + 1])
		{
			Number::propagateAdjointsMulti(tapes[i].markIt(), tapes[i].begin());
		}
	}

	for (size_t j = 0; j < nParam; ++j) for (size_t k = 0; k < nPay; ++k)
	{
		results.risks[j][k] = 0.0;
		for (size_t i = 0; i < models.size(); ++i)
		{
			if (mdlInit[i]) results.risks[j][k] += models[i]->getParameters()[j]->adjoint(k);
		}
		results.risks[j][k] /= nPath;
	}

	Number::tape->clear();

	return results;
}
