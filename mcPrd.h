#pragma once

#include <map>
#include "mcBase.h"

#define ONE_HOUR 0.000114469
#define ONE_DAY 0.003773585

template <class T>
class European : public Product<T>
{
    double              optionStrike;
    Time                optionExerciseDate;
    Time                optionSettlementDate;

    vector<Time>              timeline4Security;
    vector<SampleStructure>   defline4Security;

    vector<string>      securityLabels;

public:
    European(const double   strike, 
             const Time     exerciseDate,
             const Time     settlementDate) :
                            optionStrike(strike),
                            optionExerciseDate(exerciseDate),
                            optionSettlementDate(settlementDate),
                            securityLabels(1)
    {

        timeline4Security.push_back(exerciseDate);
        defline4Security.resize(1); 

        defline4Security[0].needNumeraire = true;

        defline4Security[0].forwardMats.push_back(settlementDate);
        defline4Security[0].discountMats.push_back(settlementDate);

        ostringstream ost;
        ost.precision(2);
        ost << fixed;
        if (settlementDate == exerciseDate)
        {
            ost << "call " << optionStrike << " " << exerciseDate;
        }
        else
        {
            ost << "call " << optionStrike << " " << exerciseDate << " " << settlementDate;
        }
        securityLabels[0] = ost.str();
    }

    European(const double   strike,
             const Time     exerciseDate) : 
             European(strike, exerciseDate, exerciseDate)
    {}

    unique_ptr<Product<T>> getClonePtr() const override { return make_unique<European<T>>(*this);}
    const vector<Time>& getTimelineRef() const override { return timeline4Security; }
    const vector<SampleStructure>& getDeflineRef() const override { return defline4Security; }
    const vector<string>& getPayoffLabelsRef() const override { return securityLabels; }

    void computePayoffs( const Scenario<T>&          path,
                         vector<T>&                  payoffs) const override
    {
        payoffs[0] = max(path[0].forwards[0] - optionStrike, 0.0)* path[0].discounts[0] / path[0].needNumeraire; 
    }
};

template <class T>
class UOC : public Product<T>
{
    double              optionStrike;
    double              securityBarrier;
    Time                maturityOfSecurity;

    double              smoothingValue4Barrier;

    vector<Time>        timeline4Security;
    vector<SampleStructure>   defline4Security;

    vector<string>      securityLabels;

public:
    UOC(const double    strike, 
        const double    barrier, 
        const Time      maturity, 
        const Time      monitorFreq,
        const double    smooth)
        : optionStrike(strike), 
          securityBarrier(barrier), 
          maturityOfSecurity(maturity),
          smoothingValue4Barrier(smooth),
          securityLabels(2)
    {
        timeline4Security.push_back(systemTime); //start from today and jump through monitoring frequency
        Time t = systemTime + monitorFreq;
            
        while (maturityOfSecurity - t > ONE_HOUR)
        {
            timeline4Security.push_back(t);
            t += monitorFreq;
        }
        timeline4Security.push_back(maturityOfSecurity);

        const size_t n = timeline4Security.size();
        defline4Security.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            defline4Security[i].needNumeraire = false; //  Numeraire is needed only on last step
            defline4Security[i].forwardMats.push_back(timeline4Security[i]); //  spot(t) = forward (t, t) 
        }
        defline4Security.back().needNumeraire = true;

        ostringstream ost;
        ost.precision(2);
        ost << fixed;
        ost << "call " << maturityOfSecurity << " " << optionStrike ;
        securityLabels[1] = ost.str();

        ost << " up and out "<< securityBarrier << " monitoring freq " << monitorFreq << " smooth " << smoothingValue4Barrier;
        securityLabels[0] = ost.str();
    }

    unique_ptr<Product<T>> getClonePtr() const override { return make_unique<UOC<T>>(*this); }
    const vector<Time>& getTimelineRef() const override { return timeline4Security; }
    const vector<SampleStructure>& getDeflineRef() const override { return defline4Security; }
    const vector<string>& getPayoffLabelsRef() const override { return securityLabels; }

    void computePayoffs( const Scenario<T>&          path,
                         vector<T>&                  payoffs) const override
    {
        //  See Savine's presentation on Fuzzy Logic, Global Derivatives 2016
        const double smooth = double(path[0].forwards[0] * smoothingValue4Barrier),
            twoSmooth = 2 * smooth,
            barSmooth = securityBarrier + smooth;

        T alive(1.0); // Start alive

        for (const auto& sample: path)
        {
            if (sample.forwards[0] > barSmooth) //  Breached
            {
                alive = T(0.0);
                break;
            }         
            if (sample.forwards[0] > securityBarrier - smooth) //  Semi-breached: apply smoothing
            {
                alive *= (barSmooth - sample.forwards[0]) / twoSmooth;
            }
        }

       payoffs[1] = max(path.back().forwards[0] - optionStrike, 0.0) / path.back().needNumeraire;
       payoffs[0] = alive * payoffs[1];
    }
};

template <class T>
class Europeans : public Product<T>
{
	//	Timeline
	vector<Time>            myMaturities;
	//  One vector of strikes per maturity
	vector<vector<double>>  myStrikes;
	vector<SampleStructure>       defline4Security;

	vector<string>          securityLabels;

public:

	//  Constructor: store data and build timeline
	Europeans(const map<Time, vector<double>>& options)
	{
		const size_t n = options.size();

		//  Timeline = each maturity is an event date
		for (const pair<Time, vector<double>>& p : options)
		{
			myMaturities.push_back(p.first);
			myStrikes.push_back(p.second);
		}

		//  Defline = num and spot(t) = forward(t,t) on every step
		defline4Security.resize(n);
		for (size_t i = 0; i < n; ++i)
		{
			defline4Security[i].needNumeraire = true;
			defline4Security[i].forwardMats.push_back(myMaturities[i]);
		}

		//  Identify the payoffs
		for (const auto& option : options)
		{
			for (const auto& strike : option.second)
			{
				ostringstream ost;
				ost.precision(2);
				ost << fixed;
				ost << "call " << option.first << " " << strike;
				securityLabels.push_back(ost.str());
			}
		}
	}

	//  access to maturities and strikes
	const vector<Time>& maturities() const
	{
		return myMaturities;
	}

	const vector<vector<double>>& strikes() const
	{
		return myStrikes;
	}

	//  Virtual copy constructor
	unique_ptr<Product<T>> getClonePtr() const override
	{
		return make_unique<Europeans<T>>(*this);
	}

	//  Timeline
	const vector<Time>& getTimelineRef() const override
	{
		return myMaturities;
	}

	//  Defline
	const vector<SampleStructure>& getDeflineRef() const override
	{
		return defline4Security;
	}

	//  Labels
	const vector<string>& getPayoffLabelsRef() const override
	{
		return securityLabels;
	}

	//  Payoffs, maturity major
	void computePayoffs(
		//  path, one entry per time step 
		const Scenario<T>&          path,
		//  pre-allocated space for resulting payoffs
		vector<T>&                  computePayoffs)
		const override
	{
		const size_t numT = myMaturities.size();

		auto payoffIt = computePayoffs.begin();
		for (size_t i = 0; i < numT; ++i)
		{
			transform(
				myStrikes[i].begin(),
				myStrikes[i].end(),
				payoffIt,
				[spot = path[i].forwards[0], num = path[i].needNumeraire]
				(const double& k)
				{
					return max(spot - k, 0.0) / num;
				}
			);

			payoffIt += myStrikes[i].size();
		}
	}
};

//  Payoff = sum { (libor(Ti, Ti+1) + cpn) 
//      * coverage(Ti, Ti+1) only if Si+1 >= Si }
template <class T>
class ContingentBond : public Product<T>
{
    Time                maturityOfSecurity;
    double              myCpn;
    double              smoothingValue4Barrier;

    vector<Time>        timeline4Security;
    vector<SampleStructure>   defline4Security;

    vector<string>      securityLabels;

    //  Pre-computed coverages
    vector<double>      myDt;

public:

    //  Constructor: store data and build timeline
    //  Timeline = system date to maturity, 
    //  with steps every payment frequency
    ContingentBond(
        const Time      maturity,
        const double    cpn,
        const Time      payFreq,
        const double    smooth)
        : 
        maturityOfSecurity(maturity),
        myCpn(cpn),
        smoothingValue4Barrier(smooth),
        securityLabels(1)
    {
        //  Timeline

        //  Today
        timeline4Security.push_back(systemTime);
        Time t = systemTime + payFreq;

        //  Payment schedule
        while (maturityOfSecurity - t > ONE_DAY)
        {
            myDt.push_back(t - timeline4Security.back());
            timeline4Security.push_back(t);
            t += payFreq;
        }

        //  Maturity
        myDt.push_back(maturityOfSecurity - timeline4Security.back());
        timeline4Security.push_back(maturityOfSecurity);

        //

        //  Defline

        //  Payoff = sum { (libor(Ti, Ti+1) + cpn) 
        //      * coverage(Ti, Ti+1) only if Si+1 >= Si }
        //  We need spot ( = forward (Ti, Ti) ) on every step,
        //  and libor (Ti, Ti+1) on on every step but the last
        //  (coverage is assumed act/365)

        const size_t n = timeline4Security.size();
        defline4Security.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            //  spot(Ti) = forward (Ti, Ti) needed on every step
            defline4Security[i].forwardMats.push_back(timeline4Security[i]);

            //  libor(Ti, Ti+1) and discount (Ti, Ti+1) needed 
            //      on every step but last
            if (i < n - 1)
            {
                defline4Security[i].liborDefs.push_back(
                    SampleStructure::RateDef(timeline4Security[i], timeline4Security[i + 1], "libor"));
            }

            //  Numeraire needed only on every step but first
            defline4Security[i].needNumeraire = i > 0;
        }

        //  Identify the product
        ostringstream ost;
        ost.precision(2);
        ost << fixed;
        ost << "contingent bond " << maturityOfSecurity << " " << myCpn;
        securityLabels[0] = ost.str();
    }

    //  Virtual copy constructor
    unique_ptr<Product<T>> getClonePtr() const override
    {
        return make_unique<ContingentBond<T>>(*this);
    }

    //  Timeline
    const vector<Time>& getTimelineRef() const override
    {
        return timeline4Security;
    }

    //  Defline
    const vector<SampleStructure>& getDeflineRef() const override
    {
        return defline4Security;
    }

    //  Labels
    const vector<string>& getPayoffLabelsRef() const override
    {
        return securityLabels;
    }

    //  Payoff
    void computePayoffs(
        //  path, one entry per time step 
        const Scenario<T>&          path,
        //  pre-allocated space for resulting payoffs
        vector<T>&                  computePayoffs)
        const override
    {
        //  We apply the smooth digital technique to stabilize risks
        //  See Savine's presentation on Fuzzy Logic, Global Derivatives 2016

        //  We apply a smoothing factor of x% of the spot both ways
        //  untemplated
        const double smooth = double(path[0].forwards[0] * smoothingValue4Barrier),
            twoSmooth = 2 * smooth;

        //  Period by period
        const size_t n = path.size() - 1;
        computePayoffs[0] = 0;
        for (size_t i = 0; i < n; ++i)
        {
            const auto& start = path[i];
            const auto& end = path[i + 1];

            const T s0 = start.forwards[0], s1 = end.forwards[0];

            //  Is asset performance positive?

            /*
                We apply smoothing here otherwise risks are unstable with bumps
                    and wrong with AAD 

                bool digital = end.forwards[0] >= start.forwards[0];
            */

            T digital;
            if (s1 - s0 > smooth)
            {
                digital = T(1.0);
            }
            else if (s1 - s0 < - smooth)
            {
                digital = T(0.0);
            }
            else // "fuzzy" barrier = interpolate
            {
                digital = (s1 - s0 + smooth) / twoSmooth;
            }

            //  ~smoothing

            computePayoffs[0] += 
                digital             //  contingency
                * ( start.libors[0] //  libor(Ti, Ti+1)
                + myCpn)            //  + coupon
                * myDt[i]           //  day count / 365
                / end.needNumeraire;    //  paid at Ti+1
        }
        computePayoffs[0] += 1.0 / path.back().needNumeraire;  //  redemption at maturity
    }
};
