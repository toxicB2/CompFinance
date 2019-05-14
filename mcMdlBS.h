#pragma once

template <class T>
class BlackScholes : public Model<T>
{
    T                   spotPrice;
    T                   currentRate;
    T                   constDividend;
    T                   flatVolatility;

    const bool          isSpotMeasure; //false risk neutral measure
    vector<Time>        timeline4Security; //Similuation timeline = today + event dates
    bool                isTodayOnTimeline;  
    //  The pruduct's defline byref
    const vector<SampleStructure>*    defline4Security;

    //  Pre-calculated on initialization
    vector<T>           stDeviations;
    vector<T>           bsDrifts;
    vector<T>           bsNumeraires;
    vector<vector<T>>   bsDiscounts;
    vector<vector<T>>   bsForwardFactors;
    vector<vector<T>>   bsLibors;//  and rates = (exp(r * (T2 - T1)) - 1) / (T2 - T1)

    vector<T*>          exportParametersPoiters;
    vector<string>      exportParameterLabels;

public:
    template <class U>
    BlackScholes(
        const U             spot,
        const U             vol,
        const bool          spotMeasure = false,
        const U             rate = U(0.0),
        const U             div = U(0.0)) : 
            spotPrice(spot),
            flatVolatility(vol),
            currentRate(rate),
            constDividend(div),
            isSpotMeasure(spotMeasure),
            exportParametersPoiters(4),
            exportParameterLabels(4) 
    {
        exportParameterLabels[0] = "spot";
        exportParameterLabels[1] = "vol";
        exportParameterLabels[2] = "rate";
        exportParameterLabels[3] = "div";

        setParamPointers();
    }

private:
    void setParamPointers() //  Has to be reset on copy
    {
        exportParametersPoiters[0] = &spotPrice;
        exportParametersPoiters[1] = &flatVolatility;
        exportParametersPoiters[2] = &currentRate;
		exportParametersPoiters[3] = &constDividend;
    }

public:

    T getSpot() const{ return spotPrice; }
    const T getVol() const { return flatVolatility; }
    const T getRate() const { return currentRate; }
    const T getDividend() const { return constDividend; }
    const vector<T*>& getParameters() override { return exportParametersPoiters;}
    const vector<string>& getParameterLabels() const override { return exportParameterLabels; }

    unique_ptr<Model<T>> getClonePtr() const override
    {
        auto clonePtr = make_unique<BlackScholes<T>>(*this);
		clonePtr->setParamPointers();
        return clonePtr;
    }

    void allocate( const vector<Time>&         productTimeline, 
                   const vector<SampleStructure>&    deflineRef) override
    {
        //  Simulation timeline = today + product timeline
        timeline4Security.clear();
        timeline4Security.push_back(systemTime);
        for (const auto& time : productTimeline){ if (time > systemTime) timeline4Security.push_back(time);}

        isTodayOnTimeline = (productTimeline[0] == systemTime);

        defline4Security = &deflineRef;

        stDeviations.resize(timeline4Security.size() - 1);
        bsDrifts.resize(timeline4Security.size() - 1);

        const size_t n = productTimeline.size();
        bsNumeraires.resize(n);       
        bsDiscounts.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            bsDiscounts[j].resize(deflineRef[j].discountMats.size());
        }

        bsForwardFactors.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            bsForwardFactors[j].resize(deflineRef[j].forwardMats.size());
        }

        bsLibors.resize(n);
        for (size_t j = 0; j < n; ++j)
        {
            bsLibors[j].resize(deflineRef[j].liborDefs.size());
        }
    }

    void init(
        const vector<Time>&         productTimeline, 
        const vector<SampleStructure>&    getDeflineRef) override
    {        
        const T mu = currentRate - constDividend;
        const size_t n = timeline4Security.size() - 1;

        for (size_t i = 0; i < n; ++i)
        {
            const double dt = timeline4Security[i + 1] - timeline4Security[i];
            stDeviations[i] = flatVolatility * sqrt(dt);
            
            if (isSpotMeasure)
            {
                //lognormal mean
                bsDrifts[i] = (mu + 0.5*flatVolatility*flatVolatility)*dt;
            }
            else
            {
                //geometric BM mean
                bsDrifts[i] = (mu - 0.5*flatVolatility*flatVolatility)*dt;
            }
        }

        const size_t m = productTimeline.size();

		for (size_t i = 0; i < m; ++i) //  loop on event dates
		{
			if (getDeflineRef[i].needNumeraire)
			{
				if (isSpotMeasure)
				{
                    //      num(t) = spot(t) / spot(0) * exp(div * t)
                    //      precalculate exp(div * t) / spot(0)
                    bsNumeraires[i] = exp(constDividend * productTimeline[i]) / spotPrice;
				}
				else
				{
                    //      numeraire is deterministic in Black-Scholes = exp(rate * t)
                    bsNumeraires[i] = exp(currentRate * productTimeline[i]);
				}
			}

			//  Discount factors
			const size_t pDF = getDeflineRef[i].discountMats.size();
			for (size_t j = 0; j < pDF; ++j)
			{
				bsDiscounts[i][j] =
					exp(-currentRate * (getDeflineRef[i].discountMats[j] - productTimeline[i]));
			}

			//  Forward factors
			const size_t pFF = getDeflineRef[i].forwardMats.size();
			for (size_t j = 0; j < pFF; ++j)
			{
				bsForwardFactors[i][j] =
					exp(mu * (getDeflineRef[i].forwardMats[j] - productTimeline[i]));
			}

			//  Libors
			const size_t pL = getDeflineRef[i].liborDefs.size();
			for (size_t j = 0; j < pL; ++j)
			{
				const double dt
					= getDeflineRef[i].liborDefs[j].end - getDeflineRef[i].liborDefs[j].start;
				bsLibors[i][j] = (exp(currentRate*dt) - 1.0) / dt;
			}
		}   
	}

	size_t getSimulationDimension() const override { return timeline4Security.size() - 1; }

private:
    inline void fillScen( const size_t                idx,    
                          const T&                    spot,  
                          ValuesForSampleStruct<T>&   scen,   
                          const SampleStructure&      def) const
    {
        if (def.needNumeraire)
        {
			scen.numeraire = bsNumeraires[idx];
            if (isSpotMeasure) scen.numeraire *= spot;
        }
        
        transform(bsForwardFactors[idx].begin(), bsForwardFactors[idx].end(), scen.forwards.begin(), 
			[&spot](const T& ff) { return spot * ff;}
                  );

        copy(bsDiscounts[idx].begin(), bsDiscounts[idx].end(), scen.discounts.begin());
        copy(bsLibors[idx].begin(), bsLibors[idx].end(), scen.libors.begin());
    }

public:

    void generatePath( const vector<double>&   gaussVec, 
                       Scenario<T>&            path) const override
    {
        T spot = spotPrice; // on tape
        size_t idx = 0;
        if (isTodayOnTimeline)
        {
            fillScen(idx, spot, path[idx], (*defline4Security)[idx]);
            ++idx;
        }

        const size_t n = timeline4Security.size() - 1;
        for (size_t i = 0; i < n; ++i)
        {
            spot = spot * exp(bsDrifts[i] + stDeviations[i] * gaussVec[i]);
            fillScen(idx, spot, path[idx], (*defline4Security)[idx]);
            ++idx;
        }
    }
};