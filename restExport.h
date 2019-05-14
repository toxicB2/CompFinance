#pragma once
// restExport.h - Contains endpoints for REST implementation
#pragma once

#ifdef RESTEXPORTS_EXPORTS
#define RESTEXPORTS_API __declspec(dllexport)
#else
#define RESTEXPORTS_API __declspec(dllimport)
#endif

NumericalParam rest2num(const double              useSobol,
	                    const double              seed1,
	                    const double              seed2,
	                    const double              numPath,
	                    const double              parallel);

extern "C" RESTEXPORTS_API double restRestartThreadPool(double xNthread);

extern "C" RESTEXPORTS_API string restPutDupire(double              spot,
	                                            vector<double>      spots,
	                                            vector<double>      times,
	                                            matrix<double>      vols,
	                                            double              maxDt,
	                                            string              xid);

extern "C" RESTEXPORTS_API string restPutBlackScholes(double              spot,
	                                                  double              vol,
	                                                  double              qSpot,
	                                                  double              rate,
	                                                  double              div,
	                                                  string              xid);

extern "C" RESTEXPORTS_API string restPutEuropean(double              strike,
	                                              double              exerciseDate,
	                                              double              settlementDate,
	                                              string              xid);

extern "C" RESTEXPORTS_API string restPutBarrier(double              strike,
	                                             double              barrier,
	                                             double              maturity,
	                                             double              monitorFreq,
	                                             double              smoothing,
	                                             string              xid);

extern "C" RESTEXPORTS_API string restPutContingent(double              coupon,
	                                                double              maturity,
	                                                double              payFreq,
	                                                double              smoothing,
	                                                string          xid);

extern "C" RESTEXPORTS_API string restPutEuropeans(vector<double>               vmaturities,
	                                               vector<double>               vstrikes,
	                                               string                       xid);

extern "C" RESTEXPORTS_API string* restPayoffIds(string  xid);

extern "C" RESTEXPORTS_API void restParameters(string  xid, string* paramLabelsOut, double* paramsTimesOut);

extern "C" RESTEXPORTS_API void restValue(string     modelid,
	                                      string     productid,
	                                      double     useSobol,
	                                      double     seed1,
	                                      double     seed2,
	                                      double     numPath,
	                                      double     parallel,
	                                      string* paramLabelsOut,
	                                      double* paramsTimesOut);

extern "C" RESTEXPORTS_API void restValueTime(string              modelid,
	                                          string              productid,
	                                          double              useSobol,
	                                          double              seed1,
	                                          double              seed2,
	                                          double              numPath,
	                                          double              parallel,

	                                          string*             paramLabelsOut,
	                                          double*             paramsTimesOut,
	                                          double              runningTime);
